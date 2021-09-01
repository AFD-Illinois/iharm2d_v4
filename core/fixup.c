/*---------------------------------------------------------------------------------

  FIXUP.C

  -Repair integration failures
  -Apply geometrical floors to density and internal energy
  -Apply ceiling to Lorentz factor
  -Apply ceiling to total fluid entropy (only if ELECTRONS enabled)
  -Apply ceiling to sigma and plasma beta inverse
  -Apply temperature ceiling
  -Replace inversion failure grid zones with values interpolated from 
   neighbouring zones

---------------------------------------------------------------------------------*/

#include "decs.h"

// Floor Codes: bit masks
#define HIT_FLOOR_GEOM_RHO 1
#define HIT_FLOOR_GEOM_U 2
#define HIT_FLOOR_B_RHO 4
#define HIT_FLOOR_B_U 8
#define HIT_FLOOR_TEMP 16
#define HIT_FLOOR_GAMMA 32
#define HIT_FLOOR_KTOT 64

// Point in m, around which to steepen floor prescription, eventually toward r^-3
#define FLOOR_R_CHAR 10

static struct FluidState *Stmp;

void fixup_ceiling(struct GridGeom *G, struct FluidState *S, int i, int j);
void fixup_floor(struct GridGeom *G, struct FluidState *S, int i, int j);

// Apply floors to density, internal energy
void fixup(struct GridGeom *G, struct FluidState *S)
{
  timer_start(TIMER_FIXUP);

  static int firstc = 1;
  if (firstc)
  {
    Stmp = calloc(1,sizeof(struct FluidState));
    firstc = 0;
  }

#pragma omp parallel for simd
  ZLOOPALL fflag[j][i] = 0;

#pragma omp parallel for collapse(2)
  ZLOOP fixup_ceiling(G, S, i, j);

  // Bulk call before bsq calculation below
  get_state_vec(G, S, CENT, 0, N2-1, 0, N1-1);

#pragma omp parallel for collapse(2)
  ZLOOP fixup_floor(G, S, i, j);

  // Some debug info about floors
#if DEBUG
  int n_geom_rho = 0, n_geom_u = 0, n_b_rho = 0, n_b_u = 0, n_temp = 0, n_gamma = 0, n_ktot = 0;

#pragma omp parallel for collapse(2) reduction(+:n_geom_rho) reduction(+:n_geom_u) \
    reduction(+:n_b_rho) reduction(+:n_b_u) reduction(+:n_temp) reduction(+:n_gamma) reduction(+:n_ktot)
  ZLOOP {
    int flag = fflag[j][i];
    if (flag & HIT_FLOOR_GEOM_RHO) n_geom_rho++;
    if (flag & HIT_FLOOR_GEOM_U) n_geom_u++;
    if (flag & HIT_FLOOR_B_RHO) n_b_rho++;
    if (flag & HIT_FLOOR_B_U) n_b_u++;
    if (flag & HIT_FLOOR_TEMP) n_temp++;
    if (flag & HIT_FLOOR_GAMMA) n_gamma++;
    if (flag & HIT_FLOOR_KTOT) n_ktot++;
  }

  LOG("FLOORS:");
  if (n_geom_rho > 0) LOGN("Hit %d GEOM_RHO", n_geom_rho);
  if (n_geom_u > 0) LOGN("Hit %d GEOM_U", n_geom_u);
  if (n_b_rho > 0) LOGN("Hit %d B_RHO", n_b_rho);
  if (n_b_u > 0) LOGN("Hit %d B_U", n_b_u);
  if (n_temp > 0) LOGN("Hit %d TEMPERATURE", n_temp);
  if (n_gamma > 0) LOGN("Hit %d GAMMA", n_gamma);
  if (n_ktot > 0) LOGN("Hit %d KTOT", n_ktot);

#endif

  LOG("End fixup");

  timer_stop(TIMER_FIXUP);
}

inline void fixup_ceiling(struct GridGeom *G, struct FluidState *S, int i, int j)
{
  // First apply ceilings:
  // 1. Limit gamma with respect to normal observer
  double gamma = mhd_gamma_calc(G, S, i, j, CENT);

  if (gamma > GAMMAMAX)
  {
    fflag[j][i] |= HIT_FLOOR_GAMMA;
    double f = sqrt((GAMMAMAX*GAMMAMAX - 1.)/(gamma*gamma - 1.));
    S->P[U1][j][i] *= f;
    S->P[U2][j][i] *= f;
    S->P[U3][j][i] *= f;
  }

  // 2. Limit KTOT
#if ELECTRONS
    // Keep to KTOTMAX by controlling u, to avoid anomalous cooling from funnel wall
    // Note: This operates on last iteration's KTOT, meaning the effective value can escape the ceiling
    if (S->P[KTOT][j][i] > KTOTMAX)
    {
      fflag[j][i] |= HIT_FLOOR_KTOT;
      S->P[UU][j][i] = KTOTMAX*pow(S->P[RHO][j][i],gam)/(gam-1.);
      S->P[KTOT][j][i] = KTOTMAX;
    }
#endif
}

inline void fixup_floor(struct GridGeom *G, struct FluidState *S, int i, int j)
{
  // Then apply floors:
  // 1. Geometric hard floors, not based on fluid relationships
  double rhoflr_geom, uflr_geom;
  if(METRIC == MKS)
  {
    double r, th, X[NDIM];
    coord(i, j, CENT, X);
    bl_coord(X, &r, &th);

    // New, steeper floor in rho
    // Previously raw r^-2, r^-1.5
    double rhoscal = pow(r, -2.) * 1 / (1 + r/FLOOR_R_CHAR);
    rhoflr_geom = RHOMIN*rhoscal;
    uflr_geom = UUMIN*pow(rhoscal, gam);

    // Impose overall minimum
    // TODO These would only be hit at by r^-3 floors for r_out = 100,000M.  Worth keeping?
    rhoflr_geom = MY_MAX(rhoflr_geom, RHOMINLIMIT);
    uflr_geom = MY_MAX(uflr_geom, UUMINLIMIT);
  } 
  else if (METRIC == MINKOWSKI)
  {
    rhoflr_geom = RHOMIN*1.e-2;
    uflr_geom = UUMIN*1.e-2;
  }

  // Record Geometric floor hits
  if (rhoflr_geom > S->P[RHO][j][i]) fflag[j][i] |= HIT_FLOOR_GEOM_RHO;
  if (uflr_geom > S->P[UU][j][i]) fflag[j][i] |= HIT_FLOOR_GEOM_U;


  // 2. Magnetic floors: impose maximum magnetization sigma = bsq/rho, inverse beta prop. to bsq/U
  double bsq = bsq_calc(S, i, j);
  double rhoflr_b = bsq/BSQORHOMAX;
  double uflr_b = bsq/BSQOUMAX;

  // Record Magnetic floor hits
  if (rhoflr_b > S->P[RHO][j][i]) fflag[j][i] |= HIT_FLOOR_B_RHO;
  if (uflr_b > S->P[UU][j][i]) fflag[j][i] |= HIT_FLOOR_B_U;

  // Evaluate highest U floor
  double uflr_max = MY_MAX(uflr_geom, uflr_b);

  // 3. Temperature ceiling: impose maximum temperature
  // Take floors on U into account
  double rhoflr_temp = MY_MAX(S->P[UU][j][i] / UORHOMAX, uflr_max / UORHOMAX);

  // Record hitting temperature ceiling
  if (rhoflr_temp > S->P[RHO][j][i]) fflag[j][i] |= HIT_FLOOR_TEMP; // Misnomer for consistency

  // Evaluate highest RHO floor
  double rhoflr_max = MY_MAX(MY_MAX(rhoflr_geom, rhoflr_b), rhoflr_temp);

  if (rhoflr_max > S->P[RHO][j][i] || uflr_max > S->P[UU][j][i]) 
  { // Apply floors

    // Initialize a dummy fluid parcel
    PLOOP
    {
      Stmp->P[ip][j][i] = 0;
      Stmp->U[ip][j][i] = 0;
    }

    // Add mass and internal energy, but not velocity
    Stmp->P[RHO][j][i] = MY_MAX(0., rhoflr_max - S->P[RHO][j][i]);
    Stmp->P[UU][j][i] = MY_MAX(0., uflr_max - S->P[UU][j][i]);

    // Get conserved variables for the parcel
    get_state(G, Stmp, i, j, CENT);
    prim_to_flux(G, Stmp, i, j, 0, CENT, Stmp->U);

    // And for the current state
    prim_to_flux(G, S, i, j, 0, CENT, S->U);

    // Add new conserved variables to current values
    PLOOP
    {
      S->U[ip][j][i] += Stmp->U[ip][j][i];
      S->P[ip][j][i] += Stmp->P[ip][j][i];
    }

    // Recover primitive variables
    pflag[j][i] = U_to_P(G, S, i, j, CENT);
  }

#if ELECTRONS
    // Reset entropy after floors
    S->P[KTOT][j][i] = (gam - 1.)*S->P[UU][j][i]/pow(S->P[RHO][j][i],gam);
#endif

}

// Replace bad points with values interpolated from neighbors
#define FLOOP for(int ip=0;ip<B1;ip++)
void fixup_utoprim(struct GridGeom *G, struct FluidState *S)
{
  timer_start(TIMER_FIXUP);

  // Flip the logic of the pflag[] so that it now indicates which cells are good
#pragma omp parallel for simd collapse(2)
  ZLOOPALL
    pflag[j][i] = !pflag[j][i];

#if DEBUG
  int nbad_utop = 0;
#pragma omp parallel for simd collapse(2) reduction (+:nbad_utop)
  ZLOOP 
  {
    // Count the 0 = bad cells
    nbad_utop += !pflag[j][i];
  }
  LOGN("Fixing %d bad cells", nbad_utop);
#endif

  // Make sure we are not using ill defined physical corner regions
  // TODO find a way to do this once, or put it in bounds at least?
  for (int j = 0; j < NG; j++)
  {
    for (int i = 0; i < NG; i++)
    {
      pflag[j][i] = 0;
      pflag[j][i+N1+NG] = 0;
      pflag[j+N2+NG][i] = 0;
      pflag[j+N2+NG][i+N1+NG] = 0;
    }
  }

#if DEBUG
  // Keep track of how many points we fix
  int nfixed_utop = 0;
#endif

  ZLOOP 
  {
    if (pflag[j][i] == 0)
    {
      double wsum = 0.;
      double sum[B1];
      FLOOP sum[ip] = 0.;
      for (int l = -1; l < 2; l++)
      {
        for (int m = -1; m < 2; m++) 
        {
          double w = 1./(abs(l) + abs(m) + 1)*pflag[j+m][i+l];
          wsum += w;
          FLOOP sum[ip] += w*S->P[ip][j+m][i+l];
        }
      }
      if(wsum < 1.e-10)
      {
        fprintf(stderr, "fixup_utoprim: No usable neighbors at %d %d\n", i, j);
        continue;
      }
      FLOOP S->P[ip][j][i] = sum[ip]/wsum;

#if DEBUG
      nfixed_utop++;
#endif

      // Make sure fixed values still abide by floors
      fixup_ceiling(G, S, i, j);
      get_state(G, S, i, j, CENT);
      fixup_floor(G, S, i, j);
    }
  }

#if DEBUG
  int nleft_utop = nbad_utop - nfixed_utop;
  if(nleft_utop > 0) fprintf(stderr,"Cells STILL BAD after fixup_utoprim: %d\n", nleft_utop);
#endif

  // Reset the pflag
#pragma omp parallel for simd
  ZLOOPALL
    pflag[j][i] = 0;

  timer_stop(TIMER_FIXUP);
}
#undef FLOOP
