/*---------------------------------------------------------------------------------

  PHYS.C

  -Compute MHD stress-energy tensor
  -Obtain conserved variables and/or fluxes from primitives
  -Compute Lorentz factor
  -Compute 4 velocities from primitives
  -Calculate components of magnetosonic velocity
  -Calculate source terms
  -Compute bsq at given grid zone

---------------------------------------------------------------------------------*/

#include "decs.h"

// MHD stress-energy tensor with first index up, second index down. A factor of sqrt(4 pi) is absorbed into the definition of b.
inline void mhd_calc(struct FluidState *S, int i, int j, int dir, double *mhd)
{
  double u, pres, w, bsq, eta, ptot;

  u = S->P[UU][j][i];
  pres = (gam - 1.)*u;
  w = pres + S->P[RHO][j][i] + u;
  bsq = bsq_calc(S, i, j);
  eta = w + bsq;
  ptot = pres + 0.5*bsq;

  DLOOP1
    mhd[mu] = eta*S->ucon[dir][j][i]*S->ucov[mu][j][i] + ptot*delta(dir, mu) - S->bcon[dir][j][i]*S->bcov[mu][j][i];
}

void prim_to_flux(struct GridGeom *G, struct FluidState *S, int i, int j, int dir, int loc, GridPrim flux)
{
  double mhd[NDIM];
  //Particl number flux
  flux[RHO][j][i] = S->P[RHO][j][i]*S->ucon[dir][j][i];
  
  mhd_calc(S, i, j, dir, mhd);
  // MHD stress-energy tensor with first index up and second index down
  flux[UU][j][i] = mhd[0] + flux[RHO][j][i];
  flux[U1][j][i] = mhd[1];
  flux[U2][j][i] = mhd[2];
  flux[U3][j][i] = mhd[3];
  //Dual of Maxwell Tensor
  flux[B1][j][i] = S->bcon[1][j][i]*S->ucon[dir][j][i] - S->bcon[dir][j][i]*S->ucon[1][j][i]; 
  flux[B2][j][i] = S->bcon[2][j][i]*S->ucon[dir][j][i] - S->bcon[dir][j][i]*S->ucon[2][j][i];
  flux[B3][j][i] = S->bcon[3][j][i]*S->ucon[dir][j][i] - S->bcon[dir][j][i]*S->ucon[3][j][i];

#if ELECTRONS
  for (int idx = KEL0; idx < NVAR; idx++)  
    flux[idx][j][i] = flux[RHO][j][i]*S->P[idx][j][i];
  flux[KTOT][j][i] = flux[RHO][j][i]*S->P[KTOT][j][i];
#endif
  
  PLOOP flux[ip][j][i] *= G->gdet[loc][j][i];
}

void prim_to_flux_vec(struct GridGeom *G, struct FluidState *S, int dir, int loc, int jstart, int jstop, int istart, int istop, GridPrim flux)
{
#pragma omp parallel
{
#pragma omp for collapse(2) nowait
  ZSLOOP(jstart, jstop, istart, istop)
  {
    double mhd[NDIM];

    flux[RHO][j][i] = S->P[RHO][j][i] * S->ucon[dir][j][i] * G->gdet[loc][j][i];
    mhd_calc(S, i, j, dir, mhd);
    // MHD stress-energy tensor w/ first index up, second index down
    flux[UU][j][i] = mhd[0] * G->gdet[loc][j][i] + flux[RHO][j][i];
    flux[U1][j][i] = mhd[1] * G->gdet[loc][j][i];
    flux[U2][j][i] = mhd[2] * G->gdet[loc][j][i];
    flux[U3][j][i] = mhd[3] * G->gdet[loc][j][i];
  }

#pragma omp for collapse(2) nowait
  ZSLOOP(jstart, jstop, istart, istop)
  {
    // Dual of Maxwell tensor
    flux[B1][j][i] = (S->bcon[1][j][i] * S->ucon[dir][j][i] - S->bcon[dir][j][i] * S->ucon[1][j][i]) * G->gdet[loc][j][i];
    flux[B2][j][i] = (S->bcon[2][j][i] * S->ucon[dir][j][i] - S->bcon[dir][j][i] * S->ucon[2][j][i]) * G->gdet[loc][j][i];
    flux[B3][j][i] = (S->bcon[3][j][i] * S->ucon[dir][j][i] - S->bcon[dir][j][i] * S->ucon[3][j][i]) * G->gdet[loc][j][i];

  }

#if ELECTRONS
#pragma omp for collapse(2)
  ZSLOOP(jstart, jstop, istart, istop)
  {
    // RHO already includes a factor of gdet!
    for (int idx = KEL0; idx < NVAR ; idx++)
      flux[idx][j][i] = flux[RHO][j][i]*S->P[idx][j][i];
    flux[KTOT][j][i] = flux[RHO][j][i]*S->P[KTOT][j][i];
  }
#endif
}
}

// Calculate magnetic field four-vector
inline void bcon_calc(struct FluidState *S, int i, int j)
{
  S->bcon[0][j][i] = S->P[B1][j][i]*S->ucov[1][j][i] + S->P[B2][j][i]*S->ucov[2][j][i] + S->P[B3][j][i]*S->ucov[3][j][i];
  for (int mu = 1; mu < 4; mu++)
    S->bcon[mu][j][i] = (S->P[B1-1+mu][j][i] + S->bcon[0][j][i]*S->ucon[mu][j][i])/S->ucon[0][j][i];
}

// Find gamma-factor wrt normal observer
inline double mhd_gamma_calc(struct GridGeom *G, struct FluidState *S, int i, int j, int loc)
{
  double qsq = G->gcov[loc][1][1][j][i]*S->P[U1][j][i]*S->P[U1][j][i]
      + G->gcov[loc][2][2][j][i]*S->P[U2][j][i]*S->P[U2][j][i]
      + G->gcov[loc][3][3][j][i]*S->P[U3][j][i]*S->P[U3][j][i]
      + 2.*(G->gcov[loc][1][2][j][i]*S->P[U1][j][i]*S->P[U2][j][i]
          + G->gcov[loc][1][3][j][i]*S->P[U1][j][i]*S->P[U3][j][i]
          + G->gcov[loc][2][3][j][i]*S->P[U2][j][i]*S->P[U3][j][i]);

#if DEBUG
  if (qsq < 0.)
  {
    if (fabs(qsq) > 1.E-10) // Then assume not just machine precision
    {
      fprintf(stderr, "gamma_calc():  failed: [%i %i] qsq = %28.18e \n", i, j, qsq);
      fprintf(stderr, "v[1-3] = %28.18e %28.18e %28.18e  \n", S->P[U1][j][i], S->P[U2][j][i], S->P[U3][j][i]);
      return 1.0;
    } 
    else
      qsq = 1.E-10; // Set floor
  }
#endif
  
  return sqrt(1. + qsq);
}

// Find contravariant four-velocity
inline void ucon_calc(struct GridGeom *G, struct FluidState *S, int i, int j, int loc)
{
  double gamma = mhd_gamma_calc(G, S, i, j, loc);
  double alpha = G->lapse[loc][j][i];
  S->ucon[0][j][i] = gamma/alpha;
  for (int mu = 1; mu < NDIM; mu++)
    S->ucon[mu][j][i] = S->P[U1+mu-1][j][i] - gamma*alpha*G->gcon[loc][0][mu][j][i];
}

// Calculate ucon, ucov, bcon, bcov from primitive variables
inline void get_state(struct GridGeom *G, struct FluidState *S, int i, int j, int loc)
{
  ucon_calc(G, S, i, j, loc);
  lower_grid(S->ucon, S->ucov, G, i, j, loc);
  bcon_calc(S, i, j);
  lower_grid(S->bcon, S->bcov, G, i, j, loc);
}

// Calculate ucon, ucov, bcon, bcov from primitive variables, over given range
// Note same range convention as ZSLOOP and other *_vec functions
void get_state_vec(struct GridGeom *G, struct FluidState *S, int loc, int jstart, int jstop, int istart, int istop)
{
#pragma omp parallel
  {
#pragma omp for collapse(2)
    ZSLOOP(jstart, jstop, istart, istop)
      ucon_calc(G, S, i, j, loc);

#pragma omp for collapse(2)
    ZSLOOP(jstart, jstop, istart, istop)
      lower_grid(S->ucon, S->ucov, G, i, j, loc);

#pragma omp for collapse(2)
    ZSLOOP(jstart, jstop, istart, istop)
      bcon_calc(S, i, j);

#pragma omp for collapse(2)
    ZSLOOP(jstart, jstop, istart, istop)
      lower_grid(S->bcon, S->bcov, G, i, j, loc);
  }
}

// Calculate components of magnetosonic velocity from primitive variables
inline void mhd_vchar(struct GridGeom *G, struct FluidState *S, int i, int j, int loc, int dir, GridDouble cmax, GridDouble cmin)
{
  double discr, vp, vm, bsq, ee, ef, va2, cs2, cms2, rho, u;
  double Acov[NDIM], Bcov[NDIM], Acon[NDIM], Bcon[NDIM];
  double Asq, Bsq, Au, Bu, AB, Au2, Bu2, AuBu, A, B, C;

  DLOOP1
    Acov[mu] = 0.;
  Acov[dir] = 1.;

  DLOOP1
    Bcov[mu] = 0.;
  Bcov[0] = 1.;

  DLOOP1
  {
    Acon[mu] = 0.;
    Bcon[mu] = 0.;
  }

  DLOOP2
  {
    Acon[mu] += G->gcon[loc][mu][nu][j][i]*Acov[nu];
    Bcon[mu] += G->gcon[loc][mu][nu][j][i]*Bcov[nu];
  }

  // Find fast magnetosonic velocity
  bsq = bsq_calc(S, i, j);
  rho = fabs(S->P[RHO][j][i]);
  u = fabs(S->P[UU][j][i]);
  ef = rho + gam*u;
  ee = bsq + ef;
  va2 = bsq/ee;
  cs2 = gam*(gam - 1.)*u/ef;

  cms2 = cs2 + va2 - cs2*va2;

  cms2 = (cms2 < 0) ? SMALL : cms2;
  cms2 = (cms2 > 1) ? 1 : cms2;

  // Require that speed of wave measured by observer q->ucon is cms2
  Asq = dot(Acon, Acov);
  Bsq = dot(Bcon, Bcov);
  Au = Bu = 0.;
  DLOOP1 
  {
    Au += Acov[mu]*S->ucon[mu][j][i];
    Bu += Bcov[mu]*S->ucon[mu][j][i];
  }
  AB = dot(Acon, Bcov);
  Au2 = Au*Au;
  Bu2 = Bu*Bu;
  AuBu = Au*Bu;

  A = Bu2 - (Bsq + Bu2)*cms2;
  B = 2.*(AuBu - (AB + AuBu)*cms2);
  C = Au2 - (Asq + Au2)*cms2;

  discr = B*B - 4.*A*C;
  discr = (discr < 0.) ? 0. : discr;
  discr = sqrt(discr);

  vp = -(-B + discr)/(2.*A);
  vm = -(-B - discr)/(2.*A);

  cmax[j][i] = (vp > vm) ? vp : vm;
  cmin[j][i] = (vp > vm) ? vm : vp;
}

// Source terms for equations of motion
inline void get_fluid_source(struct GridGeom *G, struct FluidState *S, GridPrim *dU)
{
#if WIND_TERM
  static struct FluidState *dS;
  static int firstc = 1;
  if (firstc)
  {
    dS = calloc(1,sizeof(struct FluidState)); 
    firstc = 0;
  }
#endif

#pragma omp parallel for collapse(2)
  ZLOOP 
  {
    double mhd[NDIM][NDIM];
    DLOOP1 mhd_calc(S, i, j, mu, mhd[mu]);
    // Contract mhd stress tensor with connection
    PLOOP (*dU)[ip][j][i] = 0.;
    DLOOP2
      for (int gam = 0; gam < NDIM; gam++)
        (*dU)[UU+gam][j][i] += mhd[mu][nu]*G->conn[nu][gam][mu][j][i];
    PLOOP (*dU)[ip][j][i] *= G->gdet[CENT][j][i];
  }

  // Add a small "wind" source term in RHO,UU
#if WIND_TERM
#pragma omp parallel for simd
  ZLOOP
  {
    /* need coordinates to evaluate particle addtn rate */
    double X[NDIM];
    coord(i, j, CENT, X);
    double r, th;
    bl_coord(X, &r, &th);
    double cth = cos(th);

    /* here is the rate at which we're adding particles */
    /* this function is designed to concentrate effect in the funnel in black hole evolutions */
    double drhopdt = 2.e-4*cth*cth*cth*cth/pow(1. + r*r,2);
    dS->P[RHO][j][i] = drhopdt ;
    double Tp = 10. ;  /* temp, in units of c^2, of new plasma */
    dS->P[UU][j][i] = drhopdt*Tp*3. ;
    /* Leave P[U{1,2,3}]=0 to add in particles in normal observer frame */
    /* Likewise leave P[BN]=0 */
  }

  /* add in plasma to the T^t_a component of the stress-energy tensor */
  /* notice that U already contains a factor of sqrt{-g} */
  get_state_vec(G, dS, CENT, 0, N2-1, 0, N1-1);
  prim_to_flux_vec(G, dS, 0, CENT, 0, N2-1, 0, N1-1, dS->U);

#pragma omp parallel for simd collapse(2)
  PLOOP 
    ZLOOP
      (*dU)[ip][j][i] += dS->U[ip][j][i];
#endif
}

// Returns b.b (twice magnetic pressure)
inline double bsq_calc(struct FluidState *S, int i, int j)
{
  double bsq = 0.;
  DLOOP1
    bsq += S->bcon[mu][j][i]*S->bcov[mu][j][i];
  
  return bsq;
}
