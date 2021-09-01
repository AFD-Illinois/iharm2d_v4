/*---------------------------------------------------------------------------------

  PROBLEM.C
 
  -Read and store problem-specific parameters from parameter file
  -Initialize Fishbone-Moncrief torus (FM)
  -Seed it with a SANE or MAD poloidal magnetic field

-----------------------------------------------------------------------------------*/

#include "bl_coord.h"
#include "decs.h"


// Local declarations
double lfish_calc(double rmax);

#define SANE 0
#define MAD 1

#define FML_DBL_OUT "%28.18e"
#define FML_INT_OUT "%10d"
#define STRING_OUT "%15s"

static int maxr_normalization = 0;

static int mad_type;
static double beta;
static double rin, rmax;
static double u_jitter;

// Assign problem param names and pointers
void set_problem_params()
{
  set_param("rin", &rin);
  set_param("rmax", &rmax);
  set_param("u_jitter", &u_jitter);

  set_param("mad_type", &mad_type);
  set_param("beta", &beta);
}

// Save problem specific details
void save_problem_data(FILE *fp)
{
  fprintf(fp, FML_INT_OUT, mad_type);
  fprintf(fp, STRING_OUT, "torus");
  fprintf(fp, FML_DBL_OUT, rin);
  fprintf(fp, FML_DBL_OUT, rmax);
  fprintf(fp, FML_DBL_OUT, beta);
  fprintf(fp, FML_DBL_OUT, u_jitter);
}


// Initializing matter and magnetic field
void init(struct GridGeom *G, struct FluidState *S)
{
  // Magnetic field
  double (*A)[N2 + 2*NG] = malloc(sizeof(*A) * (N1 + 2*NG));
  
  // FM parameters
  double l = lfish_calc(rmax);
  double kappa = 1.e-3;

  // Grid parameters
  Rhor = (1. + sqrt(1. - a*a));
  
  // Set grid points, metric, connection coefficients and light-crossing time
  set_grid(G);
  LOG("Grid set");

  // Initializing matter variables
  double rhomax = 0.;
  double umax = 0.;
  ZSLOOP(-1, N2, -1, N1) {
    double X[NDIM];
    coord(i, j, CENT, X);
    double r, th;
    bl_coord(X, &r, &th);

    double sth = sin(th);
    double cth = cos(th);

    // Calculate lnh
    double DD = r * r - 2. * r + a * a;
    double AA = (r * r + a * a) * (r * r + a * a) -
             DD * a * a * sth * sth;
    double SS = r * r + a * a * cth * cth;

    double thin = M_PI / 2.;
    double sthin = sin(thin);
    double cthin = cos(thin);

    double DDin = rin * rin - 2. * rin + a * a;
    double AAin = (rin * rin + a * a) * (rin * rin + a * a)
             - DDin * a * a * sthin * sthin;
    double SSin = rin * rin + a * a * cthin * cthin;

    // Equation (3.6) in https://ui.adsabs.harvard.edu/abs/1976ApJ...207..962F
    double lnh;
    if (r >= rin) {
      lnh =
          0.5 *
          log((1. +
         sqrt(1. +
              4. * (l * l * SS * SS) * DD / (AA * AA * sth * sth)))
        / (SS * DD / AA))
          - 0.5 * sqrt(1. +
           4. * (l * l * SS * SS) * DD /
           (AA * AA * sth * sth))
          - 2. * a * r * l / AA -
          (0.5 *
           log((1. +
          sqrt(1. +
               4. * (l * l * SSin * SSin) * DDin /
               (AAin * AAin * sthin * sthin))) /
         (SSin * DDin / AAin))
           - 0.5 * sqrt(1. +
            4. * (l * l * SSin * SSin) * DDin / (AAin * AAin * sthin * sthin))
           - 2. * a * rin * l / AAin);
    } else {
      lnh = 1.;
    }

    // Regions outside torus
    if (lnh < 0. || r < rin) {
      // Nominal values; real value set by fixup

      S->P[RHO][j][i] = 1.e-7 * RHOMIN;
      S->P[UU][j][i] = 1.e-7 * UUMIN;
      S->P[U1][j][i] = 0.;
      S->P[U2][j][i] = 0.;
      S->P[U3][j][i] = 0.;
    }
    /* Region inside magnetized torus; u^i is calculated in
     * Boyer-Lindquist coordinates, as per FM,
     * so it needs to be transformed at the end */
    else {
      double hm1 = exp(lnh) - 1.;
      double rho = pow(hm1 * (gam - 1.) / (kappa * gam),
               1. / (gam - 1.));
      double u = kappa * pow(rho, gam) / (gam - 1.);

      // Calculate u^phi from u_{(\phi)}, refer to Equation (3.3) in https://ui.adsabs.harvard.edu/abs/1976ApJ...207..962F
      double expm2chi = SS * SS * DD / (AA * AA * sth * sth);
      double up1 =
          sqrt((-1. +
          sqrt(1. + 4. * l * l * expm2chi)) / 2.);
      double up = 2. * a * r * sqrt(1. +
                 up1 * up1) / sqrt(AA * SS *
                 DD) +
          sqrt(SS / AA) * up1 / sth;


      S->P[RHO][j][i] = rho;
      if (rho > rhomax) rhomax = rho;
      u *= (1. + u_jitter * (ran_uniform() - 0.5));
      if (u > umax && r > rin) umax = u;
      S->P[UU][j][i] = u;
      S->P[U1][j][i] = 0.;
      S->P[U2][j][i] = 0.;
      S->P[U3][j][i] = up;

      // Convert from 4-velocity to 3-velocity
      coord_transform(G, S, i, j);
    }

    S->P[B1][j][i] = 0.;
    S->P[B2][j][i] = 0.;
    S->P[B3][j][i] = 0.;
  }

  // Normalize densities so that max(rho) = 1
#pragma omp parallel for simd
  ZLOOPALL
  {
    S->P[RHO][j][i] /= rhomax;
    S->P[UU][j][i] /= rhomax;
  }
  umax /= rhomax;
  rhomax = 1;
  fixup(G, S);
  set_bounds(G, S);

  // Find corner-centered vector potential
#pragma omp parallel for simd
  ZSLOOP(-NG+1, N2+NG-1, -NG+1, N1+NG-1)
  {
    double X[NDIM];
    coord(i, j, CORN, X);
    double r, th;
    bl_coord(X, &r, &th);

    double q;
    
    // Field in disk
    double rho_av = 0.25*(S->P[RHO][j][i] + S->P[RHO][j][i-1] + S->P[RHO][j-1][i] + S->P[RHO][j-1][i-1]);
    double uu_av = 0.25*(S->P[UU][j][i] + S->P[UU][j][i-1] + S->P[UU][j-1][i] + S->P[UU][j-1][i-1]);
    
    if (mad_type == SANE)
      q = rho_av/rhomax - 0.2;
    else if (mad_type == MAD)
      q = pow(sin(th), 3) * pow(r/rin, 3.) * exp(-r/400) * rho_av/rhomax - 0.2;
    else
    {
      printf("MAD type %i  not supported!\n", mad_type);
      exit(-1);
    }

    A[i][j] = 0.;
    if (q > 0.) A[i][j] = q;
  }

  // Calculate B-field and find bsq_max
  double bsq_max = 0.;
  double beta_min = 1e100;
  ZLOOP
  {
    double X[NDIM];
    coord(i, j, CORN, X);
    double r, th;
    bl_coord(X, &r, &th);
    
    // Flux-ct
    S->P[B1][j][i] = -(A[i][j] - A[i][j + 1] + A[i + 1][j] - A[i + 1][j + 1]) / (2. * dx[2] * G->gdet[CENT][j][i]);
    S->P[B2][j][i] = (A[i][j] + A[i][j + 1] - A[i + 1][j] - A[i + 1][j + 1]) / (2. * dx[1] * G->gdet[CENT][j][i]);
    S->P[B3][j][i] = 0.;

    get_state(G, S, i, j, CENT);
    double bsq_ij = bsq_calc(S, i, j);
    if (bsq_ij > bsq_max) bsq_max = bsq_ij;

    double beta_ij = (gam - 1.)*(S->P[UU][j][i])/(0.5*(bsq_ij+SMALL));
    if(beta_ij < beta_min) beta_min = beta_ij;
  }

  double norm = 0;
  if (!maxr_normalization)
  {
    // Ratio of max UU, beta
    double beta_act = (gam - 1.) * umax / (0.5 * bsq_max);
    
    // In plane only
    LOGN("Umax is %.10e", umax);
    LOGN("bsq_max is %.10e", bsq_max);
    LOGN("beta is %.10e", beta_act);
    norm = sqrt(beta_act / beta);
  } 
  else
  {
    // Beta_min = 100 normalization
    LOGN("Min beta in torus is %f", beta_min);
    norm = sqrt(beta_min / beta) ;
  }

  // Apply normalization
  LOGN("Normalization is %f\n", norm);
  ZLOOP 
  {
    S->P[B1][j][i] *= norm ;
    S->P[B2][j][i] *= norm ;
    // B3 is uniformly 0
  }

#if ELECTRONS
  init_electrons(G, S);
#endif

  // Enforce boundary conditions
  fixup(G, S);
  set_bounds(G, S);

  LOG("Finished init()");
}

// Find angular momentum density at pressure max radius
// See Equation (3.8) in https://ui.adsabs.harvard.edu/abs/1976ApJ...207..962F  
double lfish_calc(double r)
{
  return (((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
     ((-2. * a * r *
       (pow(a, 2) - 2. * a * sqrt(r) +
        pow(r,
      2))) / sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
      ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) + pow(a, 2) *
      (2. + r))) / sqrt(1 + (2. * a) / pow (r, 1.5) - 3. / r)))
    / (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
       (pow(a, 2) + (-2. + r) * r))
      );
}
