/*---------------------------------------------------------------------------------

  PROBLEM.C
 
  -Read and store problem-specific parameters from parameter file
  -Initialize Fishbone-Moncrief torus (FM)
  -Seed with no field as this is a test problem

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

static double rin, rmax;
static double u_jitter;

// Assign problem param names and pointers
void set_problem_params()
{
  set_param("rin", &rin);
  set_param("rmax", &rmax);
  set_param("u_jitter", &u_jitter);
}

// Save problem specific details
void save_problem_data(FILE *fp)
{
  fprintf(fp, STRING_OUT, "torus_equilibrium");
  fprintf(fp, FML_DBL_OUT, rin);
  fprintf(fp, FML_DBL_OUT, rmax);
  fprintf(fp, FML_DBL_OUT, u_jitter);
}


// Initializing matter
void init(struct GridGeom *G, struct FluidState *S)
{
    
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
    // Region inside magnetized torus; u^i is calculated in Boyer-Lindquist coordinates, as per FM, so it needs to be transformed at the end
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
