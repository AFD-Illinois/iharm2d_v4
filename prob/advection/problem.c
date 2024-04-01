/*---------------------------------------------------------------------------------

  PROBLEM.C
 
  -Read and store problem-specific parameters from parameter file
  -Initialize Orszag-Tang problem

-----------------------------------------------------------------------------------*/

#include "decs.h"

// Local declarations
#define FML_DBL_OUT "%28.18e"
#define FML_INT_OUT "%10d"
#define STRING_OUT "%15s"

static double A;
static double r_sq;
static double v_init;

// Assign problem param names and pointers
void set_problem_params()
{
	set_param("A", &A);
	set_param("r_sq", &r_sq);
	set_param("v_init", &v_init);
}

// Save problem specific details
void save_problem_data(FILE *fp)
{
	fprintf(fp, FML_DBL_OUT, A);
	fprintf(fp, FML_DBL_OUT, r_sq);
	fprintf(fp, FML_DBL_OUT, v_init);
}

// Initializing mean state and perturbations based on nmode
void init(struct GridGeom *G, struct FluidState *S)
{
   double X[NDIM];

  set_grid(G);
  LOG("Set grid");

  double r_init = sqrt(r_sq);
  // Lorentz factor
  double gamma = sqrt(1./(1. - (v_init * v_init)));

  ZLOOP {

    coord(i, j, CENT, X);

    // Local variables to improve readability later
    double r = sqrt(pow((X[1] - 0.5), 2) + pow((X[2] - 0.5), 2));

    // Velocity along each axis
    double vi = sqrt(v_init * v_init / 2.);

    S->P[RHO][j][i] = 1. + (A * exp(-pow(r, 2) / pow(r_init, 2)));
    S->P[UU][j][i]  = 0.1;
    S->P[U1][j][i]  = gamma * vi;
    S->P[U2][j][i]  = gamma * vi;
    S->P[U3][j][i]  = 0.;
    S->P[B1][j][i]  = 0.;
    S->P[B2][j][i]  = 0.;
    S->P[B3][j][i]  = 0.;

  } // ZLOOP

  tf = sqrt(2.) / v_init;

  //Enforce boundary conditions
  set_bounds(G, S);
}
