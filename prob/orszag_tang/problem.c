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

static double tscale;
static double phase;

// Assign problem param names and pointers
void set_problem_params()
{
	set_param("tscale", &tscale);
  set_param("phase", &phase);
}

// Save problem specific details
void save_problem_data(FILE *fp)
{
	fprintf(fp, FML_DBL_OUT, tscale);
	fprintf(fp, FML_DBL_OUT, phase);
}

// Initializing mean state and perturbations based on nmode
void init(struct GridGeom *G, struct FluidState *S)
{
   double X[NDIM];

  set_grid(G);
  LOG("Set grid");

  ZLOOP {

    coord(i, j, CENT, X);

    S->P[RHO][j][i] = 25./9.;
    S->P[UU][j][i]  = 5./3./(gam - 1.);
    S->P[U1][j][i]  = -sin(X[2] + phase);
    S->P[U2][j][i]  = sin(X[1] + phase);
    S->P[U3][j][i]  = 0.;
    S->P[B1][j][i]  = -sin(X[2] + phase);
    S->P[B2][j][i]  = sin(2.*(X[1] + phase));
    S->P[B3][j][i]  = 0.;

    // rescale by tscale. 0(non-relativistic) < tscale < 1(fully relativistic)
    S->P[UU][j][i] *= tscale * tscale;
    S->P[U1][j][i] *= tscale;
    S->P[U2][j][i] *= tscale;
    S->P[U3][j][i] *= tscale;
    S->P[B1][j][i] *= tscale;
    S->P[B2][j][i] *= tscale;
    S->P[B3][j][i] *= tscale;

  } // ZLOOP

  //Enforce boundary conditions
  set_bounds(G, S);
}
