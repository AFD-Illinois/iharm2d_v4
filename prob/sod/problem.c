/*---------------------------------------------------------------------------------

  PROBLEM.C
 
  -Read and store problem-specific parameters from parameter file
  -Initialize Sod shock tube problem

-----------------------------------------------------------------------------------*/

#include "decs.h"

// Local declarations
#define FML_DBL_OUT "%28.18e"
#define FML_INT_OUT "%10d"
#define STRING_OUT "%15s"

static double tscale;

// Assign problem param names and pointers
void set_problem_params()
{
	set_param("tscale", &tscale);
}

// Save problem specific details
void save_problem_data(FILE *fp)
{
	fprintf(fp, FML_DBL_OUT, tscale);
}

// Initializing mean state and perturbations based on nmode
void init(struct GridGeom *G, struct FluidState *S)
{
   double X[NDIM];

  set_grid(G);
  LOG("Set grid");

  // Make problem nonrelativistic
  tf /= tscale;
  dt /= tscale;
  DTd /= tscale;
  DTl /= tscale;  

  ZLOOP {

    coord(i, j, CENT, X);

    S->P[RHO][j][i] = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.125;
    S->P[UU][j][i]  = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.1;
    S->P[U1][j][i]  = 0.;
    S->P[U2][j][i]  = 0.;
    S->P[U3][j][i]  = 0.;
    S->P[B1][j][i]  = 0.;
    S->P[B2][j][i]  = 0.;
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
