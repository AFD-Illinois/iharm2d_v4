/*---------------------------------------------------------------------------------

  PROBLEM.C
 
  -Read and store problem-specific parameters from parameter file
  -Initialize sound wave (courtesy of bhlight)

-----------------------------------------------------------------------------------*/

#include "decs.h"
#include <complex.h>

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

	// Mean (background) state
	double rho0 = 1.;
	double u0 = 0.01;
	double u10 = 0.1;

	// Perturbations to mean state
	double drho = 1.;
	double du = 0.;
	double du1 = 0.;

	// Wave-vector
	double k = 2.*M_PI;
	double amp = 0.01;

	// Rotate solution
	double theta = M_PI/4.;
	k *= sqrt(2.);

	// Set grid point, metric, connection coefficients and light-crossing time
	set_grid(G);
	LOG("Grid set");

	ZLOOP {
		coord(i, j, CENT, X);

		double mode = amp*cos(k*cos(theta)*X[1] + k*sin(theta)*X[2]);

		S->P[RHO][j][i]	= rho0 + drho*mode;
		S->P[UU][j][i] 	= u0   + du*mode;
		S->P[U1][j][i] 	= (u10 + du1*mode)*cos(theta);
		S->P[U2][j][i] 	= (u10 + du1*mode)*sin(theta);
		S->P[U3][j][i] 	= 0.;
		S->P[B1][j][i] 	= 0.;
		S->P[B2][j][i] 	= 0.;
		S->P[B3][j][i] 	= 0.;

	}

	// Enforce boundary conditions
	set_bounds(G, S);

	LOG("Finished init()");
}
