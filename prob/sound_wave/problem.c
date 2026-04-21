/*---------------------------------------------------------------------------------

  PROBLEM.C

  -Read and store problem-specific parameters from parameter file
  -Initialize sound wave

-----------------------------------------------------------------------------------*/

#include "decs.h"
#include <complex.h>

// Local declarations
#define FML_DBL_OUT "%28.18e"
#define FML_INT_OUT "%10d"
#define STRING_OUT  "%15s"

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

// Initializing mean state and sound-wave eigenmode perturbations
void init(struct GridGeom *G, struct FluidState *S)
{
	double X[NDIM];

	// Mean (background) state  ---  Gamma = 4/3 assumed in param.dat
	double rho0 = 1.0;
	double u0   = 1.0;
	double u10  = 0.0;

	// Sound-wave eigenvector, Euclidean-normalized
	// (derived from linearized relativistic Euler about the background above)
	double drho =  0.580429504019981;   //  sqrt(63/187)
	double du   =  0.773906005359975;   //  sqrt(112/187)
	double du1  =  0.253320866973317;   //  sqrt(12/187)

	// Wave-vector (one full wavelength in the box along the rotated direction)
	double k   = 2.0*M_PI;
	double amp = 1.0e-4;

	// Rotate solution so the wave propagates at 45 deg in the x1-x2 plane
	double theta = M_PI/4.;
	k *= sqrt(2.);

	// Set grid points, metric, connection coefficients and light-crossing time
	set_grid(G);
	LOG("Grid set");

	ZLOOP {
		coord(i, j, CENT, X);

		double mode = amp*cos(k*cos(theta)*X[1] + k*sin(theta)*X[2]);

		S->P[RHO][j][i] = rho0 + drho*mode;
		S->P[UU][j][i]  = u0   + du  *mode;
		S->P[U1][j][i]  = (u10 + du1*mode)*cos(theta);
		S->P[U2][j][i]  = (u10 + du1*mode)*sin(theta);
		S->P[U3][j][i]  = 0.;
		S->P[B1][j][i]  = 0.;
		S->P[B2][j][i]  = 0.;
		S->P[B3][j][i]  = 0.;
	}

	// Enforce boundary conditions
	set_bounds(G, S);

	LOG("Finished init()");
}