/*---------------------------------------------------------------------------------

  PROBLEM.C
 
  -Read and store problem-specific parameters from parameter file
  -Initialize 1d modes

-----------------------------------------------------------------------------------*/

#include "decs.h"
#include <complex.h>

// Local declarations
#define FML_DBL_OUT "%28.18e"
#define FML_INT_OUT "%10d"
#define STRING_OUT "%15s"

static int nmode;
static int idim;

// Assign problem param names and pointers
void set_problem_params()
{
	set_param("nmode", &nmode);
	set_param("idim", &idim);
}

// Save problem specific details
void save_problem_data(FILE *fp)
{
	fprintf(fp, FML_INT_OUT, nmode);
	fprintf(fp, FML_INT_OUT, idim);
}

// Initializing mean state and perturbations based on nmode
void init(struct GridGeom *G, struct FluidState *S)
{
	double X[NDIM];
	
	// Check if appropriate dimension has been provided
	if (idim != 1 && idim != 2 && idim !=3)
	{
		fprintf(stderr, "ERROR: Primary axis must be 1, 2 or 3!\n");
		exit(1);
	}

	// Check if appropriate mode has been provided
	if (nmode != 0 && nmode != 1 && nmode != 2 && nmode != 3)
	{
		fprintf(stderr, "ERROR: Not a valid mode. Modes are:\n"
		"\t0: ENTROPY\n"
		"\t1: SLOW\n"
		"\t2: ALFVEN\n"
		"\t3: FAST\n"
		"\tYou chose %d.\n", nmode);
		exit(1);
	}

	// Mean (background) state
	double rho0 = 1.;
	double u0 = 1.;
	double B10, B20, B30; // Will be initialized based on idim

	if (idim == 1)
	{
		B10 = 1.;
		B20 = 0.;
		B30 = 0.;
	}

	if (idim == 2)
	{
		B10 = 0.;
		B20 = 1.;
		B30 = 0.;
	}

	if (idim == 3)
	{
		B10 = 0.;
		B20 = 0.;
		B30 = 1.;
	}

	// Wave-vector
	double k1 = 2.*M_PI;
	double amp = 1.e-4;

	complex omega, drho, du, du1, du2, du3, dB1, dB2, dB3, dutemp, dBtemp;

	// Eigenmode
	if (nmode == 0) // Entropy wave
	{
		omega = 2.*M_PI/tf*I;
		drho = 1.;
		du = 0.;
		du1 = 0.;
		du2 = 0.;
		du3 = 0.;
		dB1 = 0.;
		dB2 = 0.;
		dB3 = 0.;
	}

	else if (nmode == 1) // Slow waves
	{
		omega = 2.74220688339*I;
		drho = 0.580429492464;
		du = 0.773905989952;
		du1 = -0.253320198552;
		du2 = 0.;
		du3 = 0.;
		dB1 = 0.;
		dB2 = 0.;
		dB3 = 0.;
	}

	else if (nmode == 2) // Alfven waves
	{
		omega = 3.44144232573*I;
		drho = 0.;
		du = 0.;
		du1 = 0.;
		du2 = 0.480384461415;
		du3 = 0.;
		dB1 = 0.;
		dB2 = 0.877058019307;
		dB3 = 0.;
	}

	else if (nmode == 3) // Fast waves
	{
		omega = 3.44144232573*I;
		drho = 0.;
		du = 0.;
		du1 = 0.;
		du2 = 0.;
		du3 = 0.480384461415;
		dB1 = 0.;
		dB2 = 0.;
		dB3 = 0.877058019307;
	}

	// Change values for perturbations in velocities and magnetic fields for other idims
	if (idim == 2)
	{
		dutemp = du3;
		du3 = du2;
		du2 = du1;
		du1 = dutemp;
		dBtemp = dB3;
		dB3 = dB2;
		dB2 = dB1;
		dB1 = dBtemp;
	}

	if (idim == 3)
	{
		dutemp = du1;
		du1 = du2;
		du2 = du3;
		du3 = dutemp;
		dBtemp = dB1;
		dB1 = dB2;
		dB2 = dB3;
		dB3 = dBtemp;
	}

	// Change final time based on frequency
	if (nmode > 0) tf = 2.*M_PI/fabs(cimag(omega));
	DTd = tf/5.;
	DTl = tf/5.;

	// Set grid
	set_grid(G);

	LOG("Grid set");

	ZLOOP
	{
		coord(i, j, CENT, X);
		double mode = amp*cos(k1*X[idim]);		

		S->P[RHO][j][i] = rho0 + creal(drho*mode);
		S->P[UU][j][i] = u0 + creal(du*mode);
		S->P[U1][j][i] = creal(du1*mode);
		S->P[U2][j][i] = creal(du2*mode);
		S->P[U3][j][i] = creal(du3*mode);
		S->P[B1][j][i] = B10 + creal(dB1*mode);
		S->P[B2][j][i] = B20 + creal(dB2*mode);
		S->P[B3][j][i] = B30 + creal(dB3*mode);

	}

	// Enforce boundary conditions
	set_bounds(G, S);

	LOG("Finished init()");
}
