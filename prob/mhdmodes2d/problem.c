/*---------------------------------------------------------------------------------

  PROBLEM.C
 
  -Read and store problem-specific parameters from parameter file
  -Initialize MHD mode (Entropy, Slow, Alfven, Fast) in 2D

-----------------------------------------------------------------------------------*/

#include "decs.h"
#include <complex.h>

// Local declarations
#define FML_DBL_OUT "%28.18e"
#define FML_INT_OUT "%10d"
#define STRING_OUT "%15s"

static int nmode;

// Assign problem param names and pointers
void set_problem_params()
{
	set_param("nmode", &nmode);
}

// Save problem specific details
void save_problem_data(FILE *fp)
{
	fprintf(fp, FML_INT_OUT, nmode);
}

// Initializing mean state and perturbations based on nmode
void init(struct GridGeom *G, struct FluidState *S)
{
	double X[NDIM];

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
	double u0   = 1.;
	double B10  = 1.;
    double B20  = 0.;
    double B30  = 0.; 

    // Wave-vector
	double k1 = 2.*M_PI;
	double k2 = 2.*M_PI;
	double amp = 1.e-4;

	complex omega, drho, du, du1, du2, du3, dB1, dB2, dB3;

    // Default value 0
    omega = 0.;
    drho  = 0.;
    du    = 0.;
    du1   = 0.;
    du2   = 0.;
    du3   = 0.;
    dB1   = 0.;
    dB2   = 0.;
    dB3   = 0.;

    // Eigenmode
    if (nmode == 0) { // Entropy
        omega = 2.*M_PI/5.*I;
        drho = 1.;
    } else if (nmode == 1) { // Slow
        omega = 2.41024185339*I;
        drho  = 0.558104461559;
        du    = 0.744139282078;
        du1   = -0.277124827421;
        du2   = 0.0630348927707;
        dB1   = -0.164323721928;
        dB2   = 0.164323721928;
    } else if (nmode == 2) { // Alfven
        omega = 3.44144232573*I;
        du3   = 0.480384461415;
        dB3   = 0.877058019307;

    } else { // Fast
        omega = 5.53726217331*I;
        drho  = 0.476395427447;
        du    = 0.635193903263;
        du1   = -0.102965815319;
        du2   = -0.316873207561;
        dB1   = 0.359559114174;
        dB2   = -0.359559114174;
    }

    // Override tf and the dump and log intervals
    if (nmode > 0) tf = 2.*M_PI/fabs(cimag(omega));
    DTd = tf/5.;
    DTl = tf/2.;

    // Set grid points, metric, connection coefficients and light-crossing time
    set_grid(G);
    LOG("Grid set");

    // Initialize primitives
    ZLOOP {
        coord(i, j, CENT, X);

        double mode = amp*cos(k1*X[1] + k2*X[2]);
        
        S->P[RHO][j][i] = rho0 + creal(drho*mode);
        S->P[UU][j][i]  = u0 + creal(du*mode);
        S->P[U1][j][i]  = creal(du1*mode);
        S->P[U2][j][i]  = creal(du2*mode);
        S->P[U3][j][i]  = creal(du3*mode);
        S->P[B1][j][i]  = B10 + creal(dB1*mode);
        S->P[B2][j][i]  = B20 + creal(dB2*mode);
        S->P[B3][j][i]  = B30 + creal(dB3*mode);

    }

    //Enforce boundary conditions
    set_bounds(G, S);

    LOG("Finished init()");

}
