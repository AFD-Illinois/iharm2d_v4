/*---------------------------------------------------------------------------------

  PROBLEM.C
 
  -Read and store problem-specific parameters from parameter file
  -Initialize Kelvin-Helmholtz problem

-----------------------------------------------------------------------------------*/

#include "decs.h"

// Local declarations
#define FML_DBL_OUT "%28.18e"
#define FML_INT_OUT "%10d"
#define STRING_OUT "%15s"

static double tscale, rho0, Drho, P0, u_flow, a_KH, sigma, amp, z1, z2;

// Assign problem param names and pointers
void set_problem_params()
{
	set_param("tscale", &tscale);
	set_param("rho0", &rho0);
	set_param("Drho", &Drho);
	set_param("P0", &P0);
	set_param("u_flow", &u_flow);
	set_param("a_KH", &a_KH);
	set_param("sigma", &sigma);
	set_param("amp", &amp);
	set_param("z1", &z1);
	set_param("z2", &z2);
}

// Save problem specific details
void save_problem_data(FILE *fp)
{
	fprintf(fp, FML_DBL_OUT, tscale);
	fprintf(fp, FML_DBL_OUT, rho0);
	fprintf(fp, FML_DBL_OUT, Drho);
	fprintf(fp, FML_DBL_OUT, P0);
	fprintf(fp, FML_DBL_OUT, u_flow);
	fprintf(fp, FML_DBL_OUT, a_KH);
	fprintf(fp, FML_DBL_OUT, sigma);
	fprintf(fp, FML_DBL_OUT, amp);
	fprintf(fp, FML_DBL_OUT, z1);
	fprintf(fp, FML_DBL_OUT, z2);
}


// Initializing mean state and perturbations based on nmode
void init(struct GridGeom *G, struct FluidState *S)
{
   double X[NDIM];

  set_grid(G);
  LOG("Set grid");

  ZLOOP {

    coord(i, j, CENT, X);

    double zdist1 = X[2] - z1;
    double zdist2 = X[2] - z2;

    S->P[RHO][j][i] = rho0 + Drho * 0.5 * (tanh(zdist1 / a_KH) - tanh(zdist2 / a_KH));
    S->P[UU][j][i]  = P0 / (gam - 1.) * tscale * tscale;
    S->P[U1][j][i]  = u_flow * (tanh(zdist1 / a_KH) - tanh(zdist2 / a_KH) - 1.) * tscale;
    S->P[U2][j][i]  = amp * sin(2. * M_PI * X[1]) * (exp(-(zdist1 * zdist1) / (sigma * sigma)) +
                        exp(-(zdist2 * zdist2) / (sigma * sigma))) * tscale;
    S->P[U3][j][i]  = 0.;
    S->P[B1][j][i]  = 0.;
    S->P[B2][j][i]  = 0.;
    S->P[B3][j][i]  = 0.;

  } // ZLOOP

  //Enforce boundary conditions
  set_bounds(G, S);
}
