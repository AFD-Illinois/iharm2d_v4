/*---------------------------------------------------------------------------------

  PROBLEM.C
 
  -Read and store problem-specific parameters from parameter file
  -Initialize Bondi flow

-----------------------------------------------------------------------------------*/

#include "bl_coord.h"
#include "decs.h"

// Local declarations
#define FML_DBL_OUT "%28.18e"
#define FML_INT_OUT "%10d"
#define STRING_OUT "%15s"

double C4, C3, n, K;

double mdot, rs;
// Assign problem param names and pointers
void set_problem_params()
{
	set_param("mdot", &mdot);
	set_param("rs", &rs);
	set_param("Rhor", &Rhor);
}

// Save problem specific details
void save_problem_data(FILE *fp)
{
	fprintf(fp, FML_DBL_OUT, mdot);
	fprintf(fp, FML_DBL_OUT, rs);
}

// Expression to roor find to solve for temperature
double get_Tfunc(double T, double r)
{
  return pow(1.+(1.+n)*T,2.)*(1.-2./r+pow(C4/r/r/pow(T,n),2.))-C3;
}

// Obtain T numerical with Newton-Raphson
double get_T(double r)
{
  double rtol = 1.e-12;
  double ftol = 1.e-14;
  double Tmin = 0.6*(sqrt(C3) - 1.)/(n + 1);
  double Tmax = pow(C4*sqrt(2./r/r/r),1./n);
  double f0, f1, fh;
  double T0, T1, Th;
  T0 = 0.6*Tmin;
  f0 = get_Tfunc(T0, r);
  T1 = Tmax;
  f1 = get_Tfunc(T1, r);

  if (f0*f1 > 0.) {
    printf("Failed solving for T at r = %e C4 = %e C3 = %e\n", r, C4, C3);
    exit(-1);
  }

  Th = (f1*T0 - f0*T1)/(f1 - f0);
  fh = get_Tfunc(Th, r);
  double epsT = rtol*(Tmin + Tmax);
  while (fabs(Th - T0) > epsT && fabs(Th - T1) > epsT && fabs(fh) > ftol) {
    if (fh*f0 < 0.) {
      T0 = Th;
      f0 = fh;
    } else {
      T1 = Th;
      f1 = fh;
    }

    Th = (f1*T0 - f0*T1)/(f1 - f0);
    fh = get_Tfunc(Th, r);
  }

  return Th;
}

// Convert 4-velocities to primitives - utilde
void fourvel_to_prim(double ucon[NDIM], GridPrim P,
  struct GridGeom *G, int i, int j)
{
  double alpha, beta[NDIM], gamma;

  alpha = 1.0/sqrt(-G->gcon[CENT][0][0][j][i]);
  beta[1] = alpha*alpha*G->gcon[CENT][0][1][j][i];
  beta[2] = alpha*alpha*G->gcon[CENT][0][2][j][i];
  beta[3] = alpha*alpha*G->gcon[CENT][0][3][j][i];
  gamma = ucon[0]*alpha;

  P[U1][j][i] = ucon[1] + beta[1]*gamma/alpha;
  P[U2][j][i] = ucon[2] + beta[2]*gamma/alpha;
  P[U3][j][i] = ucon[3] + beta[3]*gamma/alpha;
}

// Compute 0th component of 4-velocity from normalization condition
void set_ut(double ucon[NDIM], struct of_geom *geom)
{
    double AA, BB, CC;

    AA = geom->gcov[0][0];
    BB = 2.*(geom->gcov[0][1]*ucon[1] +
            geom->gcov[0][2]*ucon[2] +
            geom->gcov[0][3]*ucon[3]);
    CC = 1. + geom->gcov[1][1]*ucon[1]*ucon[1] +
            geom->gcov[2][2]*ucon[2]*ucon[2] +
            geom->gcov[3][3]*ucon[3]*ucon[3] +
        2. *(geom->gcov[1][2]*ucon[1]*ucon[2] +
            geom->gcov[1][3]*ucon[1]*ucon[3] +
            geom->gcov[2][3]*ucon[2]*ucon[3]);

    double discr = BB*BB - 4.*AA*CC;
    ucon[0] = (-BB - sqrt(discr))/(2.*AA);
}

// Calculate primitives at zone (i,j)
void get_prim_bondi(int i, int j, GridPrim P, struct GridGeom *G)
{
    static int firstc = 1;
    if (firstc) {
        n = 1./(gam - 1.); // Half the number of degrees of freedom of gas

        // Solution constants
        double uc = sqrt(1/(2.*rs));
        double Vc = sqrt(pow(uc,2)/(1. - 3.*pow(uc,2)));
        double Tc = -n*pow(Vc,2)/((n + 1.)*(n*pow(Vc,2) - 1.));
        C4 = uc*pow(rs,2)*pow(Tc,n);
        C3 = pow(1. + (1. + n)*Tc,2)*(1. - 2./rs + pow(uc, 2));
        K  = pow(4*M_PI*C4 / mdot, 1/n);
        firstc = 0;
    }

    // Obtain r,th
    double r, th, X[NDIM];
    coord(i, j, CENT, X);
    bl_coord(X, &r, &th);

    // Moved the horizon a little further out to avoid numerical issues
    // near the event horizon
    while (r < Rhor) {
        i++;
        coord(i, j, CENT, X);
        bl_coord(X, &r, &th);
    }

    // Obtain temperature, and consequently rho, u, and ur
    double T = get_T(r);
    double ur = -C4/(pow(T,n)*pow(r,2));
    double rho = pow(K, -n)*pow(T, n);
    double u = rho*T / (gam - 1.);
    double ucon_bl[NDIM], ucon_ks[NDIM], ucon_mks[NDIM];
    struct of_geom geom_bl;

    // Set metric in Boyer-Lindquist coordinates
    blgset(i, j, &geom_bl);

    // Compute velocity primitives
    DLOOP1 {
        ucon_bl[mu] = 0.;
        ucon_ks[mu] = 0.;
        ucon_mks[mu] = 0.;
    }
    ucon_bl[1] = ur;

    set_ut(ucon_bl, &geom_bl);
    bl_to_ks(X, ucon_bl, ucon_ks);

    double dxdX[NDIM][NDIM], dXdx[NDIM][NDIM];
    set_dxdX(X, dxdX);
    invert(&dxdX[0][0], &dXdx[0][0]);
    DLOOP2 {
        ucon_mks[mu] += dXdx[mu][nu]*ucon_ks[nu];
    }

    fourvel_to_prim(ucon_mks, P, G, i, j);

    P[RHO][j][i] = rho;
    P[UU][j][i]  = u;
    P[B1][j][i]  = 0.;
    P[B2][j][i]  = 0.;
    P[B3][j][i]  = 0.;

}

// Initializing mean state and perturbations based on nmode
void init(struct GridGeom *G, struct FluidState *S)
{
    set_grid(G);

    LOG("Grid set");

    ZLOOP {
        get_prim_bondi(i, j, S->P, G);
    }

    if (DEBUG) {
        printf("a = %e Rhor = %e\n", a, Rhor);

        printf("mdot = %e\n", mdot);
        printf("rs   = %e\n", rs);
        printf("n    = %e\n", n);
        printf("C4   = %e\n", C4);
        printf("C3   = %e\n", C3);
    }

    //Enforce boundary conditions
    fixup(G, S);
    set_bounds(G, S);

    LOG("Finished init()");
}

void bound_gas_prob_x1r(int i, int j, GridPrim  P, struct GridGeom *G)
{
    get_prim_bondi(i, j, P, G);
}