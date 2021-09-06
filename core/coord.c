/*---------------------------------------------------------------------------------

  COORD.C

  -SET GRID POINTS AT CENTER, CORNER AND FACES 
  -EVALUATE BL R AND TH FROM KS
  -COMPUTE TRANSFORMATION MATRIX FOR KS->MKS OR KS->FMKS
  -COMPUTE METRIC COEFFICIENTS IN MKS/FMKS
  -COMPUTE LIGHT-CROSSING TIME
  -INITIALIZE FAILURE FLAGS TO ZERO

---------------------------------------------------------------------------------*/

/*
 *      -- given the indices i,j and location in the cell, return with
 *         the values of X1,X2 there;
 *      -- the locations are defined by :
 *          -----------------------
 *          |                     |
 *          |                     |
 *          |FACE1   CENT         |
 *          |                     |
 *          |CORN    FACE2        |
 *          ----------------------
 *
 */

#include "decs.h"

double thG_of_X(const double X[NDIM]);
void thJ_of_X(const double X[NDIM], double *y, double* thJ);
double r_of_X(const double X[NDIM]);
double th_of_X(const double X[NDIM]);

// Set coordinate values at grid loc [i,j,LOC]
inline void coord(int i, int j, int loc, double *X)
{
  X[0] = 0; // Make sure all memory passed in is initialized
  if (loc == FACE1)
  {
    X[1] = startx[1] + (i - NG) * dx[1];
    X[2] = startx[2] + (j + 0.5 - NG) * dx[2];
  }
  else if (loc == FACE2) 
  {
    X[1] = startx[1] + (i + 0.5 - NG) * dx[1];
    X[2] = startx[2] + (j - NG) * dx[2];
  } 
  else if (loc == CENT)
  {
    X[1] = startx[1] + (i + 0.5 - NG) * dx[1];
    X[2] = startx[2] + (j + 0.5 - NG) * dx[2];
  } 
  else if (loc == CORN)
  {
    X[1] = startx[1] + (i - NG) * dx[1];
    X[2] = startx[2] + (j - NG) * dx[2];
  }

#if DEBUG
  else 
  {
    fprintf(stderr, "Invalid coordinate location!\n");
    exit(-1);
  }
#endif
}

// Computes theta_G from X2
inline double thG_of_X(const double X[NDIM])
{
  return M_PI*X[2] + ((1. - hslope)/2.)*sin(2.*M_PI*X[2]);
}

// Computes theta_J from X2
inline void thJ_of_X(const double X[NDIM], double *y, double* thJ)
{
  *y = 2*X[2] - 1.;
  *thJ = poly_norm*(*y)*(1. + pow((*y)/poly_xt,poly_alpha)/(poly_alpha+1.)) +
    0.5*M_PI;
}

// Computes r from X1
inline double r_of_X(const double X[NDIM])
{
  return exp(X[1]);
}

// Computes theta from (X1,X2)
inline double th_of_X(const double X[NDIM])
{
  double thG = thG_of_X(X);

#if DEREFINE_POLES
  double y, thJ;
  thJ_of_X(X, &y, &thJ);
  return thG + exp(mks_smooth*(startx[1] - X[1]))*(thJ - thG);
#else
  return thG;
#endif
}

// Boyer-Lindquist coordinate of point X
inline void bl_coord(const double X[NDIM], double *r, double *th)
{
  *r = r_of_X(X);
  *th = th_of_X(X);

  // Avoid singularity at polar axis
#if COORDSINGFIX
  if (fabs(*th) < SINGSMALL) {
    if ((*th) >= 0)
      *th = SINGSMALL;
    if ((*th) < 0)
      *th = -SINGSMALL;
  }
  if (fabs(M_PI - (*th)) < SINGSMALL) {
    if ((*th) >= M_PI)
      *th = M_PI + SINGSMALL;
    if ((*th) < M_PI)
      *th = M_PI - SINGSMALL;
  }
#endif
}

// Computes transformation matrix for KS->MKS and KS->FMKS
inline void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM])
{
  memset(dxdX, 0, NDIM*NDIM*sizeof(double));

#if METRIC == MINKOWSKI
  for (int mu = 0; mu < NDIM; mu++)
  {
    dxdX[mu][mu] = 1.;
  }
#elif METRIC == MKS && !DEREFINE_POLES
  dxdX[0][0] = 1.;
  dxdX[1][1] = exp(X[1]);
  dxdX[2][2] = M_PI - (hslope - 1.)*M_PI*cos(2.*M_PI*X[2]);
  dxdX[3][3] = 1.;
#elif METRIC == MKS && DEREFINE_POLES
  dxdX[0][0] = 1.;
  dxdX[1][1] = exp(X[1]);
  dxdX[2][1] = -exp(mks_smooth*(startx[1]-X[1]))*mks_smooth*(
    M_PI/2. -
    M_PI*X[2] +
    poly_norm*(2.*X[2]-1.)*(1+(pow((-1.+2*X[2])/poly_xt,poly_alpha))/(1 + poly_alpha)) -
    1./2.*(1. - hslope)*sin(2.*M_PI*X[2])
    );
  dxdX[2][2] = M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]) +
    exp(mks_smooth*(startx[1]-X[1]))*(
      -M_PI +
      2.*poly_norm*(1. + pow((2.*X[2]-1.)/poly_xt,poly_alpha)/(poly_alpha+1.)) +
      (2.*poly_alpha*poly_norm*(2.*X[2]-1.)*pow((2.*X[2]-1.)/poly_xt,poly_alpha-1.))/((1.+poly_alpha)*poly_xt) -
      (1.-hslope)*M_PI*cos(2.*M_PI*X[2])
      );
  dxdX[3][3] = 1.;
#else
#error "Unsupported metric!"
#endif
}

// Computes covariant metric in KS
void gcov_func(double X[NDIM], double gcov[NDIM][NDIM])
{
  memset(gcov, 0, NDIM*NDIM*sizeof(double));

  #if METRIC == MINKOWSKI
  gcov[0][0] = -1.;
  for (int j = 1; j < NDIM; j++) {
    gcov[j][j] = 1.;
  }
  #else //Everything else is covered in set_dxdX
  double sth, cth, s2, rho2;
  double r, th;

  bl_coord(X, &r, &th);

  cth = cos(th);
  sth = sin(th);

  s2 = sth*sth;
  rho2 = r*r + a*a*cth*cth;

  gcov[0][0] = -1. + 2.*r/rho2;
  gcov[0][1] = 2.*r/rho2;
  gcov[0][3] = -2.*a*r*s2/rho2;

  gcov[1][0] = gcov[0][1];
  gcov[1][1] = 1. + 2.*r/rho2;
  gcov[1][3] = -a*s2*(1. + 2.*r/rho2);

  gcov[2][2] = rho2;

  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];
  gcov[3][3] = s2*(rho2 + a*a*s2*(1. + 2.*r/rho2));

  // Apply coordinate transformation to code coordinates X
  double dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);

  double gcov_ks[NDIM][NDIM];
  memcpy(gcov_ks, gcov, NDIM*NDIM*sizeof(double));
  memset(gcov, 0, NDIM*NDIM*sizeof(double));

  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      for (int lam = 0; lam < NDIM; lam++) {
        for (int kap = 0; kap < NDIM; kap++) {
          gcov[mu][nu] += gcov_ks[lam][kap]*dxdX[lam][mu]*dxdX[kap][nu];
        }
      }
    }
  }
  #endif // METRIC
}

// Establish X coordinates
void set_points()
{
#if METRIC == MINKOWSKI
  startx[1] = x1Min;
  startx[2] = x2Min;
  dx[1] = (x1Max - x1Min)/N1TOT;
  dx[2] = (x2Max - x2Min)/N2TOT;
#elif METRIC == MKS
  // Set Rin such that we have 5 zones completely inside the event horizon
  // If xeh = log(Rhor), xin = log(Rin), and xout = log(Rout),
  // then we want xeh = xin + 5.5 * (xout - xin) / N1TOT, or solving/replacing:
  Rin = exp((N1TOT * log(Rhor) / 5.5 - log(Rout)) / (-1. + N1TOT / 5.5));

  startx[1] = log(Rin);
  if (startx[1] < 0.0) ERROR("Not enough radial zones! Increase N1!");
  startx[2] = 0.;

  dx[1] = log(Rout/Rin)/N1TOT;
  dx[2] = 1./N2TOT;

#if DEREFINE_POLES
  poly_norm = 0.5*M_PI*1./(1. + 1./(poly_alpha + 1.)*
                           1./pow(poly_xt, poly_alpha));
#endif
#endif // METRIC
}

// Sets the grid struct G
void set_grid(struct GridGeom *G)
{

  // Set up boundaries, steps in coordinate grid
  set_points();
  dV = dx[1]*dx[2];

#if !INTEL_WORKAROUND
#pragma omp parallel for collapse(2)
#endif
  JSLOOP(-NG, N2 - 1 + NG) {
    ISLOOP(-NG, N1 - 1 + NG) {
      set_grid_loc(G, i, j, CENT);
      set_grid_loc(G, i, j, CORN);
      set_grid_loc(G, i, j, FACE1);
      set_grid_loc(G, i, j, FACE2);

      // Connection only needed at zone center
      conn_func(G, i, j);
    }
  }
}

// Makes necessary function calls to set grid at various LOC
inline void set_grid_loc(struct GridGeom *G, int i, int j, int loc)
{
  double X[NDIM];
  double gcov[NDIM][NDIM], gcon[NDIM][NDIM];

  coord(i, j, loc, X);
  gcov_func(X, gcov);
  G->gdet[loc][j][i] = gcon_func(gcov, gcon);
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      G->gcov[loc][mu][nu][j][i] = gcov[mu][nu];
      G->gcon[loc][mu][nu][j][i] = gcon[mu][nu];
    }
  }
  G->lapse[loc][j][i] = 1./sqrt(-G->gcon[loc][0][0][j][i]);
}

// Initializes flags and fails to zero
void zero_arrays()
{
  ZLOOPALL
  {
    pflag[j][i] = 0;
    fail_save[j][i] = 0;
  }
}
