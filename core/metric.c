/*---------------------------------------------------------------------------------

  METRIC.C

  -Helper functions for metric tensors 
  -Compute 4x4 matrix minor, adjoint, determinant and inverse
  -Compute connection coefficients
  -Raise and lower rank-1 tensors
  -Take dot product of a contravariant and covariant rank-1 tensor

---------------------------------------------------------------------------------*/

#include "decs.h"

double MINOR(double m[16], int r0, int r1, int r2, int c0, int c1, int c2);
void adjoint(double m[16], double adjOut[16]);
double determinant(double m[16]);

inline double gcon_func(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM])
{
  double gdet = invert(&gcov[0][0],&gcon[0][0]);
  return sqrt(fabs(gdet));
}

inline void get_gcov(struct GridGeom *G, int i, int j, int loc, double gcov[NDIM][NDIM]) {
  DLOOP2 gcov[mu][nu] = G->gcov[loc][mu][nu][j][i];
}

inline void get_gcon(struct GridGeom *G, int i, int j, int loc, double gcon[NDIM][NDIM])
{
  DLOOP2 gcon[mu][nu] = G->gcon[loc][mu][nu][j][i];
}

// Calculate connection coefficient 
inline void conn_func(struct GridGeom *G, int i, int j)
{
  double tmp[NDIM][NDIM][NDIM];
  double X[NDIM], Xh[NDIM], Xl[NDIM];
  double gh[NDIM][NDIM];
  double gl[NDIM][NDIM];
  coord(i, j, CENT, X);

  for (int mu = 0; mu < NDIM; mu++) {
    for (int kap = 0; kap < NDIM; kap++) {
      Xh[kap] = X[kap];
    }
    for (int kap = 0; kap < NDIM; kap++) {
      Xl[kap] = X[kap];
    }
    Xh[mu] += DELTA;
    Xl[mu] -= DELTA;
    gcov_func(Xh, gh);
    gcov_func(Xl, gl);

    for (int lam = 0; lam < NDIM; lam++) {
      for (int nu = 0; nu < NDIM; nu++) {
        G->conn[lam][nu][mu][j][i] = (gh[lam][nu] - gl[lam][nu])/(Xh[mu] - 
                                                                  Xl[mu]);
      }
    }
  }

  // Rearrange to find \Gamma_{lam nu mu}
  for (int lam = 0; lam < NDIM; lam++) {
    for (int nu = 0; nu < NDIM; nu++) {
      for (int mu = 0; mu < NDIM; mu++) {
        tmp[lam][nu][mu] = 0.5 * (G->conn[nu][lam][mu][j][i] + 
                                  G->conn[mu][lam][nu][j][i] - 
                                  G->conn[mu][nu][lam][j][i]);
      }
    }
  }

  // now mu nu kap

  // Raise index to get \Gamma^lam_{nu mu}
  for (int lam = 0; lam < NDIM; lam++) {
    for (int nu = 0; nu < NDIM; nu++) {
      for (int mu = 0; mu < NDIM; mu++) {
        G->conn[lam][nu][mu][j][i] = 0.;
        for (int kap = 0; kap < NDIM; kap++)
          G->conn[lam][nu][mu][j][i] += G->gcon[CENT][lam][kap][j][i]*
                                        tmp[kap][nu][mu];
      }
    }
  }
}

// Lower a contravariant rank-1 tensor to a covariant one
inline void lower_grid(GridVector vcon, GridVector vcov, struct GridGeom *G, int i,
  int j, int loc)
{
  for (int mu = 0; mu < NDIM; mu++) {
    vcov[mu][j][i] = 0.;
    for (int nu = 0; nu < NDIM; nu++) {
      vcov[mu][j][i] += G->gcov[loc][mu][nu][j][i]*vcon[nu][j][i];
    }
  }
}

// Lower the grid of contravariant rank-1 tensors to covariant ones
void lower_grid_vec(GridVector vcon, GridVector vcov, struct GridGeom *G, int jstart, int jstop, int istart, int istop, int loc)
{
#pragma omp parallel for simd collapse(3)
  DLOOP1 {
    ZSLOOP(jstart, jstop, istart, istop) vcov[mu][j][i] = 0.;
  }
#pragma omp parallel for simd collapse(4)
  DLOOP2 {
      ZSLOOP(jstart, jstop, istart, istop) vcov[mu][j][i] += G->gcov[loc][mu][nu][j][i]*vcon[nu][j][i];
  }
}

inline void raise_grid(GridVector vcov, GridVector vcon, struct GridGeom *G, int i, int j, int loc)
{
  for (int mu = 0; mu < NDIM; mu++) {
    vcon[mu][j][i] = 0.;
    for (int nu = 0; nu < NDIM; nu++) {
      vcon[mu][j][i] += G->gcon[loc][mu][nu][j][i]*vcov[nu][j][i];
    }
  }
}

// Take dot product of a contravariant and covariant rank-1 tensor
inline double dot(double vcon[NDIM], double vcov[NDIM])
{
  double dot = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    dot += vcon[mu]*vcov[mu];
  }
  return dot;
}

// Minor of a 4x4 matrix
inline double MINOR(double m[16], int r0, int r1, int r2, int c0, int c1, int c2)
{
  return m[4*r0+c0]*(m[4*r1+c1]*m[4*r2+c2] - m[4*r2+c1]*m[4*r1+c2]) -
         m[4*r0+c1]*(m[4*r1+c0]*m[4*r2+c2] - m[4*r2+c0]*m[4*r1+c2]) +
         m[4*r0+c2]*(m[4*r1+c0]*m[4*r2+c1] - m[4*r2+c0]*m[4*r1+c1]);
}

inline void adjoint(double m[16], double adjOut[16])
{
  adjOut[ 0] =  MINOR(m,1,2,3,1,2,3);
  adjOut[ 1] = -MINOR(m,0,2,3,1,2,3);
  adjOut[ 2] =  MINOR(m,0,1,3,1,2,3);
  adjOut[ 3] = -MINOR(m,0,1,2,1,2,3);

  adjOut[ 4] = -MINOR(m,1,2,3,0,2,3);
  adjOut[ 5] =  MINOR(m,0,2,3,0,2,3);
  adjOut[ 6] = -MINOR(m,0,1,3,0,2,3);
  adjOut[ 7] =  MINOR(m,0,1,2,0,2,3);

  adjOut[ 8] =  MINOR(m,1,2,3,0,1,3);
  adjOut[ 9] = -MINOR(m,0,2,3,0,1,3);
  adjOut[10] =  MINOR(m,0,1,3,0,1,3);
  adjOut[11] = -MINOR(m,0,1,2,0,1,3);

  adjOut[12] = -MINOR(m,1,2,3,0,1,2);
  adjOut[13] =  MINOR(m,0,2,3,0,1,2);
  adjOut[14] = -MINOR(m,0,1,3,0,1,2);
  adjOut[15] =  MINOR(m,0,1,2,0,1,2);
}

// Computes determinant of 4x4 tensor
inline double determinant(double m[16])
{
  return m[0]*MINOR(m,1,2,3,1,2,3) -
         m[1]*MINOR(m,1,2,3,0,2,3) +
         m[2]*MINOR(m,1,2,3,0,1,3) -
         m[3]*MINOR(m,1,2,3,0,1,2);
}

// Computes inverse of a 4x4 matrix
inline double invert(double *m, double *invOut)
{
  adjoint(m, invOut);

  double det = determinant(m);
  double inv_det = 1. / det;
  for (int i = 0; i < 16; ++i) {
    invOut[i] = invOut[i]*inv_det;
  }

  return det;
}
