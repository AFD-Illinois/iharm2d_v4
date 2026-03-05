/**
 * @file metric.c
 * @brief Metric tensor utilities: inversion, Christoffel symbols, index raising/lowering.
 *
 * @details Provides:
 * - **gcon_func()**: Invert the 4x4 covariant metric @f$g_{\mu\nu}@f$ to obtain
 *   @f$g^{\mu\nu}@f$ using the adjoint / cofactor method.  Returns @f$\sqrt{|g|}@f$.
 * - **conn_func()**: Compute all 40 independent Christoffel connection coefficients
 *   @f$\Gamma^\lambda_{\mu\nu}@f$ at a zone center via numerical finite differencing of
 *   the metric (using @c DELTA = 1e-5 as step size).
 * - **lower_grid() / raise_grid()**: Index lowering/raising over the full GridVector array.
 * - **dot()**: Scalar product @f$v^\mu w_\mu@f$ of a contravariant and covariant 4-vector.
 * - **invert()**: General 4x4 matrix inversion used by gcon_func() and coord.c.
 *
 * @note Connection coefficients are computed once at startup and cached in GridGeom::conn.
 */

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

/**
 * @brief Invert the 4x4 covariant metric tensor to obtain the contravariant metric.
 *
 * @details Uses the adjoint (cofactor matrix) method.  Also computes
 * @f$\sqrt{|det(g_{\mu\nu})|}@f$ and stores it via the return value.
 *
 * @param gcov   4x4 covariant metric input.
 * @param gcon   4x4 contravariant metric output.
 * @return       @f$\sqrt{|det(g)|}@f$ (equals gdet).
 */
inline double gcon_func(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM])
{
  double gdet = invert(&gcov[0][0],&gcon[0][0]);
  return sqrt(fabs(gdet));
}

// Make a local copy of gcov (from Grid struct)
/**
 * @brief Copy covariant metric g_mu_nu from GridGeom into a local 4x4 array.
 * @param G    Grid geometry.
 * @param i    X1 zone index.
 * @param j    X2 zone index.
 * @param loc  Grid centering location.
 * @param gcov Output 4x4 covariant metric.
 */
inline void get_gcov(struct GridGeom *G, int i, int j, int loc, double gcov[NDIM][NDIM]) {
  DLOOP2 gcov[mu][nu] = G->gcov[loc][mu][nu][j][i];
}

// Make a local copy of gcon (from Grid struct)
/**
 * @brief Copy contravariant metric g^mu^nu from GridGeom into a local 4x4 array.
 * @param G    Grid geometry.
 * @param i    X1 zone index.
 * @param j    X2 zone index.
 * @param loc  Grid centering location.
 * @param gcon Output 4x4 contravariant metric.
 */
inline void get_gcon(struct GridGeom *G, int i, int j, int loc, double gcon[NDIM][NDIM])
{
  DLOOP2 gcon[mu][nu] = G->gcon[loc][mu][nu][j][i];
}

// Calculate connection coefficient
/**
 * @brief Compute all Christoffel connection coefficients @f$\Gamma^\lambda_{\mu\nu}@f$ at zone (i, j).
 *
 * @details Uses second-order centered finite differences of the metric:
 * @f[
 *   \partial_\mu g_{\lambda\nu} \approx \frac{g_{\lambda\nu}(X + \delta_\mu) - g_{\lambda\nu}(X - \delta_\mu)}{2\,\text{DELTA}}
 * @f]
 * then combines the first derivatives into the standard expression for the Christoffel symbols
 * and raises the first index using @f$g^{\mu\nu}@f$ to obtain @f$\Gamma^\lambda_{\mu\nu}@f$.
 * Results are stored in G->conn[lam][nu][mu][j][i].
 *
 * @param G  Grid geometry (gcon at CENT must be pre-computed; conn is written here).
 * @param i  X1 zone index.
 * @param j  X2 zone index.
 */
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
/**
 * @brief Lower vector index: vcov_mu = g_mu_nu * vcon^nu at zone (i, j).
 * @param vcon  Contravariant 4-vector.
 * @param vcov  Covariant 4-vector output.
 * @param G     Grid geometry.
 * @param i     X1 zone index.
 * @param j     X2 zone index.
 * @param loc   Grid centering location.
 */
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
/** @brief Vectorized lower_grid() over a rectangular zone range (OpenMP SIMD). */
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

// Raise the grid of covariant rank-1 tensors to covariant ones
/**
 * @brief Raise vector index: vcon^mu = g^mu^nu * vcov_nu at zone (i, j).
 * @param vcov  Covariant 4-vector input.
 * @param vcon  Contravariant 4-vector output.
 * @param G     Grid geometry.
 * @param i     X1 zone index.
 * @param j     X2 zone index.
 * @param loc   Grid centering location.
 */
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
/**
 * @brief Compute the scalar product @f$v^\mu w_\mu@f$.
 * @param vcon  Contravariant 4-vector.
 * @param vcov  Covariant 4-vector.
 * @return @f$\sum_\mu v^\mu w_\mu@f$.
 */
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

// Adjoint of a 4x4 matrix
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
/**
 * @brief Compute the inverse of a 4x4 matrix using the adjoint method.
 *
 * @details Evaluates @f$A^{-1} = \text{adj}(A) / \det(A)@f$.
 * If |det(A)| < 1e-10 the inversion is singular; return value is the determinant
 * and entries may be garbage.
 *
 * @param m       Input matrix, stored as a flat 16-element (row-major) array.
 * @param invOut  Output inverse matrix, same storage convention.
 * @return        Determinant of @c m.
 */
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
