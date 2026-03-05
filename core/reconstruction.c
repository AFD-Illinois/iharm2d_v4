/**
 * @file reconstruction.c
 * @brief Spatial reconstruction algorithms for obtaining left/right interface states.
 *
 * @details Implements three reconstruction schemes selectable at compile time via the
 * RECONSTRUCTION macro in parameters.h:
 *
 * - **LINEAR** (linear_mc()): Piecewise-linear reconstruction with the Monotonized
 *   Central (MC) slope limiter.  Second-order accurate, TVD.
 *
 * - **WENO** (weno()): 5th-order Weighted Essentially Non-Oscillatory reconstruction
 *   based on Tchekhovskoy et al. (2007) and Shu (2011).  Uses three candidate
 *   stencils and smoothness indicators to weight them optimally.
 *
 * - **MP5** (mp5()): 5th-order Monotonicity Preserving reconstruction imported from
 *   the PLUTO code (Suresh & Huynh 1997).  Uses a flux-limiting approach to preserve
 *   monotonicity while achieving high-order accuracy away from extrema.
 *
 * The main entry point reconstruct() dispatches to the selected algorithm via the
 * preprocessor macro RECON_ALGO and is parallelized with OpenMP over primitive variables
 * and zones.  All schemes require a 5-point stencil, needing NG >= 3 ghost zones.
 *
 * @note WENO and MP5 require NG >= 3.  A compile-time check enforces this.
 */

/*---------------------------------------------------------------------------------

  RECONSTRUCTION.C

  -Linear, WENO and MP5 reconstruction algorithms

---------------------------------------------------------------------------------*/

#include "decs.h"

// Set macro RECON_ALGO
#if RECONSTRUCTION == LINEAR
#define RECON_ALGO linear_mc
#elif RECONSTRUCTION == WENO
#define RECON_ALGO weno
#elif RECONSTRUCTION == MP5
#define RECON_ALGO mp5
#else
#error "Reconstruction not specified!"
#endif

// Sanity checks
#if (RECONSTRUCTION == WENO || RECONSTRUCTION == MP5) && NG < 3
#error "not enough ghost zones! PPM/WENO/MP5 + NG < 3\n"
#endif

void linear_mc(double unused1, double x1, double x2, double x3, double unused2, double *lout, double *rout);
void weno(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout);

double median(double a, double b, double c);
double mp5_subcalc(double Fjm2, double Fjm1, double Fj, double Fjp1, double Fjp2);
void mp5(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout);

/**
 * @brief Linear reconstruction with Monotonized Central (MC) slope limiter.
 *
 * @details Computes the limited slope and reconstructs the left and right interface
 * values from cell-centered values x1, x2, x3 (at i-1, i, i+1):
 * @f{align*}{
 *   \Delta q_m &= 2(q_i - q_{i-1}),\quad
 *   \Delta q_p = 2(q_{i+1} - q_i),\quad
 *   \Delta q_c = \frac{1}{2}(q_{i+1} - q_{i-1})\\
 *   s &= \text{MC}(\Delta q_m, \Delta q_p, \Delta q_c)\\
 *   q_L &= q_i - \tfrac{1}{2}s,\quad q_R = q_i + \tfrac{1}{2}s
 * @f}
 *
 * @param unused1  Unused (padding for uniform 5-point API with WENO/MP5).
 * @param x1  Value at i-1.
 * @param x2  Value at i (cell center being reconstructed).
 * @param x3  Value at i+1.
 * @param unused2  Unused.
 * @param[out] lout  Left face value @f$q_{i-1/2,L}@f$.
 * @param[out] rout  Right face value @f$q_{i+1/2,R}@f$.
 */
inline void linear_mc(double unused1, double x1, double x2, double x3, double unused2, double *lout, double *rout)
{
  double Dqm,Dqp,Dqc,s;

  Dqm = 2. * (x2 - x1);
  Dqp = 2. * (x3 - x2);
  Dqc = 0.5 * (x3 - x1);

  s = Dqm * Dqp;

  if (s <= 0.)
    s = 0.;
  else {
    if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
      s = Dqm;
    else if (fabs(Dqp) < fabs(Dqc))
      s = Dqp;
    else
      s = Dqc;
  }

  // Reconstruct left, right
  *lout = x2 - 0.5*s;
  *rout = x2 + 0.5*s;
}

/**
 * @brief 5th-order WENO reconstruction (Tchekhovskoy et al. 2007, Shu 2011).
 *
 * @details Uses three candidate 3-point stencils to reconstruct both left and
 * right interface values.  The nonlinear weights are computed from smoothness
 * indicators (beta) so that the scheme reduces to the optimal 5th-order stencil
 * in smooth regions and degrades gracefully near discontinuities.
 *
 * @param x1  Value at i-2.
 * @param x2  Value at i-1.
 * @param x3  Value at i.
 * @param x4  Value at i+1.
 * @param x5  Value at i+2.
 * @param[out] lout  Left interface value @f$q_{i-1/2,L}@f$.
 * @param[out] rout  Right interface value @f$q_{i+1/2,R}@f$.
 */
inline void weno(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout)
{
  // S11 1, 2, 3
  double vr[3], vl[3];
  vr[0] =  (3./8.)*x1 - (5./4.)*x2 + (15./8.)*x3;
  vr[1] = (-1./8.)*x2 + (3./4.)*x3 + (3./8.)*x4;
  vr[2] =  (3./8.)*x3 + (3./4.)*x4 - (1./8.)*x5;

  vl[0] =  (3./8.)*x5 - (5./4.)*x4 + (15./8.)*x3;
  vl[1] = (-1./8.)*x4 + (3./4.)*x3 + (3./8.)*x2;
  vl[2] =  (3./8.)*x3 + (3./4.)*x2 - (1./8.)*x1;

  // Smoothness indicators, T07 A18 or S11 8
  double beta[3];
  beta[0] = (13./12.)*pow(x1 - 2.*x2 + x3, 2) +
            (1./4.)*pow(x1 - 4.*x2 + 3.*x3, 2);
  beta[1] = (13./12.)*pow(x2 - 2.*x3 + x4, 2) +
            (1./4.)*pow(x4 - x2, 2);
  beta[2] = (13./12.)*pow(x3 - 2.*x4 + x5, 2) +
            (1./4.)*pow(x5 - 4.*x4 + 3.*x3, 2);

  // Nonlinear weights S11 9
  double den, wtr[3], Wr, wr[3], wtl[3], Wl, wl[3], eps;
  eps=1.e-26;

  den = eps + beta[0]; den *= den; wtr[0] = (1./16.)/den;
  den = eps + beta[1]; den *= den; wtr[1] = (5./8. )/den;
  den = eps + beta[2]; den *= den; wtr[2] = (5./16.)/den;
  Wr = wtr[0] + wtr[1] + wtr[2];
  wr[0] = wtr[0]/Wr ;
  wr[1] = wtr[1]/Wr ;
  wr[2] = wtr[2]/Wr ;

  den = eps + beta[2]; den *= den; wtl[0] = (1./16.)/den;
  den = eps + beta[1]; den *= den; wtl[1] = (5./8. )/den;
  den = eps + beta[0]; den *= den; wtl[2] = (5./16.)/den;
  Wl = wtl[0] + wtl[1] + wtl[2];
  wl[0] = wtl[0]/Wl;
  wl[1] = wtl[1]/Wl;
  wl[2] = wtl[2]/Wl;

  *lout = vl[0]*wl[0] + vl[1]*wl[1] + vl[2]*wl[2];
  *rout = vr[0]*wr[0] + vr[1]*wr[1] + vr[2]*wr[2];
}

// MP5 reconstruction from PLUTO
// Imported by Mani Chandra
#define MINMOD(a, b) ((a)*(b) > 0.0 ? (fabs(a) < fabs(b) ? (a):(b)):0.0)
inline double median(double a, double b, double c)
{
  return (a + MINMOD(b - a, c - a));
}
#define ALPHA (4.0)
#define EPSM (1.e-12)
inline double mp5_subcalc(double Fjm2, double Fjm1, double Fj, double Fjp1, double Fjp2)
{
  double f, d2, d2p, d2m;
  double dMMm, dMMp;
  double scrh1,scrh2, Fmin, Fmax;
  double fAV, fMD, fLC, fUL, fMP;

  f  = 2.0*Fjm2 - 13.0*Fjm1 + 47.0*Fj + 27.0*Fjp1 - 3.0*Fjp2;
  f /= 60.0;

  fMP = Fj + MINMOD(Fjp1 - Fj, ALPHA*(Fj - Fjm1));

  if ((f - Fj)*(f - fMP) <= EPSM)
    return f;

  d2m = Fjm2 + Fj   - 2.0*Fjm1;              // Eqn. 2.19
  d2  = Fjm1 + Fjp1 - 2.0*Fj;
  d2p = Fj   + Fjp2 - 2.0*Fjp1;              // Eqn. 2.19

  scrh1 = MINMOD(4.0*d2 - d2p, 4.0*d2p - d2);
  scrh2 = MINMOD(d2, d2p);
  dMMp  = MINMOD(scrh1,scrh2);               // Eqn. 2.27

  scrh1 = MINMOD(4.0*d2m - d2, 4.0*d2 - d2m);
  scrh2 = MINMOD(d2, d2m);
  dMMm  = MINMOD(scrh1,scrh2);               // Eqn. 2.27

  fUL = Fj + ALPHA*(Fj - Fjm1);              // Eqn. 2.8
  fAV = 0.5*(Fj + Fjp1);                     // Eqn. 2.16
  fMD = fAV - 0.5*dMMp;                      // Eqn. 2.28
  fLC = 0.5*(3.0*Fj - Fjm1) + 4.0/3.0*dMMm;  // Eqn. 2.29

  scrh1 = fmin(Fj, Fjp1); scrh1 = fmin(scrh1, fMD);
  scrh2 = fmin(Fj, fUL);    scrh2 = fmin(scrh2, fLC);
  Fmin  = fmax(scrh1, scrh2);                // Eqn. (2.24a)

  scrh1 = fmax(Fj, Fjp1); scrh1 = fmax(scrh1, fMD);
  scrh2 = fmax(Fj, fUL);    scrh2 = fmax(scrh2, fLC);
  Fmax  = fmin(scrh1, scrh2);                // Eqn. 2.24b

  f = median(f, Fmin, Fmax);                 // Eqn. 2.26
  return f;
}

/**
 * @brief MP5 reconstruction from PLUTO (Suresh & Huynh 1997). Imported by Mani Chandra.
 * @param x1  Value at i-2.
 * @param x2  Value at i-1.
 * @param x3  Value at i.
 * @param x4  Value at i+1.
 * @param x5  Value at i+2.
 * @param[out] lout  Left interface value.
 * @param[out] rout  Right interface value.
 */
inline void mp5(double x1, double x2, double x3, double x4, double x5, double *lout,
  double *rout)
{
  *rout = mp5_subcalc(x1, x2, x3, x4, x5);
  *lout = mp5_subcalc(x5, x4, x3, x2, x1);
}
#undef MINMOD

/**
 * @brief Reconstruct left and right primitive states at cell interfaces.
 *
 * @details Iterates over all primitive variables and zones (using OpenMP), applying
 * the compile-time-selected reconstruction algorithm (RECON_ALGO) in the given
 * coordinate direction to produce face-centered left (Pl) and right (Pr) states.
 * For direction 1 (X1), uses the stencil [i-2, i-1, i, i+1, i+2].
 * For direction 2 (X2), uses the stencil [j-2, j-1, j, j+1, j+2].
 *
 * @param S    Source fluid state (its P array provides the cell-centered values).
 * @param Pl   Output left-reconstructed primitives at cell interfaces.
 * @param Pr   Output right-reconstructed primitives at cell interfaces.
 * @param dir  Direction of reconstruction: 1 (X1) or 2 (X2).
 */
void reconstruct(struct FluidState *S, GridPrim Pl, GridPrim Pr, int dir)
{
  timer_start(TIMER_RECON);
  if (dir == 1)
#pragma omp parallel for collapse(2)
    PLOOP 
      JSLOOP(-1, N2)
        ISLOOP(-1, N1)
          RECON_ALGO(S->P[ip][j][i-2], S->P[ip][j][i-1], S->P[ip][j][i], S->P[ip][j][i+1], S->P[ip][j][i+2], &(Pl[ip][j][i]), &(Pr[ip][j][i]));
  else if (dir == 2)
#pragma omp parallel for collapse(2)
    PLOOP
      JSLOOP(-1, N2)
        ISLOOP(-1, N1)
          RECON_ALGO(S->P[ip][j-2][i], S->P[ip][j-1][i], S->P[ip][j][i], S->P[ip][j+1][i], S->P[ip][j+2][i], &(Pl[ip][j][i]), &(Pr[ip][j][i]));
  
  timer_stop(TIMER_RECON);
}
