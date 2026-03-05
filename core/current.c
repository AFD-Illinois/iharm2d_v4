/**
 * @file current.c
 * @brief Compute the electromagnetic 4-current j^μ from the fluid state.
 *
 * @details Evaluates the 4-current density by taking finite differences of
 * the contravariant Maxwell tensor @f$\sqrt{-g}\,F^{\mu\nu}@f$:
 * @f[
 *   4\pi\,j^\mu = \frac{1}{\sqrt{-g}}
 *     \left(\partial_0 (\sqrt{-g} F^{0\mu}) + \partial_i (\sqrt{-g} F^{i\mu})\right)
 * @f]
 * The time derivative uses the two saved primitive states (current and previous
 * step), i.e. the result is centred at the half-timestep.  The spatial
 * derivatives use second-order centred differences of the time-centred
 * (averaged) state.
 *
 * Helper routines:
 * - **gFcon_calc()**: Computes @f$\sqrt{-g}\,F^{\mu\nu}@f$ from the ideal MHD
 *   expression @f$F^{\mu\nu} = -\epsilon^{\mu\nu\kappa\lambda}\,u_\kappa b_\lambda / g@f$.
 * - **antisym()**: Returns the value of the completely antisymmetric Levi-Civita
 *   symbol @f$\epsilon^{abcd}@f$ as +1, -1, or 0.
 * - **pp()**: Determines the parity of a permutation (Hardy's algorithm).
 *
 * @note The 4-current output is stored in S->jcon[mu][j][i] and written to
 * dump files for post-processing (e.g., computing resistive dissipation).
 */

/*---------------------------------------------------------------------------------

  CURRENT.C

  -Calculate 4-current from fluid variables

-----------------------------------------------------------------------------------*/

#include "decs.h"

double gFcon_calc(struct GridGeom *G, struct FluidState *S, int mu, int nu, int i, int j);
int antisym(int a, int b, int c, int d);
int pp(int n, int *P);


static struct FluidState *Sa;

/**
 * @brief Compute the contravariant 4-current j^μ over all active zones.
 *
 * @details Uses a staggered time-centring scheme:
 * - Time derivative terms use S (end of step) and Ssave (beginning of step).
 * - Spatial derivative terms use Sa = (S + Ssave)/2 (time-centred average).
 * get_state_vec() is called on all three states before the differentiation loop
 * so that ucov and bcov are consistent.
 * The final j^μ at each zone is:
 * @f[
 *   j^\mu = \frac{1}{\sqrt{4\pi}\,\sqrt{-g}}
 *     \left[\frac{\sqrt{-g}F^{0\mu}_{\rm new} - \sqrt{-g}F^{0\mu}_{\rm old}}{\Delta t}
 *         + \frac{\sqrt{-g}F^{1\mu}(i+1) - \sqrt{-g}F^{1\mu}(i-1)}{2\Delta X^1}
 *         + \frac{\sqrt{-g}F^{2\mu}(j+1) - \sqrt{-g}F^{2\mu}(j-1)}{2\Delta X^2}
 *     \right]
 * @f]
 *
 * @param G       Grid geometry.
 * @param S       Fluid state at current time (end of step); jcon is written here.
 * @param Ssave   Fluid state at previous time (beginning of step).
 * @param dtsave  Physical time step Δt used to form the time derivative.
 */
// Calculate the current
void current_calc(struct GridGeom *G, struct FluidState *S, struct FluidState *Ssave, double dtsave)
{
  timer_start(TIMER_CURRENT);

  static int first_run = 1;
  if (first_run)
  {
    // We only need the primitives, but this is fast
    Sa = calloc(1,sizeof(struct FluidState));
    first_run = 0;
  }

  // Calculate time-centered P value of the primitives
#pragma omp parallel for simd collapse(2)
  PLOOP
  {
    ZLOOPALL
    {
      Sa->P[ip][j][i] = 0.5*(S->P[ip][j][i] + Ssave->P[ip][j][i]);
    }
  }

  // Keep all get_state calls outside the loop so it doesn't modify S{a,save}
  get_state_vec(G, S, CENT, -1, N2, -1, N1);
  get_state_vec(G, Ssave, CENT, -1, N2, -1, N1);
  get_state_vec(G, Sa, CENT, -1, N2, -1, N1);

  // Initialize the 4-current to zero
#pragma omp parallel for simd collapse(2)
  DLOOP1 ZLOOPALL S->jcon[mu][j][i] = 0.;

  // Calculate j^{\mu} using centered differences for active zones
#pragma omp parallel for collapse(2)
  ZLOOP 
  {
    double gF0p[NDIM], gF0m[NDIM], gF1p[NDIM], gF1m[NDIM], gF2p[NDIM], gF2m[NDIM];

    // Get sqrt{-g}*F^{mu nu} at neighboring points

    // X0
    DLOOP1 
    {
      gF0p[mu] = gFcon_calc(G, S,  0, mu, i, j);
      gF0m[mu] = gFcon_calc(G, Ssave, 0, mu, i, j);
    }

    // X1
    DLOOP1
    {
      gF1p[mu] = gFcon_calc(G, Sa, 1, mu, i+1, j);
      gF1m[mu] = gFcon_calc(G, Sa, 1, mu, i-1, j);
    }

    // X2
    DLOOP1
    {
      gF2p[mu] = gFcon_calc(G, Sa, 2, mu, i, j+1);
      gF2m[mu] = gFcon_calc(G, Sa, 2, mu, i, j-1);
    }

    // Difference: D_mu f^{mu nu} = 4 \pi j^nu
    DLOOP1 
    {
      // extra factor of sqrt(4*pi)*j given harm's b_unit
      S->jcon[mu][j][i] = (1./(sqrt(4.*M_PI)*G->gdet[CENT][j][i]))*(
                           (gF0p[mu] - gF0m[mu])/dtsave +
                           (gF1p[mu] - gF1m[mu])/(2.*dx[1]) +
                           (gF2p[mu] - gF2m[mu])/(2.*dx[2]));
    }
  }
  timer_stop(TIMER_CURRENT);
}

/**
 * @brief Compute @f$\sqrt{-g}\,F^{\mu\nu}@f$ at zone (i, j).
 *
 * @details Uses the ideal MHD relation @f$F^{\mu\nu} = b^\mu u^\nu - b^\nu u^\mu@f$
 * (equivalent to the Levi-Civita form) via:
 * @f[
 *   \sqrt{-g}\,F^{\mu\nu} = -\frac{1}{\sqrt{-g}}\sum_{\kappa,\lambda}
 *     \epsilon^{\mu\nu\kappa\lambda}\,u_\kappa\,b_\lambda \cdot \sqrt{-g}
 * @f]
 * Returns 0 immediately if mu == nu (antisymmetry).
 *
 * @param G   Grid geometry (gdet at CENT used for normalisation).
 * @param S   Fluid state (ucov, bcov must be pre-calculated via get_state()).
 * @param mu  First index.
 * @param nu  Second index.
 * @param i   X1 zone index.
 * @param j   X2 zone index.
 * @return    @f$\sqrt{-g}\,F^{\mu\nu}@f$ at zone (i, j).
 */
// Return mu, nu component of contravariant Maxwell tensor at grid zone i, j multiplied by gdet
inline double gFcon_calc(struct GridGeom *G, struct FluidState *S, int mu, int nu, int i, int j)
{
  double Fcon;

  if (mu == nu) return 0.;

  Fcon = 0.;
  for (int kap = 0; kap < NDIM; kap++)
  {
    for (int lam = 0; lam < NDIM; lam++)
    {
      Fcon += (-1./G->gdet[CENT][j][i])*antisym(mu,nu,kap,lam)*S->ucov[kap][j][i]*S->bcov[lam][j][i];
    }
  }

  return Fcon*G->gdet[CENT][j][i];
}

/**
 * @brief Evaluate the 4D Levi-Civita symbol @f$\epsilon^{abcd}@f$.
 *
 * @details Returns +1 if (a,b,c,d) is an even permutation of (0,1,2,3),
 * -1 if an odd permutation, and 0 if any two indices are equal.
 * Returns 100 if any index is out of [0,3].
 *
 * @param a  First index (0–3).
 * @param b  Second index.
 * @param c  Third index.
 * @param d  Fourth index.
 * @return +1, -1, 0, or 100 (invalid).
 */
// Completely antisymmetric 4D symbol
inline int antisym(int a, int b, int c, int d)
{
  // Check for valid permutation
  if (a < 0 || a > 3) return 100;
  if (b < 0 || b > 3) return 100;
  if (c < 0 || c > 3) return 100;
  if (d < 0 || d > 3) return 100;

  // Entries different? 
  if (a == b) return 0;
  if (a == c) return 0;
  if (a == d) return 0;
  if (b == c) return 0;
  if (b == d) return 0;
  if (c == d) return 0;

  // Determine parity of permutation
  int p[4] = {a, b, c, d};

  return pp(4, p);
}

// Due to Norm Hardy; good for general n
inline int pp(int n, int P[n])
{
  int x;
  int p = 0;
  int v[n];

  for (int j = 0; j < n; j++) v[j] = 0;

  for (int j = 0; j < n; j++) {
    if (v[j]) p++; 
    else 
    {
      x = j;
      do 
      {
        x = P[x];
        v[x] = 1;
      } while (x != j);
    }
  }

  if (p % 2 == 0) return 1;
  else return -1;
}
