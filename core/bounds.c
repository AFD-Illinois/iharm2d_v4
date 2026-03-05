/**
 * @file bounds.c
 * @brief Physical boundary conditions for the fluid primitives and fluxes.
 *
 * @details Implements set_bounds() which fills ghost zone primitive variables
 * using one of four strategies per boundary (configured in parameters.h):
 *
 * - **OUTFLOW**: Copy the outermost active zone into all ghost zones. For MKS,
 *   rescale the magnetic field by the ratio of sqrt(-g) values to approximately
 *   conserve magnetic flux.
 * - **PERIODIC**: Map ghost zones to the corresponding real zone on the other side.
 * - **POLAR** (X2 only): Reflect the fluid state across the axis.  U2 and B2
 *   change sign because they are the θ-direction components of their fields:
 *   the θ-direction reverses under the reflection θ → −θ, so these components
 *   are odd functions of θ (they vanish at the axis).  The r- and φ-components
 *   (U1, U3, B1, B3) are even in θ and do not change sign.
 * - **USER** (X1R only): Delegate to the problem-specific function bound_gas_prob_x1r().
 *
 * For MKS the radial boundaries additionally enforce:
 * - **Inflow check** (inflow_check()): If X1L_INFLOW=0 (X1R_INFLOW=0),
 *   the ghost-zone velocity is adjusted so that u^r = 0 when the fluid would flow
 *   inward (outward) through the inner (outer) boundary.
 * - **Flux fix** (fix_flux()): Called by advance_fluid() before the conservative
 *   update. Zeros the X2 fluxes at the poles, reflects B2 flux across them,
 *   and clips the mass flux to enforce the inflow/outflow condition.
 *
 * @note set_bounds() is called after every U_to_P() inversion and after fixup()
 * during both the predictor and corrector substeps.
 */

/*---------------------------------------------------------------------------------

  BOUNDS.C

  -Implements physical boundary conditions
  -Ensure no inflow at radial boundaries
  -Ensure radial mass flux at radial boundaries is zero
  -B2 flux at X1 and X3 faces at polar boundaries is reflected for ghost zones
  -All X2 fluxes at polar boundaries are zeroed

---------------------------------------------------------------------------------*/

#include "decs.h"

// Sanity checks: grid dimensions, supported boundary conditions
#if N2 > 1 && N2 < NG
#error "N2 must be >= NG"
#endif


#if X1L_BOUND != PERIODIC && X1L_BOUND != OUTFLOW
#error "Unsupported X1L_BOUND"
#endif
#if X1R_BOUND != PERIODIC && X1R_BOUND != OUTFLOW && X1R_BOUND != USER
#error "Unsupported X1R_BOUND"
#endif

#if X2L_BOUND != PERIODIC && X2L_BOUND != OUTFLOW && X2L_BOUND != POLAR
#error "Unsupported X2L_BOUND"
#endif
#if X2R_BOUND != PERIODIC && X2R_BOUND != OUTFLOW && X2R_BOUND != POLAR
#error "Unsupported X2R_BOUND"
#endif

void inflow_check(struct GridGeom *G, struct FluidState *S, int i, int j, int type);

// Apply boundary conditions along X1 and X2
/**
 * @brief Apply all configured boundary conditions to the primitive variable array.
 *
 * @details Fills ghost zones on all four boundaries of the 2D grid in order:
 * X1 left, X1 right, X2 left, X2 right.  For each boundary, the behavior is
 * determined by the compile-time flags X1L_BOUND, X1R_BOUND, X2L_BOUND, X2R_BOUND.
 * An optional inflow-check pass (MKS only) ensures no fluid enters through radial
 * boundaries if X1L_INFLOW or X1R_INFLOW is set to 0.
 *
 * @param G  Grid geometry.
 * @param S  Fluid state (P[] is read and ghost zones are written).
 */
void set_bounds(struct GridGeom *G, struct FluidState *S)
{
  timer_start(TIMER_BOUND);


#pragma omp parallel for
  JLOOP 
  {
    ISLOOP(-NG, -1)
    {
#if N1 < NG
      int iactive = NG;
      PLOOP S->P[ip][j][i] = S->P[ip][j][iactive];
      pflag[j][i] = pflag[j][iactive];
#elif X1L_BOUND == OUTFLOW
			int iz = 0 + NG;
			PLOOP S->P[ip][j][i] = S->P[ip][j][iz];
			pflag[j][i] = pflag[j][iz];

			double rescale = G->gdet[CENT][j][iz]/G->gdet[CENT][j][i];
			S->P[B1][j][i] *= rescale;
			S->P[B2][j][i] *= rescale;
			S->P[B3][j][i] *= rescale;
#elif X1L_BOUND == PERIODIC
			int iz = N1 + i;
			PLOOP S->P[ip][j][i] = S->P[ip][j][iz];
			pflag[j][i] = pflag[j][iz];
#endif
    }
  }

#if METRIC == MKS
  if(X1L_INFLOW == 0)
  {
    // Make sure there is no inflow at the inner boundary

#pragma omp parallel for
    JLOOP
      ISLOOP(-NG, -1)
        inflow_check(G, S, i, j, 0);
  }
#endif


#pragma omp parallel for
  JLOOP
  {
    ISLOOP(N1, N1 - 1 + NG)
    {
#if N1 < NG
      int iactive = N1 - 1 + NG;
      PLOOP S->P[ip][j][i] = S->P[ip][j][iactive];
      pflag[j][i] = pflag[j][iactive];
#elif X1R_BOUND == OUTFLOW
      int iz = N1 - 1 + NG;
      PLOOP S->P[ip][j][i] = S->P[ip][j][iz];
      pflag[j][i] = pflag[j][iz];

      double rescale = G->gdet[CENT][j][iz]/G->gdet[CENT][j][i];
      S->P[B1][j][i] *= rescale;
      S->P[B2][j][i] *= rescale;
      S->P[B3][j][i] *= rescale;
#elif X1R_BOUND == USER
      bound_gas_prob_x1r(i, j, S->P, G);
#elif X1R_BOUND == PERIODIC
			int iz = i - N1;
			PLOOP S->P[ip][j][i] = S->P[ip][j][iz];
			pflag[j][i] = pflag[j][iz];
#endif
    }
  }

#if METRIC == MKS
  if(X1R_INFLOW == 0)
  {
    // Make sure there is no inflow at the outer boundary

#pragma omp parallel for
    JLOOP
      ISLOOP(N1, N1 - 1 + NG)
        inflow_check(G, S, i, j, 1);
  }
#endif


#pragma omp parallel for
  ILOOPALL
  {
    JSLOOP(-NG, -1)
    {
#if N2 < NG
      int jactive = NG;
      PLOOP S->P[ip][j][i] = S->P[ip][jactive][i];
      pflag[j][i] = pflag[jactive][i];
#elif X2L_BOUND == OUTFLOW
      int jz = 0 + NG ;
      PLOOP S->P[ip][j][i] = S->P[ip][jz][i];
      pflag[j][i] = pflag[jz][i];
#elif X2L_BOUND == POLAR
      // Reflect the zone past NG by NG-j
      int jrefl = NG + (NG - j) - 1;
      PLOOP S->P[ip][j][i] = S->P[ip][jrefl][i];
      pflag[j][i] = pflag[jrefl][i];
      S->P[U2][j][i] *= -1.;
      S->P[B2][j][i] *= -1.;
#elif X2L_BOUND == PERIODIC
			int jz = N2 + j;
			PLOOP S->P[ip][j][i] = S->P[ip][jz][i];
			pflag[j][i] = pflag[jz][i];
#endif
    }
  }


#pragma omp parallel for
  ILOOPALL
  {
    JSLOOP(N2, N2-1+NG)
    {
#if N2 < NG
      int jactive = N2 - 1 + NG;
      PLOOP S->P[ip][j][i] = S->P[ip][jactive][i];
      pflag[j][i] = pflag[jactive][i];
#elif X2R_BOUND == OUTFLOW
      int jz = N2 - 1 + NG;
      PLOOP S->P[ip][j][i] = S->P[ip][jz][i];
      pflag[j][i] = pflag[jz][i];
#elif X2R_BOUND == POLAR
      // As j grows beyond N2+NG, reflect the zone that far previous
      int jrefl = (N2 + NG) + (N2 + NG - j) - 1;
      PLOOP S->P[ip][j][i] = S->P[ip][jrefl][i];
      pflag[j][i] = pflag[jrefl][i];
      S->P[U2][j][i] *= -1.;
      S->P[B2][j][i] *= -1.;
#elif X2R_BOUND == PERIODIC
			int jz = j - N2;
			PLOOP S->P[ip][j][i] = S->P[ip][jz][i];
			pflag[j][i] = pflag[jz][i];
#endif
    }
  }

  timer_stop(TIMER_BOUND);
}

#if METRIC == MKS
/**
 * @brief Enforce a no-inflow (or no-outflow) condition on a single ghost zone.
 *
 * @details If the radial 4-velocity u^r points the wrong way for the given
 * boundary type, resets U1 so that u^r = 0 (the fluid is at rest in the
 * normal-observer frame radially):
 *  - type == 0 (inner boundary): Zero u^r when fluid would flow inward (u^r > 0).
 *  - type == 1 (outer boundary): Zero u^r when fluid would flow outward (u^r < 0).
 *
 * The algorithm strips gamma from the 3-velocities, sets U1 = beta^r/alpha
 * (the value that makes u^r = 0), then re-applies the new Lorentz factor to
 * keep the primitive velocities normalised.
 *
 * @param G     Grid geometry.
 * @param S     Fluid state (P[U1-U3] is read and possibly modified).
 * @param i     X1 zone index (ghost zone).
 * @param j     X2 zone index.
 * @param type  0 = inner (X1 left) boundary; 1 = outer (X1 right) boundary.
 */
void inflow_check(struct GridGeom *G, struct FluidState *S, int i, int j, int type)
{
  double alpha, beta1, vsq;
  ucon_calc(G, S, i, j, CENT);

  if (((S->ucon[1][j][i] > 0.) && (type == 0)) ||
      ((S->ucon[1][j][i] < 0.) && (type == 1)))
  {
    // Find gamma and remove it from S->Pitives
    double gamma = mhd_gamma_calc(G, S, i, j, CENT);
    S->P[U1][j][i] /= gamma;
    S->P[U2][j][i] /= gamma;
    S->P[U3][j][i] /= gamma;
    alpha = G->lapse[CENT][j][i];
    beta1 = G->gcon[CENT][0][1][j][i]*alpha*alpha;

    // Reset radial velocity so radial 4-velocity is zero
    S->P[U1][j][i] = beta1/alpha;

    // Now find new gamma and put it back in
    vsq = 0.;
    for (int mu = 1; mu < NDIM; mu++) 
    {
      for (int nu = 1; nu < NDIM; nu++)
      {
        vsq += G->gcov[CENT][mu][nu][j][i]*S->P[U1+mu-1][j][i]*S->P[U1+nu-1][j][i];
      }
    }
    if (fabs(vsq) < 1.e-13)
      vsq = 1.e-13;
    if (vsq >= 1.)
      vsq = 1. - 1./(GAMMAMAX*GAMMAMAX);
    gamma = 1./sqrt(1. - vsq);
    S->P[U1][j][i] *= gamma;
    S->P[U2][j][i] *= gamma;
    S->P[U3][j][i] *= gamma;
  }
}

/**
 * @brief Apply flux corrections at radial and polar boundaries before the conservative update.
 *
 * @details Called by advance_fluid() before the divergence-form conservative update.
 * Performs three independent fixups:
 *
 * 1. **Inflow clip at inner boundary** (if X1L_INFLOW == 0):
 *    F^{X1}_{RHO}[j][NG] = min(F, 0) — ensures mass can only flow outward through the inner face.
 *
 * 2. **Outflow clip at outer boundary** (if X1R_INFLOW == 0):
 *    F^{X1}_{RHO}[j][N1+NG] = max(F, 0) — ensures mass can only flow inward through the outer face.
 *
 * 3. **Polar axis:**
 *    - All X2 fluxes at both poles are zeroed (no transport perpendicular to the pole in spherical KS).
 *    - The X1 flux of B2 is anti-reflected across the poles (one ghost zone) to maintain the
 *      field topology expected by constrained transport.
 *
 * @param F  Fluid flux struct; X1[] and X2[] arrays are modified in-place.
 */
void fix_flux(struct FluidFlux *F)
{
  if (X1L_INFLOW == 0)
  {

#pragma omp parallel for
    JLOOPALL
      F->X1[RHO][j][0+NG] = MY_MIN(F->X1[RHO][j][0+NG], 0.);
  }

  if (X1R_INFLOW == 0)
  {

#pragma omp parallel for
    JLOOPALL
      F->X1[RHO][j][N1+NG] = MY_MAX(F->X1[RHO][j][N1+NG], 0.);
  }


#pragma omp parallel for
  ILOOPALL
  {
    F->X1[B2][-1+NG][i] = -F->X1[B2][0+NG][i];
    PLOOP F->X2[ip][0+NG][i] = 0.;
  }

#pragma omp parallel for
  ILOOPALL
  {
    F->X1[B2][N2+NG][i] = -F->X1[B2][N2-1+NG][i];
    PLOOP F->X2[ip][N2+NG][i] = 0.;
  }
}
#endif // METRIC
