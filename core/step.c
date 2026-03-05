/**
 * @file step.c
 * @brief Second-order predictor-corrector time integration.
 *
 * @details Implements a two-stage Runge-Kutta (RK2) algorithm to advance the
 * GRMHD equations from time @f$t^n@f$ to @f$t^{n+1} = t^n + \Delta t@f$.
 *
 * **Algorithm overview:**
 * 1. **Save** primitives at @f$t^n@f$ into Ssave (needed for 4-current calculation).
 * 2. **Predictor** – advance_fluid(S, S, Stmp, dt/2):
 *    - Reconstruct primitives at faces from @f$P^n@f$.
 *    - Solve the LLF Riemann problem to obtain fluxes @f$F^{n+1/2}@f$.
 *    - Apply flux-CT to maintain @f$\nabla \cdot B = 0@f$.
 *    - Update conserved variables: @f$U^* = U^n + \frac{\Delta t}{2}(\nabla \cdot F + S)@f$.
 *    - Invert @f$U^*@f$ to @f$P^*@f$ (U_to_P).
 *    - Heat electrons (if ELECTRONS is enabled).
 *    - Apply floors/ceilings (fixup), then boundary conditions (set_bounds).
 * 3. **Corrector** – advance_fluid(S, Stmp, S, dt):
 *    - Same sequence using @f$P^*@f$ as the source state and @f$P^n@f$ as the base.
 *    - Result is @f$P^{n+1}@f$ stored back in S.
 * 4. **4-current** – computed from @f$(P^n + P^{n+1})/2@f$ if a dump is imminent.
 * 5. **Timestep update** – new dt = min(ndt, SAFE * dt_old).
 *
 * @note The function advance_fluid() is declared @c inline and is local to this file.
 */

/*---------------------------------------------------------------------------------

  STEP.C

  -Advances simulation by one timestep

---------------------------------------------------------------------------------*/

#include "decs.h"

// Declarations
double advance_fluid(struct GridGeom *G, struct FluidState *Si, struct FluidState *Ss, struct FluidState *Sf, double Dt);

/**
 * @brief Advance the simulation state by one full timestep.
 *
 * @details Executes the predictor-corrector (RK2) integration scheme.
 * Allocates Stmp and Ssave on the first call (static; persists across calls).
 * After the corrector, increments the global time @c t by @c dt and
 * updates @c dt for the next step using the CFL estimate from advance_fluid().
 *
 * @param G  Pointer to the precomputed grid geometry.
 * @param S  Pointer to the fluid state (updated in-place to @f$P^{n+1}@f$).
 */
void step(struct GridGeom *G, struct FluidState *S)
{
  static struct FluidState *Stmp;
  static struct FluidState *Ssave;
  static int first_call = 1;
  if (first_call)
  {
    Stmp = calloc(1,sizeof(struct FluidState));
    Ssave = calloc(1,sizeof(struct FluidState));
    first_call = 0;
  }

  // Need both P_n and P_n+1 to calculate current

#pragma omp parallel for simd collapse(2)
  PLOOP ZLOOPALL Ssave->P[ip][j][i] = S->P[ip][j][i];

  LOGN("Step %d",nstep);
  FLAG("Start step");

  // Predictor setup
  advance_fluid(G, S, S, Stmp, 0.5*dt);
  FLAG("Advance Fluid Tmp");
#if ELECTRONS
  heat_electrons(G, S, Stmp);
  FLAG("Heat Electrons Tmp");
#endif

  // Fixup routines: smooth over outlier zones
  fixup(G, Stmp);
  FLAG("Fixup Tmp");
#if ELECTRONS
  fixup_electrons(Stmp);
  FLAG("Fixup e- Tmp");
#endif
  set_bounds(G, Stmp);
  FLAG("First bounds Tmp");
  fixup_utoprim(G, Stmp);
  FLAG("Fixup U_to_P Tmp");
  set_bounds(G, Stmp);
  FLAG("Second bounds Tmp");

  // Corrector step
  double ndt = advance_fluid(G, S, Stmp, S, dt);
  FLAG("Advance Fluid Full");

#if ELECTRONS
  heat_electrons(G, Stmp, S);
  FLAG("Heat Electrons Full");
#endif

  fixup(G, S);
  FLAG("Fixup Full");
#if ELECTRONS
  fixup_electrons(S);
  FLAG("Fixup e- Full");
#endif
  set_bounds(G, S);
  FLAG("First bounds Full");
  fixup_utoprim(G, S);
  FLAG("Fixup U_to_P Full");
  set_bounds(G, S);
  FLAG("Second bounds Full");

  // Increment time
  t += dt;

  // If we're dumping this step, update the current
  if (t >= tdump) {
    current_calc(G, S, Ssave, dt);
  }

  // Set next timestep
  if (ndt > SAFE * dt)
    dt = SAFE * dt;
  else
    dt = ndt;
}

/**
 * @brief Perform a single sub-step of the time integration.
 *
 * @details Updates the fluid state from Si (base) using fluxes evaluated from Ss
 * (source), and stores the result in Sf.  The update equation is:
 * @f[
 *   U_f = U_i + \Delta t \left(\frac{F^{X1}_{i,j} - F^{X1}_{i+1,j}}{\Delta X^1}
 *       + \frac{F^{X2}_{i,j} - F^{X2}_{i,j+1}}{\Delta X^2} + S_{ij}\right)
 * @f]
 * After updating conserved variables, inverts back to primitives via U_to_P().
 *
 * @param G   Grid geometry.
 * @param Si  Initial fluid state (provides base conserved variables U_i).
 * @param Ss  Source fluid state (primitives used for reconstruction and flux evaluation).
 * @param Sf  Output fluid state (receives updated primitives P_f and conserved U_f).
 * @param Dt  Sub-step size in code time units.
 * @return    CFL-limited next timestep estimate.
 */
inline double advance_fluid(struct GridGeom *G, struct FluidState *Si, struct FluidState *Ss, struct FluidState *Sf, double Dt)
{
  static GridPrim *dU;
  static struct FluidFlux *F;

  static int firstc = 1;
  if (firstc)
  {
    dU = calloc(1,sizeof(GridPrim));
    F = calloc(1,sizeof(struct FluidFlux));
    firstc = 0;
  }


#pragma omp parallel for simd collapse(2)
  PLOOP ZLOOPALL Sf->P[ip][j][i] = Si->P[ip][j][i];

  double ndt = get_flux(G, Ss, F);

#if METRIC == MKS
  fix_flux(F);
#endif

  //Constrained transport for B
  flux_ct(F);

  // Flux diagnostic globals
  diag_flux(F);

  // Update Si to Sf
  timer_start(TIMER_UPDATE_U);
  get_state_vec(G, Ss, CENT, 0, N2 - 1, 0, N1 - 1);
  get_fluid_source(G, Ss, dU);

  get_state_vec(G, Si, CENT, 0, N2 - 1, 0, N1 - 1);
  prim_to_flux_vec(G, Si, 0, CENT, 0, N2 - 1, 0, N1 - 1, Si->U);

#pragma omp parallel for collapse(2)
  PLOOP ZLOOP 
  {
    Sf->U[ip][j][i] = Si->U[ip][j][i] +
      Dt*((F->X1[ip][j][i] - F->X1[ip][j][i+1])/dx[1] +
          (F->X2[ip][j][i] - F->X2[ip][j+1][i])/dx[2] +
          (*dU)[ip][j][i]);
  }
  timer_stop(TIMER_UPDATE_U);

  timer_start(TIMER_U_TO_P);
#pragma omp parallel for collapse(2)
  ZLOOP
    pflag[j][i] = U_to_P(G, Sf, i, j, CENT);
  
  timer_stop(TIMER_U_TO_P);

#pragma omp parallel for simd
  ZLOOPALL
    fail_save[j][i] = pflag[j][i];

  return ndt;
}
