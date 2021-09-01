/*---------------------------------------------------------------------------------

  STEP.C

  -Advances simulation by one timestep

---------------------------------------------------------------------------------*/

#include "decs.h"

// Declarations
double advance_fluid(struct GridGeom *G, struct FluidState *Si, struct FluidState *Ss, struct FluidState *Sf, double Dt);

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
  // Work around ICC 18.0.2 bug in assigning to pointers to structs
#if INTEL_WORKAROUND
  memcpy(&(Ssave->P),&(S->P),sizeof(GridPrim));
#else
#pragma omp parallel for simd collapse(2)
  PLOOP ZLOOPALL Ssave->P[ip][j][i] = S->P[ip][j][i];
#endif
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

  // Work around ICC 18.0.2 bug in assigning to pointers to structs
#if INTEL_WORKAROUND
  memcpy(&(Sf->P),&(Si->P),sizeof(GridPrim));
#else
#pragma omp parallel for simd collapse(2)
  PLOOP ZLOOPALL Sf->P[ip][j][i] = Si->P[ip][j][i];
#endif

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
