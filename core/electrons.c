/**
 * @file electrons.c
 * @brief Electron thermodynamics: entropy initialization, viscous heating, and fixups.
 *
 * @details Implements the sub-grid electron heating models described in
 * Sadowski et al. (2017, MNRAS 456, 1837) and related papers.  All code in
 * this file is conditionally compiled with `#if ELECTRONS`.
 *
 * The electron state is tracked via specific entropy variables K_el per model:
 * @f[ K_{\rm el} = \frac{(\gamma_e - 1)\,u_e}{\rho^{\gamma_e}} @f]
 * and a total-gas entropy K_tot:
 * @f[ K_{\rm tot} = \frac{(\gamma - 1)\,u}{\rho^\gamma} @f]
 *
 * **Functions:**
 * - **init_electrons()**: Sets K_tot and all K_el[] from the initial primitive
 *   state using a constant electron fraction fel0 of the total internal energy.
 * - **heat_electrons()**: After each time step, for each zone computes the
 *   heating fraction f_el (from one of four models) and updates K_el via the
 *   entropy jump as fluid elements are compressed/shocked.
 * - **heat_electrons_1zone()**: Per-zone kernel for heat_electrons().
 * - **get_fels()**: Computes the electron heating fraction f_el for a given
 *   zone and model index (KAWAZURA, WERNER, ROWAN, or SHARMA).
 * - **fixup_electrons()**: Clamps K_el to lie within [kelmin, kelmax] to
 *   enforce T_i/T_e bounds (tptemin ≤ T_i/T_e ≤ tptemax) and replace NANs.
 * - **fixup_electrons_1zone()**: Per-zone kernel for fixup_electrons().
 *
 * @note When ALLMODELS=1, all four heating models are evolved simultaneously
 * as separate primitive variables KEL0, KEL1, KEL2, KEL3.  When ALLMODELS=0,
 * only KAWAZURA is used.
 */

/*---------------------------------------------------------------------------------

  ELECTRONS.C

  -Initialize electron and gas entropies
  -Assign electron and total entropies based on https://academic.oup.com/mnras/article/454/2/1848/2892599

---------------------------------------------------------------------------------*/

#include "decs.h"

#if ELECTRONS

// Defined as in decs.h, CONSTANT not included in ALLMODELS version
// KAWAZURA is run by default if ALLMODELS=0 
#define KAWAZURA  9
#define WERNER    10
#define ROWAN     11
#define SHARMA    12
#define CONSTANT 5 //tbh, this is never considered 

void fixup_electrons_1zone(struct FluidState *S, int i, int j);
void heat_electrons_1zone(struct GridGeom *G, struct FluidState *Sh, struct FluidState *S, int i, int j);
double get_fels(struct GridGeom *G, struct FluidState *S, int i, int j, int model);

/**
 * @brief Initialize electron and total entropy variables from the initial primitives.
 *
 * @details For every zone, sets:
 * - K_tot = (γ-1)*u / ρ^γ  (total gas entropy).
 * - All K_el (KEL0..NVAR-1) = (γ_e-1)*fel0*u / ρ^γ_e
 *   (all models start at a constant fraction fel0 of the total internal energy).
 * Calls set_bounds() at the end to propagate initial values into ghost zones.
 *
 * @param G  Grid geometry.
 * @param S  Fluid state (P[RHO], P[UU] read; P[KTOT] and P[KEL*] written).
 */
// Initialize electron entropies like (k = P / rho^gam)
void init_electrons(struct GridGeom *G, struct FluidState *S)
{
  ZLOOPALL {
    // Set electron internal energy to constant fraction of internal energy
    double uel = fel0*S->P[UU][j][i];

    // Initialize entropies
    S->P[KTOT][j][i] = (gam-1.)*S->P[UU][j][i]*pow(S->P[RHO][j][i],-gam);

    // Initialize model entropy(ies)
    for (int idx = KEL0; idx < NVAR ; idx++) {
      S->P[idx][j][i] = (game-1.)*uel*pow(S->P[RHO][j][i],-game);
    }
  }

  // Necessary?  Usually called right afterward
  set_bounds(G, S);
}

/**
 * @brief Apply electron heating to all active zones.
 *
 * @details OpenMP-parallelised loop over ZLOOP that calls heat_electrons_1zone()
 * for each active zone.
 *
 * @param G   Grid geometry.
 * @param Ss  Fluid state at start of step (used to compute heating fraction f_el).
 * @param Sf  Fluid state at end of step (K_el and K_tot are updated here).
 */
// Heat electrons over the physical domain
void heat_electrons(struct GridGeom *G, struct FluidState *Ss, struct FluidState *Sf)
{
  timer_start(TIMER_ELECTRON_HEAT);

#pragma omp parallel for collapse(2)
  ZLOOP {
    heat_electrons_1zone(G, Ss, Sf, i, j);
  }

  timer_stop(TIMER_ELECTRON_HEAT);
}

// Heat electrons at given zone
inline void heat_electrons_1zone(struct GridGeom *G, struct FluidState *Ss, struct FluidState *Sf, int i, int j)
{
  // Actual entropy at final time
  double kHarm = (gam-1.)*Sf->P[UU][j][i]/pow(Sf->P[RHO][j][i],gam);

  // Evolve model entropy(ies)
  for (int idx = KEL0; idx < NVAR ; idx++) {
    double fel = get_fels(G, Ss, i, j, idx);
    Sf->P[idx][j][i] += (game-1.)/(gam-1.)*pow(Ss->P[RHO][j][i],gam-game)*fel*(kHarm - Sf->P[KTOT][j][i]);
  }

  // Reset total entropy
  Sf->P[KTOT][j][i] = kHarm;
}

/**
 * @brief Compute the electron heating fraction f_el at a single zone.
 *
 * @details Implements four published sub-grid models, selected by the `model`
 * parameter (which equals the primitive index KEL0, KEL1, etc.):
 * - **KAWAZURA** (model == 9): Eq. (2) of Kawazura et al. (2019, PNAS).
 *   Uses T_p/T_e ratio and plasma β to determine Q_i/Q_e.
 * - **WERNER** (model == 10): Eq. (3) of Werner et al. (2018, MNRAS).
 *   Uses magnetisation σ = b²/ρ.
 * - **ROWAN** (model == 11): Eq. (34) of Rowan et al. (2017, ApJ).
 *   Uses β and σ relative to a maximum β.
 * - **SHARMA** (model == 12): Section 4 of Sharma et al. (2007, ApJ).
 *   Uses T_e/T_p ratio (inverse of KAWAZURA convention).
 * When SUPPRESS_HIGHB_HEAT is defined, f_el is set to 0 for σ > 1.
 *
 * @param G      Grid geometry.
 * @param S      Fluid state (P[RHO], P[UU], P[model] read).
 * @param i      X1 zone index.
 * @param j      X2 zone index.
 * @param model  Primitive index corresponding to the heating model.
 * @return       Electron heating fraction 0 ≤ f_el ≤ 0.5.
 */
// New function for ALLMODELS runs.
inline double get_fels(struct GridGeom *G, struct FluidState *S, int i, int j, int model)
{
  get_state(G, S, i, j, CENT);
  double bsq = bsq_calc(S, i, j);
  double fel = 0.0;
if (model == KAWAZURA) {
	// Equation (2) in http://www.pnas.org/lookup/doi/10.1073/pnas.1812491116
  double Tpr = (gamp-1.)*S->P[UU][j][i]/S->P[RHO][j][i];
  double uel = 1./(game-1.)*S->P[model][j][i]*pow(S->P[RHO][j][i],game);
  double Tel = (game-1.)*uel/S->P[RHO][j][i];
  if(Tel <= 0.) Tel = SMALL;
  if(Tpr <= 0.) Tpr = SMALL;

  double Trat = fabs(Tpr/Tel);
  double pres = S->P[RHO][j][i]*Tpr; // Proton pressure
  double beta = pres/bsq*2;
  if(beta > 1.e20) beta = 1.e20;
  
  double QiQe = 35./(1. + pow(beta/15.,-1.4)*exp(-0.1/Trat));
  fel = 1./(1. + QiQe);
} else if (model == WERNER) {
	// Equation (3) in http://academic.oup.com/mnras/article/473/4/4840/4265350
  double sigma = bsq/S->P[RHO][j][i];
  fel = 0.25*(1+pow(((sigma/5.)/(2+(sigma/5.))), .5));
} else if (model == ROWAN) {
	// Equation (34) in https://iopscience.iop.org/article/10.3847/1538-4357/aa9380
  double pres = (gamp-1.)*S->P[UU][j][i]; // Proton pressure
  double pg = (gam-1)*S->P[UU][j][i];
  double beta = pres/bsq*2;
  double sigma = bsq/(S->P[RHO][j][i]+S->P[UU][j][i]+pg);
  double betamax = 0.25/sigma;
  fel = 0.5*exp(-pow(1-beta/betamax, 3.3)/(1+1.2*pow(sigma, 0.7)));
} else if (model == SHARMA) {
	// Equation for \delta on  pg. 719 (Section 4) in https://iopscience.iop.org/article/10.1086/520800
  double Tpr = (gamp-1.)*S->P[UU][j][i]/S->P[RHO][j][i];
  double uel = 1./(game-1.)*S->P[model][j][i]*pow(S->P[RHO][j][i],game);
  double Tel = (game-1.)*uel/S->P[RHO][j][i];
  if(Tel <= 0.) Tel = SMALL;
  if(Tpr <= 0.) Tpr = SMALL;

  double Trat_inv = fabs(Tel/Tpr); //Inverse of the temperature ratio in KAWAZURA
  double QeQi = 0.33 * pow(Trat_inv, 0.5);
	fel = 1./(1.+1./QeQi);
}

// Avoid electron heating at high sigma
#if SUPPRESS_HIGHB_HEAT
  if(bsq/S->P[RHO][j][i] > 1.) fel = 0;
#endif

  return fel;
}

/**
 * @brief Clamp electron entropy variables to physical bounds in all active zones.
 *
 * @details Parallelised over ZLOOP; calls fixup_electrons_1zone() per zone to
 * enforce tptemin ≤ T_p/T_e ≤ tptemax and replace any NaN values.
 *
 * @param S  Fluid state (P[KEL*] modified in-place).
 */
// Fix electron if either (i) the heating prescription provied NaN for some reason, or
// (ii) the ratio of ion to electron temperature was too large, or
// (iii) the ratio of ion to electron temperature was too small
void fixup_electrons(struct FluidState *S)
{
  timer_start(TIMER_ELECTRON_FIXUP);

#pragma omp parallel for collapse(2)
  ZLOOP {
    fixup_electrons_1zone(S, i, j);
  }

  timer_stop(TIMER_ELECTRON_FIXUP);
}

// Fix for a given zone
inline void fixup_electrons_1zone(struct FluidState *S, int i, int j)
{
  double kelmax = S->P[KTOT][j][i]*pow(S->P[RHO][j][i],gam-game)/(tptemin*(gam-1.)/(gamp-1.) + (gam-1.)/(game-1.));
  double kelmin = S->P[KTOT][j][i]*pow(S->P[RHO][j][i],gam-game)/(tptemax*(gam-1.)/(gamp-1.) + (gam-1.)/(game-1.));

  // Replace NANs with cold electrons
  for (int idx = KEL0; idx < NVAR ; idx++) {
    if (isnan(S->P[idx][j][i])) S->P[idx][j][i] = kelmin;
	// Enforce maximum Tp/Te
    S->P[idx][j][i] = MY_MAX(S->P[idx][j][i], kelmin);
	// Enforce minimum Tp/Te
    S->P[idx][j][i] = MY_MIN(S->P[idx][j][i], kelmax);
  }
}
#endif // ELECTRONS
