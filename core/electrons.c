/*---------------------------------------------------------------------------------

  ELECTRONS.C

  -Initialize electron and gas entropies
  -Assign electron and total entropies based on https://academic.oup.com/mnras/article/454/2/1848/2892599

---------------------------------------------------------------------------------*/

#include "decs.h"

#if ELECTRONS

// TODO put these in options with a default in decs.h
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

// TODO merge these
void heat_electrons(struct GridGeom *G, struct FluidState *Ss, struct FluidState *Sf)
{
  timer_start(TIMER_ELECTRON_HEAT);

#pragma omp parallel for collapse(2)
  ZLOOP {
    heat_electrons_1zone(G, Ss, Sf, i, j);
  }

  timer_stop(TIMER_ELECTRON_HEAT);
}

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

#if SUPPRESS_HIGHB_HEAT
  if(bsq/S->P[RHO][j][i] > 1.) fel = 0;
#endif

  return fel;
}

void fixup_electrons(struct FluidState *S)
{
  timer_start(TIMER_ELECTRON_FIXUP);

#pragma omp parallel for collapse(2)
  ZLOOP {
    fixup_electrons_1zone(S, i, j);
  }

  timer_stop(TIMER_ELECTRON_FIXUP);
}

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
