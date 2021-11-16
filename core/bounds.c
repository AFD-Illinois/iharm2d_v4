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
void set_bounds(struct GridGeom *G, struct FluidState *S)
{
  timer_start(TIMER_BOUND);

#if !INTEL_WORKAROUND
#pragma omp parallel for
#endif
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
#if !INTEL_WORKAROUND
#pragma omp parallel for
#endif
    JLOOP
      ISLOOP(-NG, -1)
        inflow_check(G, S, i, j, 0);
  }
#endif

#if !INTEL_WORKAROUND
#pragma omp parallel for
#endif
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
#if !INTEL_WORKAROUND
#pragma omp parallel for
#endif   
    JLOOP
      ISLOOP(N1, N1 - 1 + NG)
        inflow_check(G, S, i, j, 1);
  }
#endif

#if !INTEL_WORKAROUND
#pragma omp parallel for
#endif
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

#if !INTEL_WORKAROUND
#pragma omp parallel for
#endif
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

void fix_flux(struct FluidFlux *F)
{
  if (X1L_INFLOW == 0)
  {
#if !INTEL_WORKAROUND
#pragma omp parallel for
#endif
    JLOOPALL
      F->X1[RHO][j][0+NG] = MY_MIN(F->X1[RHO][j][0+NG], 0.);
  }

  if (X1R_INFLOW == 0)
  {
#if !INTEL_WORKAROUND
#pragma omp parallel for
#endif
    JLOOPALL
      F->X1[RHO][j][N1+NG] = MY_MAX(F->X1[RHO][j][N1+NG], 0.);
  }

#if !INTEL_WORKAROUND
#pragma omp parallel for
#endif
  ILOOPALL
  {
    F->X1[B2][-1+NG][i] = -F->X1[B2][0+NG][i];
    PLOOP F->X2[ip][0+NG][i] = 0.;
  }

#if !INTEL_WORKAROUND
#pragma omp parallel for
#endif
  ILOOPALL
  {
    F->X1[B2][N2+NG][i] = -F->X1[B2][N2-1+NG][i];
    PLOOP F->X2[ip][N2+NG][i] = 0.;
  }
}
#endif // METRIC
