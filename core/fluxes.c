/*---------------------------------------------------------------------------------

  FLUXES.C

  -Compute fluid fluxes
  -Apply constrained transport

---------------------------------------------------------------------------------*/

#include "decs.h"

void lr_to_flux(struct GridGeom *G, struct FluidState *Sl, struct FluidState *Sr, int dir, int loc, GridPrim *flux, GridVector *ctop);
double ndt_min(GridVector *ctop);

double ndt_min(GridVector *ctop)
{
  timer_start(TIMER_CMAX);
  double ndt_min = 1e20;

#if DEBUG
  int min_x1, min_x2;
#endif

#pragma omp parallel for collapse(2) reduction(min:ndt_min)
  ZLOOP
  {
    double ndt_zone = 0;
    for (int mu = 1; mu < NDIM-1; mu++)
      ndt_zone += 1/(cour*dx[mu]/(*ctop)[mu][j][i]);
    ndt_zone = 1/ndt_zone;

    if(ndt_zone < ndt_min)
    {
      ndt_min = ndt_zone;
#if DEBUG
      min_x1 = i;
      min_x2 = j;
#endif
    }
  }

#if DEBUG
  fprintf(stderr, "Timestep set by %d %d\n",min_x1, min_x2);
#endif

  timer_stop(TIMER_CMAX);
  return ndt_min;
}

double get_flux(struct GridGeom *G, struct FluidState *S, struct FluidFlux *F)
{
  static struct FluidState *Sl, *Sr;
  static GridVector *ctop;
  double cmax[NDIM], ndts[NDIM];
  memset(cmax, 0, NDIM*sizeof(double));
  memset(ndts, 0, NDIM*sizeof(double));

  static int firstc = 1;

  if (firstc)
  {
    Sl  = calloc(1,sizeof(struct FluidState));
    Sr  = calloc(1,sizeof(struct FluidState));
    ctop = calloc(1,sizeof(GridVector));
    firstc = 0;
  }

  // reconstruct X1
  reconstruct(S, Sl->P, Sr->P, 1);

  // lr_to_flux X1
  lr_to_flux(G, Sl, Sr, 1, FACE1, &(F->X1), ctop);

  // reconstruct X2
  reconstruct(S, Sl->P, Sr->P, 2);

  // lr_to_flux X2
  lr_to_flux(G, Sl, Sr, 2, FACE2, &(F->X2), ctop);

  return ndt_min(ctop);
}

// Note that the sense of L/R flips from zone to interface during function call
void lr_to_flux(struct GridGeom *G, struct FluidState *Sr, struct FluidState *Sl, int dir, int loc, GridPrim *flux, GridVector *ctop)
{
  timer_start(TIMER_LR_TO_F);

  static GridPrim *fluxL, *fluxR;
  static GridDouble *cmaxL, *cmaxR, *cminL, *cminR, *cmax, *cmin;

  static int firstc = 1;

  if (firstc)
  {
    fluxL = calloc(1,sizeof(GridPrim));
    fluxR = calloc(1,sizeof(GridPrim));
    cmaxL = calloc(1,sizeof(GridDouble));
    cmaxR = calloc(1,sizeof(GridDouble));
    cminL = calloc(1,sizeof(GridDouble));
    cminR = calloc(1,sizeof(GridDouble));
    cmax = calloc(1,sizeof(GridDouble));
    cmin = calloc(1,sizeof(GridDouble));

    firstc = 0;
  }

  // Properly offset left face
  // These are un-macro'd to bundle OpenMP thread tasks rather than memory accesses
  PLOOP 
  {
    if (dir == 1) 
    {
#pragma omp parallel for
      ZSLOOP_REVERSE(-1, N2, -1, N1)
        Sl->P[ip][j][i] = Sl->P[ip][j][i - 1];
    } 
    else if (dir == 2) 
    {
#pragma omp parallel for
      for (int i = (N1) + NG; i >= (-1) + NG; i--)
        for (int j = (N2) + NG; j >= (-1) + NG; j--)
          Sl->P[ip][j][i] = Sl->P[ip][j - 1][i];
    } 
  }

  timer_start(TIMER_LR_STATE);

  get_state_vec(G, Sl, loc, -1, N2, -1, N1);
  get_state_vec(G, Sr, loc, -1, N2, -1, N1);

  timer_stop(TIMER_LR_STATE);

  timer_start(TIMER_LR_PTOF);

  prim_to_flux_vec(G, Sl, 0,   loc, -1, N2, -1, N1, Sl->U);
  prim_to_flux_vec(G, Sl, dir, loc, -1, N2, -1, N1, *fluxL);

  prim_to_flux_vec(G, Sr, 0,   loc, -1, N2, -1, N1, Sr->U);
  prim_to_flux_vec(G, Sr, dir, loc, -1, N2, -1, N1, *fluxR);

  timer_stop(TIMER_LR_PTOF);

  timer_start(TIMER_LR_VCHAR);
#pragma omp parallel
  {
#pragma omp for nowait
    ZSLOOP(-1, N2, -1, N1)
        mhd_vchar(G, Sl, i, j, loc, dir, *cmaxL, *cminL);
#pragma omp for
    ZSLOOP(-1, N2, -1, N1)
        mhd_vchar(G, Sr, i, j, loc, dir, *cmaxR, *cminR);
  }
  timer_stop(TIMER_LR_VCHAR);

  timer_start(TIMER_LR_CMAX);
#pragma omp parallel for
  ZSLOOP(-1, N2, -1, N1) 
  {
    (*cmax)[j][i] = fabs(MY_MAX(MY_MAX(0., (*cmaxL)[j][i]), (*cmaxR)[j][i]));
    (*cmin)[j][i] = fabs(MY_MAX(MY_MAX(0., -(*cminL)[j][i]), -(*cminR)[j][i]));
    (*ctop)[dir][j][i] = MY_MAX((*cmax)[j][i], (*cmin)[j][i]);
    
    if (isnan(1./(*ctop)[dir][j][i])) 
    {
      printf("ctop is 0 or NaN at zone: %i %i (%i) ", i, j, dir);
#if METRIC == MKS
      double X[NDIM];
      double r, th;
      coord(i, j, CENT, X);
      bl_coord(X, &r, &th);
      printf("(r,th = %f %f)\n", r, th);
#endif
      printf("\n");
      exit(-1);
    }
  }
  timer_stop(TIMER_LR_CMAX);

  timer_start(TIMER_LR_FLUX);
#pragma omp parallel for simd collapse(2)
  PLOOP
    ZSLOOP(-1, N2, -1, N1)
      (*flux)[ip][j][i] = 0.5*((*fluxL)[ip][j][i] + (*fluxR)[ip][j][i] - (*ctop)[dir][j][i]*(Sr->U[ip][j][i] - Sl->U[ip][j][i]));
  timer_stop(TIMER_LR_FLUX);
  timer_stop(TIMER_LR_TO_F);
}

void flux_ct(struct FluidFlux *F)
{
  timer_start(TIMER_FLUX_CT);

  static GridDouble emf;
  static int firstc = 1;
  if (firstc)
  {
    firstc = 0;
  }

#pragma omp parallel
  {
#pragma omp for simd
    ZSLOOP(0, N2, 0, N1)
      emf[j][i] =  0.25*(F->X1[B2][j][i] + F->X1[B2][j-1][i] - F->X2[B1][j][i] - F->X2[B1][j][i-1]);

    // Rewrite EMFs as fluxes, after Toth
#pragma omp for simd nowait
    ZSLOOP(0, N2 - 1, 0, N1)
    {
      F->X1[B1][j][i] =  0.;
      F->X1[B2][j][i] =  0.5*(emf[j][i] + emf[j+1][i]);
    }
#pragma omp for simd nowait
    ZSLOOP(0, N2, 0, N1 - 1)
    {
      F->X2[B1][j][i] = -0.5*(emf[j][i] + emf[j][i+1]);
      F->X2[B2][j][i] =  0.;
    }
  } // omp parallel
  timer_stop(TIMER_FLUX_CT);
}
