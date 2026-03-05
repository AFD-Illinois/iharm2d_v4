/**
 * @file fluxes.c
 * @brief Numerical flux computation, constrained transport, and CFL timestep.
 *
 * @details This module is responsible for:
 *
 * - **Reconstruction + Riemann solve** (get_flux()): For each coordinate direction,
 *   reconstructs left/right primitive states at cell faces via reconstruct(), then
 *   calls lr_to_flux() to compute the Local Lax-Friedrichs (LLF/Rusanov) flux:
 *   @f[ F_{i+1/2} = \frac{1}{2}\left(F_L + F_R - c_\text{max}(U_R - U_L)\right) @f]
 *   where @f$c_\text{max}@f$ is the maximum fast magnetosonic wave speed at the face.
 *
 * - **Constrained transport** (flux_ct()): Implements the Toth (2000) flux-CT scheme
 *   to maintain the divergence-free condition @f$\nabla \cdot B = 0@f$ to machine precision.
 *   Computes the electric field (EMF) at cell corners from the face-centered magnetic flux,
 *   then rewrites the @f$B@f$-component fluxes using those EMFs.
 *
 * - **CFL timestep** (ndt_min()): Computes the next timestep as the global minimum of
 *   @f$ \Delta t = \text{cour} / \sum_i c_{\text{max},i} / \Delta X^i @f$ over all zones.
 *
 * @see reconstruction.c for the spatial interpolation methods.
 * @see phys.c for mhd_vchar() which supplies the wave speeds.
 */

/*---------------------------------------------------------------------------------

  FLUXES.C

  -Compute fluid fluxes, maximum wavespeed and timestep
  -Apply constrained transport

---------------------------------------------------------------------------------*/

#include "decs.h"

void lr_to_flux(struct GridGeom *G, struct FluidState *Sl, struct FluidState *Sr, int dir, int loc, GridPrim *flux, GridVector *ctop);
double ndt_min(GridVector *ctop);

/**
 * @brief Compute the CFL-limited next timestep from the maximum wave speed grid.
 *
 * @details Performs a parallel reduction over all active zones to find the zone that
 * constrains the timestep the most.  For each zone the local constraint is:
 * @f[ \Delta t_{\rm zone} = \frac{1}{\sum_{i=1}^{2} c_{{\rm top},i} / (\text{cour} \cdot \Delta X^i)} @f]
 *
 * @param ctop  Grid of maximum wave speeds in each direction: (*ctop)[mu][j][i].
 * @return      CFL-limited timestep.
 */
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

/**
 * @brief Reconstruct primitives at faces and compute LLF numerical fluxes.
 *
 * @details For each coordinate direction (X1 and X2):
 * 1. Calls reconstruct() to build left and right primitive states at all faces.
 * 2. Calls lr_to_flux() to evaluate the LLF Riemann flux and fill ctop.
 * 3. Returns the minimum CFL timestep from ndt_min().
 *
 * @param G  Grid geometry.
 * @param S  Current fluid state (source for reconstruction).
 * @param F  Output flux struct; filled with X1 and X2 fluxes.
 * @return   CFL-limited next timestep.
 */
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

/**
 * @brief Compute the Local Lax-Friedrichs (LLF) flux from left and right reconstructed states.
 *
 * @details Given the reconstructed left (Sl) and right (Sr) primitive state arrays,
 * this function:
 * 1. Shifts the left-state array so that Sl[i] refers to the left interface of zone i.
 * 2. Computes 4-velocities and magnetic 4-vectors for both states (get_state_vec()).
 * 3. Computes physical fluxes FL and FR and conserved variables UL, UR (prim_to_flux_vec()).
 * 4. Computes maximum and minimum characteristic speeds cmax, cmin (mhd_vchar()).
 * 5. Assembles the LLF flux:
 *    @f[ F = \frac{1}{2}\left(F_L + F_R - c_{\rm top}(U_R - U_L)\right) @f]
 *    where @f$c_{\rm top} = \max(|c_{\rm max}|, |c_{\rm min}|)@f$.
 *
 * @note The argument order (Sr before Sl) is intentional – the caller passes
 *       the "right" state first, matching the convention used in reconstruction.c.
 *
 * @param G     Grid geometry.
 * @param Sr    Right reconstructed state (before shifting).
 * @param Sl    Left reconstructed state (before shifting).
 * @param dir   Direction (1 = X1, 2 = X2).
 * @param loc   Grid centering location (FACE1 or FACE2).
 * @param flux  Output array for the assembled fluxes.
 * @param ctop  Output grid of maximum wave speeds (used for CFL).
 */
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

/**
 * @brief Apply the flux-CT (constrained transport) scheme to enforce @f$\nabla \cdot B = 0@f$.
 *
 * @details Implements the method of Toth (2000, JCP 161, 605).
 * 1. Computes the electromotive force (EMF) at each cell corner:
 *    @f[ \mathcal{E}_{i,j} = \frac{1}{4}\left(F^{X1}_{B2,i,j} + F^{X1}_{B2,i,j-1}
 *        - F^{X2}_{B1,i,j} - F^{X2}_{B1,i-1,j}\right) @f]
 * 2. Rewrites the B-component fluxes so that the discrete curl of the face-integrated
 *    EMF exactly cancels numerical divergence:
 *    - @f$F^{X1}_{B1} = 0@f$ (no X1 transport of B1 through X1 faces).
 *    - @f$F^{X1}_{B2}@f$ replaced by corner-averaged EMF.
 *    - @f$F^{X2}_{B1}@f$ replaced by negative corner-averaged EMF.
 *    - @f$F^{X2}_{B2} = 0@f$.
 *
 * @param F  Flux struct to modify in-place.
 */
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
