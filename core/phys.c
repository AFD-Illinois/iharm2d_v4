/**
 * @file phys.c
 * @brief GRMHD physics: stress-energy tensor, primitive-to-flux transforms, wave speeds, and source terms.
 *
 * @details This is the core physics module.  It implements the ideal GRMHD equations
 * in the 3+1 (ADM) split with a gamma-law equation of state P = (gamma-1)*u.
 *
 * **Key quantities:**
 * - Total enthalpy: @f$w = \rho + u + P + b^2 = \rho + u + (\gamma-1)u + b^2@f$
 * - Total pressure: @f$P_\text{tot} = P + b^2/2@f$
 * - Stress-energy tensor:
 *   @f[ T^\mu{}_\nu = w\, u^\mu u_\nu + P_\text{tot}\, \delta^\mu_\nu - b^\mu b_\nu @f]
 * - Conserved variables (all multiplied by @f$\sqrt{-g}@f$):
 *   - @f$D = \rho u^t \sqrt{-g}@f$
 *   - @f$\tau = -(T^t{}_t + D)@f$
 *   - @f$S_i = T^t{}_i \sqrt{-g}@f$
 *   - @f$B^i = B^i \sqrt{-g}@f$ (directly advected)
 * - Geometric source term: @f$\delta_T = T^\mu{}_\nu \Gamma^\nu_{\mu\lambda} \sqrt{-g}@f$
 *
 * **Magnetic 4-vector:**
 * @f[
 *   b^0 = B^i u_i / u^0, \quad b^i = (B^i + b^0 u^i) / u^0
 * @f]
 *
 * **Wave speeds:** The fast magnetosonic speed is computed from the combined
 * Alfvén+sound speed formula of Gammie, McKinney & Toth (2003).
 *
 * **Optional wind term:** When WIND_TERM=1, a small mass/energy injection is added
 * near the polar axis to prevent floors from repeatedly activating there.
 */

/*---------------------------------------------------------------------------------

  PHYS.C

  -Compute MHD stress-energy tensor
  -Obtain conserved variables and/or fluxes from primitives
  -Compute Lorentz factor
  -Compute 4 velocities from primitives
  -Calculate components of magnetosonic velocity
  -Calculate source terms
  -Compute bsq at given grid zone

---------------------------------------------------------------------------------*/

#include "decs.h"

// MHD stress-energy tensor with first index up, second index down.
// A factor of sqrt(4 pi) is absorbed into the definition of b.
/**
 * @brief Compute one row of the MHD stress-energy tensor T^dir_mu at zone (i,j).
 *
 * @details Evaluates
 * @f[ T^{\text{dir}}{}_\mu = w\, u^{\text{dir}} u_\mu + P_\text{tot}\,\delta^{\text{dir}}_\mu - b^{\text{dir}} b_\mu @f]
 * where @f$w = \rho + u + P + b^2@f$ and @f$P_\text{tot} = P + b^2/2@f$.
 * The 4-vectors ucon, ucov, bcon, bcov must already be populated in S (via get_state()).
 *
 * @param S    Fluid state with precomputed 4-vectors.
 * @param i    X1 zone index.
 * @param j    X2 zone index.
 * @param dir  Direction of the first (raised) index.
 * @param[out] mhd  Output array of length NDIM storing @f$T^{\text{dir}}{}_\mu@f$.
 */
inline void mhd_calc(struct FluidState *S, int i, int j, int dir, double *mhd)
{
  double u, pres, w, bsq, eta, ptot;

  u = S->P[UU][j][i];
  pres = (gam - 1.)*u;
  w = pres + S->P[RHO][j][i] + u;
  bsq = bsq_calc(S, i, j);
  eta = w + bsq;
  ptot = pres + 0.5*bsq;

  DLOOP1
    mhd[mu] = eta*S->ucon[dir][j][i]*S->ucov[mu][j][i] + ptot*delta(dir, mu) - S->bcon[dir][j][i]*S->bcov[mu][j][i];
}

/**
 * @brief Compute all primitive-to-flux contributions at a single zone (i, j).
 *
 * @details Fills flux[var][j][i] with the @c dir-direction flux of each variable,
 * multiplied by sqrt(-g) = G->gdet[loc][j][i].  For dir=0 this gives the conserved
 * variable vector U.  Electron entropy variables (if enabled) are advected passively
 * as specific quantities (KEL * rho * u^dir).
 *
 * @param G    Grid geometry.
 * @param S    Fluid state with precomputed ucon, ucov, bcon, bcov.
 * @param i    X1 zone index.
 * @param j    X2 zone index.
 * @param dir  Flux direction (0 = conserved, 1 = X1, 2 = X2).
 * @param loc  Grid centering location.
 * @param flux Output flux array (GridPrim).
 */
void prim_to_flux(struct GridGeom *G, struct FluidState *S, int i, int j, int dir, int loc, GridPrim flux)
{
  double mhd[NDIM];
  //Particle number flux
  flux[RHO][j][i] = S->P[RHO][j][i]*S->ucon[dir][j][i];
  
  mhd_calc(S, i, j, dir, mhd);
  // MHD stress-energy tensor with first index up and second index down
  flux[UU][j][i] = mhd[0] + flux[RHO][j][i];
  flux[U1][j][i] = mhd[1];
  flux[U2][j][i] = mhd[2];
  flux[U3][j][i] = mhd[3];
  //Dual of Maxwell Tensor
  flux[B1][j][i] = S->bcon[1][j][i]*S->ucon[dir][j][i] - S->bcon[dir][j][i]*S->ucon[1][j][i]; 
  flux[B2][j][i] = S->bcon[2][j][i]*S->ucon[dir][j][i] - S->bcon[dir][j][i]*S->ucon[2][j][i];
  flux[B3][j][i] = S->bcon[3][j][i]*S->ucon[dir][j][i] - S->bcon[dir][j][i]*S->ucon[3][j][i];

#if ELECTRONS
  for (int idx = KEL0; idx < NVAR; idx++)  
    flux[idx][j][i] = flux[RHO][j][i]*S->P[idx][j][i];
  flux[KTOT][j][i] = flux[RHO][j][i]*S->P[KTOT][j][i];
#endif
  
  // Multiply with gdet (Jacobian) to account for the coordinate system
  PLOOP flux[ip][j][i] *= G->gdet[loc][j][i];
}

/**
 * @brief OpenMP-parallelized prim_to_flux() over a rectangular zone range.
 *
 * @details Equivalent to calling prim_to_flux() for every zone in the range
 * [jstart,jstop] x [istart,istop] but parallelized with collapsed OpenMP loops.
 *
 * @param G      Grid geometry.
 * @param S      Fluid state with precomputed 4-vectors.
 * @param dir    Flux direction (0 = conserved, 1 = X1, 2 = X2).
 * @param loc    Grid centering location.
 * @param jstart First j index (active-zone relative, i.e., 0 = first active zone).
 * @param jstop  Last j index (inclusive).
 * @param istart First i index.
 * @param istop  Last i index (inclusive).
 * @param flux   Output flux array.
 */
void prim_to_flux_vec(struct GridGeom *G, struct FluidState *S, int dir, int loc, int jstart, int jstop, int istart, int istop, GridPrim flux)
{
#pragma omp parallel
{
#pragma omp for collapse(2) nowait
  ZSLOOP(jstart, jstop, istart, istop)
  {
    double mhd[NDIM];

    flux[RHO][j][i] = S->P[RHO][j][i] * S->ucon[dir][j][i] * G->gdet[loc][j][i];
    mhd_calc(S, i, j, dir, mhd);
    // MHD stress-energy tensor w/ first index up, second index down
    flux[UU][j][i] = mhd[0] * G->gdet[loc][j][i] + flux[RHO][j][i];
    flux[U1][j][i] = mhd[1] * G->gdet[loc][j][i];
    flux[U2][j][i] = mhd[2] * G->gdet[loc][j][i];
    flux[U3][j][i] = mhd[3] * G->gdet[loc][j][i];
  }

#pragma omp for collapse(2) nowait
  ZSLOOP(jstart, jstop, istart, istop)
  {
    // Dual of Maxwell tensor
    flux[B1][j][i] = (S->bcon[1][j][i] * S->ucon[dir][j][i] - S->bcon[dir][j][i] * S->ucon[1][j][i]) * G->gdet[loc][j][i];
    flux[B2][j][i] = (S->bcon[2][j][i] * S->ucon[dir][j][i] - S->bcon[dir][j][i] * S->ucon[2][j][i]) * G->gdet[loc][j][i];
    flux[B3][j][i] = (S->bcon[3][j][i] * S->ucon[dir][j][i] - S->bcon[dir][j][i] * S->ucon[3][j][i]) * G->gdet[loc][j][i];

  }

#if ELECTRONS
#pragma omp for collapse(2)
  ZSLOOP(jstart, jstop, istart, istop)
  {
    // RHO already includes a factor of gdet!
    for (int idx = KEL0; idx < NVAR ; idx++)
      flux[idx][j][i] = flux[RHO][j][i]*S->P[idx][j][i];
    flux[KTOT][j][i] = flux[RHO][j][i]*S->P[KTOT][j][i];
  }
#endif
}
}

/**
 * @brief Compute the contravariant magnetic 4-vector b^mu from primitive B fields.
 *
 * @details Uses the ideal MHD relation and the already-computed 4-velocity u^mu:
 * @f[
 *   b^0 = B^i u_i / u^0, \quad b^k = (B^k + b^0 u^k) / u^0
 * @f]
 * The 4-velocity ucov must already be set in S.
 *
 * @param S  Fluid state (reads P[B1-B3], ucon, ucov; writes bcon).
 * @param i  X1 zone index.
 * @param j  X2 zone index.
 */
inline void bcon_calc(struct FluidState *S, int i, int j)
{
  S->bcon[0][j][i] = S->P[B1][j][i]*S->ucov[1][j][i] + S->P[B2][j][i]*S->ucov[2][j][i] + S->P[B3][j][i]*S->ucov[3][j][i];
  for (int mu = 1; mu < 4; mu++)
    S->bcon[mu][j][i] = (S->P[B1-1+mu][j][i] + S->bcon[0][j][i]*S->ucon[mu][j][i])/S->ucon[0][j][i];
}

/**
 * @brief Compute the Lorentz factor gamma with respect to the normal observer.
 *
 * @details Evaluates @f$\gamma = \sqrt{1 + q^2}@f$ where
 * @f$q^2 = g_{ij} v^i v^j@f$ (spatial metric contracted with 3-velocity).
 * If @f$q^2 < 0@f$ due to numerical error, applies a floor of 1e-10 in DEBUG mode.
 *
 * @param G   Grid geometry.
 * @param S   Fluid state (reads P[U1-U3]).
 * @param i   X1 zone index.
 * @param j   X2 zone index.
 * @param loc Grid centering location.
 * @return    Lorentz factor gamma >= 1.
 */
inline double mhd_gamma_calc(struct GridGeom *G, struct FluidState *S, int i, int j, int loc)
{
  double qsq = G->gcov[loc][1][1][j][i]*S->P[U1][j][i]*S->P[U1][j][i]
      + G->gcov[loc][2][2][j][i]*S->P[U2][j][i]*S->P[U2][j][i]
      + G->gcov[loc][3][3][j][i]*S->P[U3][j][i]*S->P[U3][j][i]
      + 2.*(G->gcov[loc][1][2][j][i]*S->P[U1][j][i]*S->P[U2][j][i]
          + G->gcov[loc][1][3][j][i]*S->P[U1][j][i]*S->P[U3][j][i]
          + G->gcov[loc][2][3][j][i]*S->P[U2][j][i]*S->P[U3][j][i]);

#if DEBUG
  if (qsq < 0.)
  {
    if (fabs(qsq) > 1.E-10) // Then assume not just machine precision
    {
      fprintf(stderr, "gamma_calc():  failed: [%i %i] qsq = %28.18e \n", i, j, qsq);
      fprintf(stderr, "v[1-3] = %28.18e %28.18e %28.18e  \n", S->P[U1][j][i], S->P[U2][j][i], S->P[U3][j][i]);
      return 1.0;
    } 
    else
      qsq = 1.E-10; // Set floor
  }
#endif
  
  return sqrt(1. + qsq);
}

/**
 * @brief Compute the contravariant 4-velocity u^mu from primitive velocities.
 *
 * @details The primitive velocities U^i = v^i are the 3-velocity components
 * in the coordinate frame (= gamma * u^i - gamma * alpha * g^0i).
 * This function recovers the full u^mu using:
 * @f[
 *   u^0 = \gamma / \alpha, \quad u^i = V^i - \gamma \alpha g^{0i}
 * @f]
 *
 * @param G   Grid geometry.
 * @param S   Fluid state (reads P[U1-U3]; writes ucon).
 * @param i   X1 zone index.
 * @param j   X2 zone index.
 * @param loc Grid centering location.
 */
inline void ucon_calc(struct GridGeom *G, struct FluidState *S, int i, int j, int loc)
{
  double gamma = mhd_gamma_calc(G, S, i, j, loc);
  double alpha = G->lapse[loc][j][i];
  S->ucon[0][j][i] = gamma/alpha;
  for (int mu = 1; mu < NDIM; mu++)
    S->ucon[mu][j][i] = S->P[U1+mu-1][j][i] - gamma*alpha*G->gcon[loc][0][mu][j][i];
}

// Calculate ucon, ucov, bcon, bcov from primitive variables
inline void get_state(struct GridGeom *G, struct FluidState *S, int i, int j, int loc)
{
  ucon_calc(G, S, i, j, loc);
  lower_grid(S->ucon, S->ucov, G, i, j, loc);
  bcon_calc(S, i, j);
  lower_grid(S->bcon, S->bcov, G, i, j, loc);
}

// Calculate ucon, ucov, bcon, bcov from primitive variables, over given range
// Note same range convention as ZSLOOP and other *_vec functions
/**
 * @brief OpenMP-parallelized get_state() over a rectangular zone range.
 *
 * @param G      Grid geometry.
 * @param S      Fluid state.
 * @param loc    Grid centering location.
 * @param jstart First j index (active-zone relative).
 * @param jstop  Last j index (inclusive).
 * @param istart First i index.
 * @param istop  Last i index (inclusive).
 */
void get_state_vec(struct GridGeom *G, struct FluidState *S, int loc, int jstart, int jstop, int istart, int istop)
{
#pragma omp parallel
  {
#pragma omp for collapse(2)
    ZSLOOP(jstart, jstop, istart, istop)
      ucon_calc(G, S, i, j, loc);

#pragma omp for collapse(2)
    ZSLOOP(jstart, jstop, istart, istop)
      lower_grid(S->ucon, S->ucov, G, i, j, loc);

#pragma omp for collapse(2)
    ZSLOOP(jstart, jstop, istart, istop)
      bcon_calc(S, i, j);

#pragma omp for collapse(2)
    ZSLOOP(jstart, jstop, istart, istop)
      lower_grid(S->bcon, S->bcov, G, i, j, loc);
  }
}

/**
 * @brief Compute the fast magnetosonic wave speeds at zone (i, j) in direction dir.
 *
 * @details Follows Gammie, McKinney & Toth (2003, ApJ 589, 444), Section 4.
 * The combined fast magnetosonic + Alfvén speed squared is:
 * @f[ c_\text{ms}^2 = c_s^2 + v_A^2 - c_s^2 v_A^2 @f]
 * where @f$c_s^2 = \gamma(\gamma-1)u / (\rho + \gamma u)@f$ and
 * @f$v_A^2 = b^2 / (\rho + u + P + b^2)@f$.
 * The wave speeds are then roots of a quadratic in the observer frame.
 * If the discriminant is negative (numerical artifact) it is clamped to zero.
 *
 * @param G    Grid geometry.
 * @param S    Fluid state with precomputed 4-vectors.
 * @param i    X1 zone index.
 * @param j    X2 zone index.
 * @param loc  Grid centering location.
 * @param dir  Wave propagation direction (1 = X1, 2 = X2).
 * @param[out] cmax  Grid array; receives the forward (max) wave speed at [j][i].
 * @param[out] cmin  Grid array; receives the reverse (min) wave speed at [j][i].
 */
inline void mhd_vchar(struct GridGeom *G, struct FluidState *S, int i, int j, int loc, int dir, GridDouble cmax, GridDouble cmin)
{
  double discr, vp, vm, bsq, ee, ef, va2, cs2, cms2, rho, u;
  double Acov[NDIM], Bcov[NDIM], Acon[NDIM], Bcon[NDIM];
  double Asq, Bsq, Au, Bu, AB, Au2, Bu2, AuBu, A, B, C;

  DLOOP1
    Acov[mu] = 0.;
  Acov[dir] = 1.;

  DLOOP1
    Bcov[mu] = 0.;
  Bcov[0] = 1.;

  DLOOP1
  {
    Acon[mu] = 0.;
    Bcon[mu] = 0.;
  }

  DLOOP2
  {
    Acon[mu] += G->gcon[loc][mu][nu][j][i]*Acov[nu];
    Bcon[mu] += G->gcon[loc][mu][nu][j][i]*Bcov[nu];
  }

  // Find fast magnetosonic velocity
  bsq = bsq_calc(S, i, j);
  rho = fabs(S->P[RHO][j][i]);
  u = fabs(S->P[UU][j][i]);
  ef = rho + gam*u;
  ee = bsq + ef;
  va2 = bsq/ee;
  cs2 = gam*(gam - 1.)*u/ef;

  cms2 = cs2 + va2 - cs2*va2;

  cms2 = (cms2 < 0) ? SMALL : cms2;
  cms2 = (cms2 > 1) ? 1 : cms2;

  // Require that speed of wave measured by observer q->ucon is cms2
  Asq = dot(Acon, Acov);
  Bsq = dot(Bcon, Bcov);
  Au = Bu = 0.;
  DLOOP1 
  {
    Au += Acov[mu]*S->ucon[mu][j][i];
    Bu += Bcov[mu]*S->ucon[mu][j][i];
  }
  AB = dot(Acon, Bcov);
  Au2 = Au*Au;
  Bu2 = Bu*Bu;
  AuBu = Au*Bu;

  A = Bu2 - (Bsq + Bu2)*cms2;
  B = 2.*(AuBu - (AB + AuBu)*cms2);
  C = Au2 - (Asq + Au2)*cms2;

  discr = B*B - 4.*A*C;
  discr = (discr < 0.) ? 0. : discr;
  discr = sqrt(discr);

  vp = -(-B + discr)/(2.*A);
  vm = -(-B - discr)/(2.*A);

  cmax[j][i] = (vp > vm) ? vp : vm;
  cmin[j][i] = (vp > vm) ? vm : vp;
}

/**
 * @brief Compute the geometric source terms for the GRMHD equations.
 *
 * @details The source terms arise from the Christoffel connection coefficients
 * in curved spacetime.  For the energy-momentum equations:
 * @f[ \delta U^\nu = T^\mu{}_\lambda \Gamma^\lambda_{\mu\nu} \sqrt{-g} \, \Delta V @f]
 * The baryon number equation has no source.  Magnetic field equations are also sourceless.
 * If WIND_TERM is enabled, an additional mass/energy injection rate is added near
 * the coordinate axis to stabilize the funnel region.
 *
 * @param G    Grid geometry (provides Christoffel symbols and gdet).
 * @param S    Fluid state with precomputed 4-vectors.
 * @param dU   Output source term array (GridPrim pointer); accumulated in-place by advance_fluid().
 */
inline void get_fluid_source(struct GridGeom *G, struct FluidState *S, GridPrim *dU)
{
#if WIND_TERM
  static struct FluidState *dS;
  static int firstc = 1;
  if (firstc)
  {
    dS = calloc(1,sizeof(struct FluidState)); 
    firstc = 0;
  }
#endif

#pragma omp parallel for collapse(2)
  ZLOOP 
  {
    double mhd[NDIM][NDIM];
    DLOOP1 mhd_calc(S, i, j, mu, mhd[mu]);
    // Contract mhd stress tensor with connection
    PLOOP (*dU)[ip][j][i] = 0.;
    DLOOP2
      for (int gam = 0; gam < NDIM; gam++)
        (*dU)[UU+gam][j][i] += mhd[mu][nu]*G->conn[nu][gam][mu][j][i];
    PLOOP (*dU)[ip][j][i] *= G->gdet[CENT][j][i];
  }

  // Add a small "wind" source term in RHO,UU
#if WIND_TERM
#pragma omp parallel for simd
  ZLOOP
  {
    /* need coordinates to evaluate particle addtn rate */
    double X[NDIM];
    coord(i, j, CENT, X);
    double r, th;
    bl_coord(X, &r, &th);
    double cth = cos(th);

    /* here is the rate at which we're adding particles */
    /* this function is designed to concentrate effect in the funnel in black hole evolutions */
    double drhopdt = 2.e-4*cth*cth*cth*cth/pow(1. + r*r,2);
    dS->P[RHO][j][i] = drhopdt ;
    double Tp = 10. ;  /* temp, in units of c^2, of new plasma */
    dS->P[UU][j][i] = drhopdt*Tp*3. ;
    /* Leave P[U{1,2,3}]=0 to add in particles in normal observer frame */
    /* Likewise leave P[BN]=0 */
  }

  /* add in plasma to the T^t_a component of the stress-energy tensor */
  /* notice that U already contains a factor of sqrt{-g} */
  get_state_vec(G, dS, CENT, 0, N2-1, 0, N1-1);
  prim_to_flux_vec(G, dS, 0, CENT, 0, N2-1, 0, N1-1, dS->U);

#pragma omp parallel for simd collapse(2)
  PLOOP 
    ZLOOP
      (*dU)[ip][j][i] += dS->U[ip][j][i];
#endif
}

/**
 * @brief Compute @f$b^\mu b_\mu@f$ (twice the magnetic pressure) at zone (i, j).
 *
 * @details Sums @f$b^\mu b_\mu = b^0 b_0 + b^1 b_1 + b^2 b_2 + b^3 b_3@f$.
 * The 4-vectors bcon and bcov must already be set in S (via get_state()).
 * Returns max(bsq, SMALL) to prevent division by zero in floor calculations.
 *
 * @param S  Fluid state with precomputed bcon, bcov.
 * @param i  X1 zone index.
 * @param j  X2 zone index.
 * @return   @f$b^\mu b_\mu \ge \text{SMALL}@f$.
 */
inline double bsq_calc(struct FluidState *S, int i, int j)
{
  double bsq = 0.;
  DLOOP1
    bsq += S->bcon[mu][j][i]*S->bcov[mu][j][i];
  
  return MY_MAX(bsq, SMALL);
}
