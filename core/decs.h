/**
 * @file decs.h
 * @brief Global declarations: macros, compile-time parameters, type definitions, data structures, and function prototypes.
 *
 * @details This is the primary header for iharm2d_v4. It is included by every translation unit.
 * It defines:
 * - Compile-time numerical and physical parameters (floors, ceilings, grid constants).
 * - Primitive and conserved variable index macros (RHO, UU, U1-U3, B1-B3, optional electron vars).
 * - Grid-centering location codes (FACE1, FACE2, CORN, CENT).
 * - Boundary condition codes (OUTFLOW, PERIODIC, POLAR, USER).
 * - Metric codes (MINKOWSKI, MKS).
 * - Reconstruction algorithm codes (LINEAR, WENO, MP5).
 * - Diagnostic and failure-mode codes.
 * - Timer codes for performance profiling.
 * - Type aliases for grid arrays (GridInt, GridDouble, GridVector, GridPrim).
 * - The three core data structures: GridGeom, FluidState, FluidFlux.
 * - Loop macros (ZLOOP, ILOOP, JLOOP, PLOOP, DLOOP1, DLOOP2, etc.).
 * - External declarations for all global variables defined in defs.h.
 * - Forward declarations of every function in the codebase.
 *
 * @note The problem-specific header @c parameters.h is included here, so all
 *       compile-time grid sizes (N1TOT, N2TOT), metric choice, reconstruction
 *       choice, and boundary condition choices are resolved at compile time.
 */

#pragma once

// Import C libraries
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <errno.h> //Errors for syscalls

#include <omp.h>

// Include problem-specific header (parameter) file
#include "parameters.h"

// Constants
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923132169164
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.4142135623730950488016887242
#endif

/*---------------------------------------------------------------------------------
                            COMPILE-TIME PARAMETERS
---------------------------------------------------------------------------------*/

/** @defgroup compile_params Compile-time Parameters
 *  @brief Constants and flags resolved at compile time.
 *  @{
 */

/** @brief Code version string written to every dump header. */
#define VERSION "iharm2d_v4-alpha-1.0"

/** @brief Number of active zones in the X1 direction (= N1TOT from parameters.h). */
#define N1       (N1TOT)
/** @brief Number of active zones in the X2 direction (= N2TOT from parameters.h). */
#define N2       (N2TOT)
/** @brief Maximum of N1 and N2; used to size temporary 1D arrays. */
#define NMAX     (N1 > N2 ? N1 : N2)

/** @brief Total number of spacetime dimensions (always 4). */
#define NDIM       (4)
/** @brief Number of ghost zones on each boundary (3 required for WENO/MP5). */
#define NG         (3)

/** @defgroup floors Floor and Ceiling Values
 *  @brief Absolute lower bounds on density and internal energy used in fixup.c.
 *  @{
 */
/** @brief Absolute minimum rest-mass density (applied after all other floors). */
#define RHOMINLIMIT (1.e-20)
/** @brief Absolute minimum internal energy density. */
#define UUMINLIMIT  (1.e-20)
/** @brief Geometric floor coefficient for rest-mass density: rhoflr = RHOMIN * r^-2 / (1 + r/10). */
#define RHOMIN      (1.e-6)
/** @brief Geometric floor coefficient for internal energy density. */
#define UUMIN       (1.e-8)
/** @} */

/** @brief A numerically small value (~ machine epsilon) used to prevent division by zero. */
#define SMALL (1.e-20)
/** @brief Finite-difference step size used when computing connection coefficients numerically. */
#define DELTA 1.e-5

/** @brief Maximum allowed ratio u/rho (temperature ceiling). Overridable at compile time. */
#ifndef UORHOMAX
#define UORHOMAX (100.)
#endif
/** @brief Maximum allowed ratio b^2/rho (magnetization ceiling for density floor). Overridable. */
#ifndef BSQORHOMAX
#define BSQORHOMAX (100.)
#endif
/** @brief Maximum allowed ratio b^2/u. Defaults to BSQORHOMAX * UORHOMAX. */
#ifndef BSQOUMAX
#define BSQOUMAX (BSQORHOMAX * UORHOMAX)
#endif
/** @brief Enable (1) or disable (0) the "wind" mass-injection source term.
 *  Disabled by default; test problems always set this to 0. */
#ifndef WIND_TERM
#define WIND_TERM 0
#endif

/** @brief Maximum allowed Lorentz factor gamma for the fluid. */
#define GAMMAMAX (50.)

/** @brief Maximum fractional increase in timestep per step (safety factor for adaptive CFL). */
#define SAFE  (1.3)

/** @brief If 1, shift the polar axis coordinate slightly off the coordinate singularity. */
#define COORDSINGFIX 1
/** @brief Half-width of the coordinate singularity avoidance zone near the pole. */
#define SINGSMALL (1.E-20)

/** @brief Enable (1) extra debug output and sanity checks. Set via -DDEBUG=1 at compile time. */
#ifndef DEBUG
#define DEBUG 0
#endif
/** @brief Enable (1) performance timer instrumentation. */
#ifndef TIMERS
#define TIMERS 1
#endif
/** @brief If 1, fix the timestep to the initial value (no adaptive CFL). */
#ifndef STATIC_TIMESTEP
#define STATIC_TIMESTEP 0
#endif

/** @brief Default length for character string buffers. */
#define STRLEN (2048)

/** @} */ // end compile_params

/*---------------------------------------------------------------------------------
             Reconstruction Algorithm Codes
---------------------------------------------------------------------------------*/
/** @defgroup recon_codes Reconstruction Algorithm Codes
 *  @{
 */
/** @brief Linear reconstruction with Monotonized Central (MC) slope limiter. */
#define LINEAR (0)
/** @brief Piecewise Parabolic Method (placeholder; not fully implemented). */
#define PPM    (1)
/** @brief 5th-order Weighted Essentially Non-Oscillatory (WENO) reconstruction. */
#define WENO   (2)
/** @brief 5th-order Monotonicity Preserving (MP5) reconstruction. */
#define MP5    (3)
/** @} */

/*---------------------------------------------------------------------------------
             Primitive and Conserved Variable Index Macros
---------------------------------------------------------------------------------*/
/** @defgroup prim_indices Primitive Variable Indices
 *  @brief Index macros for the NVAR-dimensional primitive (and conserved) variable array.
 *
 *  Primitives: P = [rho, u, U^1, U^2, U^3, B^1, B^2, B^3], possibly followed by electron entropy variables.
 *  @{
 */
#define RHO (0) /**< Rest-mass density rho. */
#define UU  (1) /**< Internal energy density u. */
#define U1  (2) /**< Contravariant velocity component v^1 (= gamma * u^1 - gamma * alpha * g^01). */
#define U2  (3) /**< Contravariant velocity component v^2. */
#define U3  (4) /**< Contravariant velocity component v^3. */
#define B1  (5) /**< Contravariant magnetic field B^1 (coordinate frame). */
#define B2  (6) /**< Contravariant magnetic field B^2. */
#define B3  (7) /**< Contravariant magnetic field B^3. */
#if ELECTRONS
/** @brief Total fluid entropy K_tot = (gamma-1) * u / rho^gamma (only with ELECTRONS). */
#define KTOT (8)
#if ALLMODELS
/** @brief Kawazura (2019) electron entropy variable. */
#define KEL0 (9)
/** @brief Werner et al. (2018) electron entropy variable. */
#define KEL1 (10)
/** @brief Rowan et al. (2017) electron entropy variable. */
#define KEL2 (11)
/** @brief Sharma et al. (2007) electron entropy variable. */
#define KEL3 (12)
/** @brief Total number of primitive variables (with all four electron models). */
#define NVAR (13)
#else
/** @brief Electron entropy for the single active heating model (ALLMODELS=0). */
#define KEL0  (9)
/** @brief Total number of primitive variables (single electron model). */
#define NVAR (10)
#endif
#else
/** @brief Total number of primitive variables (pure MHD, no electrons). */
#define NVAR (8)
#endif
/** @} */

/*---------------------------------------------------------------------------------
             Grid Centering Location Codes
---------------------------------------------------------------------------------*/
/** @defgroup loc_codes Grid Centering Location Codes
 *  @brief Used as the @c loc argument to coord(), get_state(), prim_to_flux(), etc.
 *
 *  The staggered-grid layout within a cell is:
 *  @verbatim
 *  -----------------------
 *  |                     |
 *  |                     |
 *  |FACE1   CENT         |
 *  |                     |
 *  |CORN    FACE2        |
 *  -----------------------
 *  @endverbatim
 *  @{
 */
#define FACE1 (0) /**< Left (X1-direction) face center. */
#define FACE2 (1) /**< Bottom (X2-direction) face center. */
#define CORN  (2) /**< Corner (lower-left vertex of cell). */
#define CENT  (3) /**< Zone center. */
#define NPG   (4) /**< Total number of grid positions stored per zone. */
/** @} */

/*---------------------------------------------------------------------------------
             Boundary Condition Codes
---------------------------------------------------------------------------------*/
/** @defgroup bc_codes Boundary Condition Codes
 *  @brief Values for X1L_BOUND, X1R_BOUND, X2L_BOUND, X2R_BOUND in parameters.h.
 *  @{
 */
#define OUTFLOW  (0) /**< Outflow (copy outermost active zone into ghost zones, rescale B). */
#define PERIODIC (1) /**< Periodic (wrap across opposite boundary). */
#define POLAR    (2) /**< Polar reflection (reflect across axis, flip U2 and B2 signs). */
#define USER     (3) /**< Problem-specific boundary condition provided by bound_gas_prob_x1r(). */
/** @} */

/*---------------------------------------------------------------------------------
             Metric Codes
---------------------------------------------------------------------------------*/
/** @defgroup metric_codes Metric Codes
 *  @brief Values for the METRIC compile-time flag.
 *  @{
 */
#define MINKOWSKI (0) /**< Flat Minkowski (Cartesian) spacetime. */
#define MKS       (1) /**< Modified Kerr-Schild coordinates for Kerr black hole. */
/** @} */

/*---------------------------------------------------------------------------------
             IO and Diagnostic Mode Codes
---------------------------------------------------------------------------------*/
/** @brief Normal (numbered) dump or restart file. */
#define IO_REGULAR (1)
/** @brief Abort dump/restart (written when 'abort' file is detected). */
#define IO_ABORT   (2)

/** @defgroup diag_codes Diagnostic Call Codes
 *  @brief Passed to diag() to control output behavior.
 *  @{
 */
#define DIAG_INIT  (0) /**< Called once at startup: open log file, write initial dump. */
#define DIAG_DUMP  (1) /**< Called when a full dump is due (t >= tdump). */
#define DIAG_LOG   (2) /**< Called when a log entry is due (t >= tlog). */
#define DIAG_FINAL (3) /**< Called at the end of the simulation. */
#define DIAG_ABORT (4) /**< Called when an 'abort' file is detected. */
/** @} */

/*---------------------------------------------------------------------------------
             Failure Mode Codes
---------------------------------------------------------------------------------*/
/** @defgroup fail_codes Failure Mode Codes
 *  @brief Values stored in pflag[] and fail_save[] to record the type of inversion failure.
 *  @{
 */
#define FAIL_UTOPRIM     (0) /**< Conservative-to-primitive inversion failed. */
#define FAIL_VCHAR_DISCR (1) /**< Negative discriminant in magnetosonic speed calculation. */
#define FAIL_COEFF_NEG   (2) /**< Negative coefficient in U_to_P Newton-Raphson step. */
#define FAIL_COEFF_SUP   (3) /**< Coefficient exceeded maximum in U_to_P. */
#define FAIL_GAMMA       (4) /**< Lorentz factor exceeded GAMMAMAX. */
#define FAIL_METRIC      (5) /**< Metric inversion failed. */
/** @} */

/*---------------------------------------------------------------------------------
             Timer Codes
---------------------------------------------------------------------------------*/
/** @defgroup timer_codes Performance Timer Codes
 *  @brief Used with timer_start() / timer_stop() to instrument code sections.
 *  @{
 */
#define TIMER_RECON    (1)  /**< Spatial reconstruction (reconstruct()). */
#define TIMER_LR_TO_F  (2)  /**< Left/right states to LLF flux (lr_to_flux()). */
#define TIMER_CMAX     (3)  /**< CFL timestep calculation (ndt_min()). */
#define TIMER_FLUX_CT  (4)  /**< Constrained transport EMF update (flux_ct()). */
#define TIMER_UPDATE_U (5)  /**< Conservative variable update (advance_fluid interior). */
#define TIMER_U_TO_P   (6)  /**< Conservative-to-primitive inversion (U_to_P()). */
#define TIMER_FIXUP    (7)  /**< Floor/ceiling fixups (fixup(), fixup_utoprim()). */
#define TIMER_BOUND    (8)  /**< Boundary condition application (set_bounds()). */
#define TIMER_BOUND_COMMS (9)  /**< Boundary communication (unused without MPI; reserved). */
#define TIMER_DIAG        (10) /**< Diagnostic computation (diag()). */
#define TIMER_LR_STATE    (11) /**< get_state() calls inside lr_to_flux(). */
#define TIMER_LR_PTOF     (12) /**< prim_to_flux() calls inside lr_to_flux(). */
#define TIMER_LR_VCHAR    (13) /**< mhd_vchar() calls inside lr_to_flux(). */
#define TIMER_LR_CMAX     (14) /**< Combined cmax/cmin computation in lr_to_flux(). */
#define TIMER_LR_FLUX     (15) /**< LLF flux assembly in lr_to_flux(). */
#define TIMER_IO          (16) /**< Dump file I/O (dump_backend()). */
#define TIMER_RESTART     (17) /**< Restart file I/O (restart_write_backend()). */
#define TIMER_CURRENT     (18) /**< 4-current calculation (current_calc()). */
#define TIMER_ALL         (19) /**< Total wallclock time per step (wraps entire step()). */
#if ELECTRONS
#define TIMER_ELECTRON_FIXUP (21) /**< Electron entropy fixup (fixup_electrons()). */
#define TIMER_ELECTRON_HEAT  (22) /**< Electron heating (heat_electrons()). */
#define NUM_TIMERS           (23) /**< Total number of timers when ELECTRONS is enabled. */
#else
#define NUM_TIMERS     (20) /**< Total number of timers (pure MHD). */
#endif
/** @} */

/*---------------------------------------------------------------------------------
                                  GLOBAL TYPES
---------------------------------------------------------------------------------*/

/** @defgroup grid_types Grid Array Type Aliases
 *  @brief Typedef'd fixed-size 2D and 3D array types for the full grid (including ghost zones).
 *
 *  All grid arrays include NG ghost zones on each side.
 *  The index convention is [j][i] = [X2][X1] (row = j, column = i).
 *  @{
 */
/** @brief 2D integer array over the full grid (active + ghost zones). */
typedef int    GridInt[N2+2*NG][N1+2*NG];
/** @brief 2D double-precision array over the full grid. */
typedef double GridDouble[N2+2*NG][N1+2*NG];
/** @brief 3D array storing a 4-vector component at each grid point: [mu][j][i]. */
typedef double GridVector[NDIM][N2+2*NG][N1+2*NG];
/** @brief 4D array storing all NVAR primitive (or conserved) variables at each grid point: [var][j][i]. */
typedef double GridPrim[NVAR][N2+2*NG][N1+2*NG];
/** @} */

/**
 * @brief Stores precomputed metric quantities at all four staggered grid locations.
 *
 * @details All metric tensors are computed once at startup by set_grid() and stored
 * in this struct so that expensive metric evaluations are not repeated each step.
 * The leading index of gcov, gcon, gdet, and lapse is the grid location (FACE1, FACE2, CORN, CENT).
 * Connection coefficients Gamma^lambda_mu_nu are only needed at CENT.
 */
struct GridGeom {
  /** @brief Covariant metric tensor g_mu_nu at each grid location: [loc][mu][nu][j][i]. */
  double gcov[NPG][NDIM][NDIM][N2+2*NG][N1+2*NG];
  /** @brief Contravariant metric tensor g^mu^nu at each grid location: [loc][mu][nu][j][i]. */
  double gcon[NPG][NDIM][NDIM][N2+2*NG][N1+2*NG];
  /** @brief Metric determinant sqrt(-g) at each grid location: [loc][j][i]. */
  double gdet[NPG][N2+2*NG][N1+2*NG];
  /** @brief Lapse function alpha = (-g^00)^{-1/2} at each grid location: [loc][j][i]. */
  double lapse[NPG][N2+2*NG][N1+2*NG];
  /** @brief Connection coefficients Gamma^lam_nu_mu at zone centers: [lam][nu][mu][j][i]. */
  double conn[NDIM][NDIM][NDIM][N2+2*NG][N1+2*NG];
};

/**
 * @brief Stores the full fluid state at a given time level.
 *
 * @details This struct holds:
 * - The NVAR-component primitive variable vector P.
 * - The NVAR-component conserved variable vector U (= sqrt(-g) * T^t_mu, D, B^i).
 * - Auxiliary 4-vectors derived from P: ucon, ucov, bcon, bcov, jcon.
 *
 * Primitives P[]:
 * - P[RHO]: rest-mass density rho
 * - P[UU]:  internal energy density u
 * - P[U1-U3]: three-velocity components (Gamma * v^i - Gamma * alpha * g^0i)
 * - P[B1-B3]: coordinate-frame magnetic field components B^i
 * - P[KTOT], P[KEL0...]: electron entropy variables (only with ELECTRONS)
 *
 * Conserved variables U[] (= fluxes in the time direction, times sqrt(-g)):
 * - U[RHO]: baryon number D = sqrt(-g) * rho * u^t
 * - U[UU]:  energy tau = sqrt(-g) * (T^t_t - rho * u^t)  [convention: -T^t_t - D]
 * - U[U1-U3]: momentum S_i = sqrt(-g) * T^t_i
 * - U[B1-B3]: conserved magnetic field sqrt(-g) * B^i
 */
struct FluidState {
  GridPrim P;      /**< Primitive variables [var][j][i]. */
  GridPrim U;      /**< Conserved variables [var][j][i]. */
  GridVector ucon; /**< Contravariant 4-velocity u^mu [mu][j][i]. */
  GridVector ucov; /**< Covariant 4-velocity u_mu [mu][j][i]. */
  GridVector bcon; /**< Contravariant magnetic 4-vector b^mu [mu][j][i]. */
  GridVector bcov; /**< Covariant magnetic 4-vector b_mu [mu][j][i]. */
  GridVector jcon; /**< Contravariant 4-current j^mu [mu][j][i] (written by current_calc()). */
};

/**
 * @brief Stores numerical fluxes at cell faces in both coordinate directions.
 *
 * @details flux_ct() modifies these fluxes to enforce the divergence-free constraint
 * on the magnetic field via constrained transport.
 */
struct FluidFlux {
  GridPrim X1; /**< Fluxes through the X1 (left/right) faces: [var][j][i]. */
  GridPrim X2; /**< Fluxes through the X2 (top/bottom) faces: [var][j][i]. */
};

/** @brief Per-zone flag recording which U_to_P inversion failed (0 = success). */
extern GridInt pflag;
/** @brief Saved copy of pflag before fixup_utoprim() so dump files record the raw state. */
extern GridInt fail_save;
/** @brief Per-zone bit-mask recording which floor/ceiling conditions were hit. */
extern GridInt fflag;

/*---------------------------------------------------------------------------------
                              GLOBAL VARIABLES SECTION
---------------------------------------------------------------------------------*/

/** @defgroup global_vars Global Variables
 *  @brief All global variables are declared here (extern) and defined in defs.h.
 *  @{
 */

/** @brief Black hole spin parameter a (|a| <= 1). */
extern double a;
/** @brief Adiabatic index of the gas (gam = 4/3 for relativistic, 5/3 for non-relativistic). */
extern double gam;
/** @brief Boyer-Lindquist radius of the event horizon r_+ = 1 + sqrt(1 - a^2). */
extern double Rhor;
/** @brief Ion-to-electron temperature ratio T_p/T_e (used in electron heating). */
extern double tp_over_te;

/** @brief Inner radial boundary in Boyer-Lindquist coordinates (set to place 5 zones inside EH). */
extern double Rin;
/** @brief Outer radial boundary in Boyer-Lindquist coordinates. */
extern double Rout;
/** @brief MKS concentration parameter: hslope = 1 gives uniform spacing in theta, < 1 concentrates toward equator. */
extern double hslope;
/** @brief Courant safety factor (CFL number, typically 0.9). */
extern double cour;
/** @brief Cell volume element dV = dx[1] * dx[2] (zone area in 2D). */
extern double dV;
/** @brief Coordinate spacing: dx[1] = Delta X1, dx[2] = Delta X2. */
extern double dx[NDIM];
/** @brief Starting coordinate values of the active grid: startx[1] = log(Rin) for MKS. */
extern double startx[NDIM];
/** @brief Cartesian extent of the domain for MINKOWSKI metric. */
extern double x1Min, x1Max, x2Min, x2Max;
/** @brief Current timestep size (adaptive CFL). */
extern double dt;
/** @brief Light-crossing time of the smallest zone (informational). */
extern double dt_light;
/** @brief Current simulation time. */
extern double t;
/** @brief Final simulation time (run ends when t >= tf). */
extern double tf;
/** @brief Total number of completed timesteps. */
extern int nstep;
/** @brief 1 if simulation was initialized from a restart file, 0 otherwise. */
extern int is_restart;

/** @brief Cadence for dump output (time interval between dumps). */
extern double DTd;
/** @brief Cadence for flux output (unused in current version; placeholder). */
extern double DTf;
/** @brief Cadence for log output (time interval between log entries). */
extern double DTl;
/** @brief Cadence for restart writes (write restart every DTr steps). */
extern int DTr;
/** @brief Cadence for performance reports (print performance every DTp steps). */
extern int DTp;
/** @brief Sequential dump file counter (used to generate dump filenames). */
extern int dump_cnt;
/** @brief Simulation time at which the next dump should be written. */
extern double tdump;
/** @brief Simulation time at which the next log entry should be written. */
extern double tlog;

/** @brief Mass flux at the inner boundary (left face of first active X1 zone). */
extern double mdot;
/** @brief Energy flux at the inner boundary. */
extern double edot;
/** @brief Angular momentum flux at the inner boundary. */
extern double ldot;
/** @brief Mass flux at the nominal event horizon (i = NG + 5). */
extern double mdot_eh;
/** @brief Energy flux at the event horizon. */
extern double edot_eh;
/** @brief Angular momentum flux at the event horizon. */
extern double ldot_eh;
/** @brief Zone index (i) of the maximum divB location (diagnostic). */
extern int icurr;
/** @brief Zone index (j) of the maximum divB location (diagnostic). */
extern int jcurr;

/** @brief Number of active OpenMP threads (set once at startup). */
extern int nthreads;

#if ELECTRONS
/** @brief Maximum allowed total entropy K_tot (ceiling applied in fixup_ceiling()). */
#define KTOTMAX (3.)
/** @brief Electron mass in CGS units [g]. */
#define ME (9.1093826e-28  )
/** @brief Proton mass in CGS units [g]. */
#define MP (1.67262171e-24 )
/** @brief Electron adiabatic index (game, typically 4/3). */
extern double game;
/** @brief Proton (ion) adiabatic index (gamp, typically 5/3). */
extern double gamp;
/** @brief Initial fraction of internal energy in electrons (0 < fel0 < 1). */
extern double fel0;
/** @brief Minimum allowed ion-to-electron temperature ratio T_p/T_e. */
extern double tptemin;
/** @brief Maximum allowed ion-to-electron temperature ratio T_p/T_e. */
extern double tptemax;
#endif

/** @brief FMKS polynomial normalization factor for pole-derefining coordinate mapping. */
extern double poly_norm;
/** @brief FMKS polynomial transition parameter x_t. */
extern double poly_xt;
/** @brief FMKS polynomial steepness exponent alpha. */
extern double poly_alpha;
/** @brief FMKS smoothing scale for the transition between MKS and FMKS near the pole. */
extern double mks_smooth;

/** @} */ // end global_vars

/*---------------------------------------------------------------------------------
                                  MACROS
---------------------------------------------------------------------------------*/

/** @defgroup loop_macros Loop Macros
 *  @brief Convenience macros for iterating over grid zones, primitives, and spacetime indices.
 *
 *  Active zones are at indices [NG, N+NG-1]; ghost zones extend the range to [0, N+2*NG-1].
 *  The convention is that the inner loop is over i (X1 direction) and the outer loop over j (X2).
 *  @{
 */

/** @brief Loop over all active X1 zones: i in [NG, N1+NG). */
#define ILOOP	\
  for (int i = 0 + NG; i < N1 + NG; i++)
/** @brief Loop over all X1 zones including ghost zones: i in [0, N1+2*NG). */
#define ILOOPALL \
  for (int i = 0; i < N1 + 2*NG; i++)
/** @brief Loop over all active X2 zones: j in [NG, N2+NG). */
#define JLOOP	\
  for (int j = 0 + NG; j < N2 + NG; j++)
/** @brief Loop over all X2 zones including ghost zones: j in [0, N2+2*NG). */
#define JLOOPALL \
  for (int j = 0; j < N2 + 2*NG; j++)
/** @brief Nested loop over all active zones: JLOOP then ILOOP (j is outer). */
#define ZLOOP	\
  JLOOP ILOOP
/** @brief Nested loop over all zones including ghosts. */
#define ZLOOPALL \
  JLOOPALL ILOOPALL
/** @brief Transposed active zone loop: ILOOP then JLOOP (i is outer, for transposed output). */
#define ZLOOP_OUT \
  ILOOP JLOOP
/** @brief Transposed full-grid loop. */
#define ZLOOP_TRANSPOSE \
  ILOOPALL JLOOPALL

/** @brief Loop over a rectangular subsection of active zones in X1: i from istart to istop (inclusive). */
#define ISLOOP(istart, istop) \
  for (int i = (istart) + NG; i <= (istop) + NG; i++)
/** @brief Loop over a rectangular subsection of active zones in X2: j from jstart to jstop (inclusive). */
#define JSLOOP(jstart, jstop) \
  for (int j = (jstart) + NG; j <= (jstop) + NG; j++)
/** @brief Nested subsection loop: JSLOOP then ISLOOP. */
#define ZSLOOP(jstart, jstop, istart, istop) \
  for (int j = (jstart) + NG; j <= (jstop) + NG; j++) \
  for (int i = (istart) + NG; i <= (istop) + NG; i++)
/** @brief Reverse-order subsection loop (high to low indices). */
#define ZSLOOP_REVERSE(jstart, jstop, istart, istop) \
  for (int j = (jstop) + NG; j >= (jstart) + NG; j--) \
  for (int i = (istop) + NG; i >= (istart) + NG; i--)
/** @brief Transposed subsection loop: ISLOOP then JSLOOP. */
#define ZSLOOP_OUT(jstart, jstop, istart, istop) \
  ISLOOP(istart,istop) JSLOOP(jstart,jstop)

/** @brief Loop over all NVAR primitive variable indices. Declares int ip. */
#define PLOOP for(int ip = 0; ip < NVAR; ip++)

/** @brief Loop over a single spacetime index mu in [0, NDIM). */
#define DLOOP1 for (int mu = 0; mu < NDIM; mu++)
/** @brief Nested loop over two spacetime indices mu, nu in [0, NDIM). */
#define DLOOP2 for (int mu = 0; mu < NDIM; mu++)	\
               for (int nu = 0; nu < NDIM; nu++)

/** @} */ // end loop_macros

/** @brief Minimum of two values. */
#define MY_MIN(fval1,fval2) ( ((fval1) < (fval2)) ? (fval1) : (fval2))
/** @brief Maximum of two values. */
#define MY_MAX(fval1,fval2) ( ((fval1) > (fval2)) ? (fval1) : (fval2))
/** @brief Sign of a value (+1 or -1). */
#define MY_SIGN(fval) ( ((fval) <0.) ? -1. : 1. )
/** @brief Kronecker delta: 1 if i==j, 0 otherwise. */
#define delta(i,j) ((i == j) ? 1. : 0.)

/** @brief Log a message to stderr (only if DEBUG is set). */
#define LOG(msg) if(DEBUG) {fprintf(stderr, msg); fprintf(stderr,"\n");}
/** @brief Log a formatted message with one argument to stderr (only if DEBUG is set). */
#define LOGN(fmt,x) if(DEBUG) {fprintf(stderr, fmt, x); fprintf(stderr,"\n");}
/** @brief Print a fatal error message and exit. */
#define ERROR(msg) {fprintf(stderr, msg); fprintf(stderr,"\n"); exit(-1);}
/** @brief Conditional debug flag; can be used to print intermediate variable values in DEBUG mode. */
#define FLAG(msg) if(DEBUG) {LOG(msg);}

/*---------------------------------------------------------------------------------
                                FUNCTION DECLARATIONS
---------------------------------------------------------------------------------*/

/** @defgroup func_decls Function Declarations
 *  @brief Forward declarations of all functions, grouped by source file.
 *  @{
 */

// bl_coord.c
/** @brief Convert code coordinates X to Boyer-Lindquist r and theta. */
void bl_coord(const double X[NDIM], double *r, double *th);

// bounds.c
/** @brief Apply all boundary conditions to the full grid. */
void set_bounds(struct GridGeom *G, struct FluidState *S);
/** @brief Fix boundary fluxes (zero polar X2 fluxes, clip inflow). */
void fix_flux(struct FluidFlux *F);

// coord.c
/** @brief Return coordinate vector X at grid location (i, j, loc). */
void coord(int i, int j, int loc, double *X);
/** @brief Compute BL (r, th) from code coordinates X (same as bl_coord() in coord.c). */
void bl_coord(const double X[NDIM], double *r, double *th);
/** @brief Compute covariant metric g_mu_nu at a point X in code coordinates. */
void gcov_func(double *X, double gcov[NDIM][NDIM]);
/** @brief Compute the Jacobian dX^mu / dX^KS at a point X. */
void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM]);
/** @brief Set coordinate system parameters (startx, dx) from run-time parameters. */
void set_points();
/** @brief Initialize the full grid geometry struct G. */
void set_grid(struct GridGeom *G);
/** @brief Set grid geometry at a single location (i, j, loc). */
void set_grid_loc(struct GridGeom *G, int i, int j, int loc);
/** @brief Zero pflag and fail_save arrays. */
void zero_arrays();

// current.c
/** @brief Compute the 4-current j^mu from the Maxwell tensor using centered differences. */
void current_calc(struct GridGeom *G, struct FluidState *S, struct FluidState *Ssave, double dtsave);

// diag.c
/** @brief Compute mass/energy/angular momentum fluxes from the flux struct. */
void diag_flux(struct FluidFlux *F);
/** @brief Main diagnostic routine: computes and outputs global diagnostics, calls dump(). */
void diag(struct GridGeom *G, struct FluidState *S, int call_code);
/** @brief Compute the flux-CT divergence of B at zone (i, j). */
double flux_ct_divb(struct GridGeom *G, struct FluidState *S, int i, int j);

// electrons.c
#if ELECTRONS
/** @brief Initialize electron entropy variables from the initial fluid state. */
void init_electrons(struct GridGeom *G, struct FluidState *S);
/** @brief Heat electrons by distributing the dissipated entropy according to the chosen heating model. */
void heat_electrons(struct GridGeom *G, struct FluidState *Sh, struct FluidState *S);
/** @brief Apply floors/ceilings to electron entropy variables. */
void fixup_electrons(struct FluidState *S);
#endif

// fixup.c
/** @brief Apply floors and ceilings to all zones in the domain. */
void fixup(struct GridGeom *G, struct FluidState *S);
/** @brief Replace zones where U_to_P() failed with values interpolated from neighbors. */
void fixup_utoprim(struct GridGeom *G, struct FluidState *S);

// fluxes.c
/** @brief Reconstruct primitives at faces and compute numerical fluxes via LLF Riemann solver. */
double get_flux(struct GridGeom *G, struct FluidState *S, struct FluidFlux *F);
/** @brief Perform constrained transport (Toth 2000) to enforce div(B) = 0. */
void flux_ct(struct FluidFlux *F);

// io.c
/** @brief Write a regular numbered dump file. */
void dump(struct GridGeom *G, struct FluidState *S);
/** @brief Backend dump writer (handles both regular and abort dumps). */
void dump_backend(struct GridGeom *G, struct FluidState *S, int type);
/** @brief Write the grid coordinate and metric file (written once at dump_cnt == 0). */
void dump_grid(struct GridGeom *G);

// metric.c
/** @brief Compute contravariant metric g^mu^nu from covariant metric g_mu_nu; return sqrt(|det g|). */
double gcon_func(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM]);
/** @brief Copy covariant metric from GridGeom struct to a local 4x4 array. */
void get_gcov(struct GridGeom *G, int i, int j, int loc, double gcov[NDIM][NDIM]);
/** @brief Copy contravariant metric from GridGeom struct to a local 4x4 array. */
void get_gcon(struct GridGeom *G, int i, int j, int loc, double gcon[NDIM][NDIM]);
/** @brief Compute Christoffel connection coefficients Gamma^lam_mu_nu at zone (i, j) via finite differences. */
void conn_func(struct GridGeom *G, int i, int j);
/** @brief Lower a contravariant vector v^mu to covariant v_mu at zone (i, j, loc). */
void lower_grid(GridVector vcon, GridVector vcov, struct GridGeom *G, int i, int j, int loc);
/** @brief Raise a covariant vector v_mu to contravariant v^mu at zone (i, j, loc). */
void raise_grid(GridVector vcov, GridVector vcon, struct GridGeom *G, int i, int j, int loc);
/** @brief Compute the inner product v^mu * v_mu. */
double dot(double vcon[NDIM], double vcov[NDIM]);
/** @brief Compute the inverse of a 4x4 matrix; return the determinant. */
double invert(double *m, double *inv);

void pack_write_axiscalar(double in[N2+2*NG][N1+2*NG], const char* name, size_t datatype);
void pack_write_Gtensor(double in[NDIM][NDIM][N2+2*NG][N1+2*NG], const char* name, size_t datatype);

// params.c
/** @brief Register all core parameters in the parameter table. */
void set_core_params();
/** @brief Register a single parameter key-pointer pair. */
void set_param(char *key, void *data);
/** @brief Read parameter values from an ASCII param.dat file. */
void read_params(char *pfname);

// phys.c
/** @brief Compute primitive-to-conserved/flux transform at a single zone (i, j). */
void prim_to_flux(struct GridGeom *G, struct FluidState *S, int i, int j, int dir, int loc, GridPrim flux);
/** @brief Vectorized (OpenMP) version of prim_to_flux() over a rectangular zone range. */
void prim_to_flux_vec(struct GridGeom *G, struct FluidState *S, int dir, int loc, int jstart, int jstop, int istart, int istop, GridPrim flux);
/** @brief Compute magnetic 4-vector b^mu from primitive B fields and 4-velocity. */
void bcon_calc(struct FluidState *S, int i, int j);
/** @brief Compute the MHD stress-energy tensor T^dir_mu at zone (i, j). */
void mhd_calc(struct FluidState *S, int i, int j, int dir, double *mhd);
/** @brief Compute the geometric source terms -T^mu_nu * Gamma^nu_alpha_mu. */
void get_fluid_source(struct GridGeom *G, struct FluidState *S, GridPrim *dU);
/** @brief Return b^mu * b_mu (twice the magnetic pressure). */
double bsq_calc(struct FluidState *S, int i, int j);
/** @brief Compute all auxiliary 4-vectors (ucon, ucov, bcon, bcov) at a single zone. */
void get_state(struct GridGeom *G, struct FluidState *S, int i, int j, int loc);
/** @brief Vectorized get_state() over a rectangular zone range. */
void get_state_vec(struct GridGeom *G, struct FluidState *S, int loc, int jstart, int jstop, int istart, int istop);
/** @brief Compute the contravariant 4-velocity u^mu from primitives at (i, j, loc). */
void ucon_calc(struct GridGeom *G, struct FluidState *S, int i, int j, int loc);
/** @brief Compute the Lorentz factor gamma w.r.t. the normal observer. */
double mhd_gamma_calc(struct GridGeom *G, struct FluidState *S, int i, int j, int loc);
/** @brief Compute the maximum and minimum fast magnetosonic wave speeds at (i, j, loc) in direction dir. */
void mhd_vchar(struct GridGeom *G, struct FluidState *Sr, int i, int j, int loc, int dir, GridDouble cmax, GridDouble cmin);

// problem.c
/** @brief Register problem-specific parameters in the parameter table. */
void set_problem_params();
/** @brief Initialize the fluid state for the chosen problem. */
void init(struct GridGeom *G, struct FluidState *S);
/** @brief Problem-specific outer (X1R) boundary condition (e.g., Bondi analytic values). */
void bound_gas_prob_x1r(int i, int j, GridPrim  P, struct GridGeom *G);
/** @brief Write problem-specific header data to a dump file. */
void save_problem_data(FILE *fp);

// random.c
/** @brief Initialize the Mersenne Twister PRNG with integer seed. */
void init_rng(int seed);
/** @brief Return a pseudo-random double in [0, 1). */
double ran_uniform();

// reconstruction.c
/** @brief Reconstruct left (Pl) and right (Pr) primitive states at cell faces in direction dir. */
void reconstruct(struct FluidState *S, GridPrim Pl, GridPrim Pr, int dir);

// restart.c
/** @brief Write a numbered restart file. */
void restart_write(struct FluidState *S);
/** @brief Backend restart writer (handles regular and abort restarts). */
void restart_write_backend(struct FluidState *S, int type);
/** @brief Read a restart file from disk into a FluidState. */
void restart_read(char *fname, struct FluidState *S);
/** @brief Check for an existing restart file; if found, initialize the simulation from it. Returns 1 if restarted. */
int restart_init(struct GridGeom *G, struct FluidState *S);

// step.c
/** @brief Advance the simulation by one full timestep using a 2nd-order predictor-corrector scheme. */
void step(struct GridGeom *G, struct FluidState *S);

// timing.c
/** @brief Initialize all performance timers to zero. */
void time_init();
/** @brief Start a named performance timer. */
void timer_start(int timerCode);
/** @brief Stop a named performance timer and accumulate elapsed time. */
void timer_stop(int timerCode);
/** @brief Print per-step timer breakdown and zone-cycles-per-second metric. */
void report_performance();

// u_to_p.c
/** @brief Invert conservative variables U to primitive variables P at zone (i, j, loc). Returns 0 on success, nonzero on failure. */
int U_to_P(struct GridGeom *G, struct FluidState *S, int i, int j, int loc);

/** @} */ // end func_decls
