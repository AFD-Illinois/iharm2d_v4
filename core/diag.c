/**
 * @file diag.c
 * @brief Simulation diagnostics: accretion rates, energy integrals, div-B, and log output.
 *
 * @details Provides two public diagnostic functions and one helper:
 *
 * - **diag_flux()**: Accumulates the mass accretion rate (mdot), energy flux
 *   (edot), and angular-momentum flux (ldot) through the innermost active
 *   radial face (i=NG) and through i=NG+5 (event horizon proxy) by summing
 *   the X1 fluxes over all j zones.
 *
 * - **diag()**: Called at initialisation (DIAG_INIT), each log step
 *   (DIAG_LOG), dump steps (DIAG_DUMP), abort (DIAG_ABORT), and end-of-run
 *   (DIAG_FINAL).  Computes grid-integrated conserved mass, angular momentum,
 *   internal energy, mass-weighted radius, magnetic flux Φ at r=5, and the
 *   EHT (2019) synchrotron proxy luminosity.  Writes a row of diagnostics
 *   to `dumps/log.out` and (at appropriate call codes) triggers dump() or
 *   dump_backend().
 *
 * - **flux_ct_divb()**: Returns the cell-centred estimate of |∇·B| using the
 *   standard Constrained-Transport stencil (four-zone difference).
 *
 * @note All reduction loops are OpenMP-parallelised with appropriate reduction
 * clauses.
 */

#include "decs.h"

/**
 * @brief Compute flux-based accretion diagnostics at the inner boundary and EH proxy.
 *
 * @details Integrates the X1 fluxes of RHO, UU (minus RHO), and U3 over all
 * j zones at:
 * - The innermost active face (i = NG): sets mdot, edot, ldot.
 * - The proxy event-horizon face (i = NG+5): sets mdot_eh, edot_eh, ldot_eh.
 * Each flux is multiplied by dx[2] to convert to a line integral.
 * The results are stored in the corresponding global variables.
 *
 * @param F  Fluid flux struct; X1 face fluxes are read (not modified).
 */
// Evaluate flux based diagnostics; put results in global variables
void diag_flux(struct FluidFlux *F)
{
  mdot = edot = ldot = 0.;
  mdot_eh = edot_eh = ldot_eh = 0.;
  int iEH = NG + 5;
#pragma omp parallel for \
  reduction(+:mdot) reduction(+:edot) reduction(+:ldot) \
  reduction(+:mdot_eh) reduction(+:edot_eh) reduction(+:ldot_eh)
    JSLOOP(0, N2 - 1)
    {
      mdot += -F->X1[RHO][j][NG]*dx[2];
      edot += (F->X1[UU][j][NG] - F->X1[RHO][j][NG])*dx[2];
      ldot += F->X1[U3][j][NG]*dx[2];
      mdot_eh += -F->X1[RHO][j][iEH]*dx[2];
      edot_eh += (F->X1[UU][j][iEH] - F->X1[RHO][j][iEH])*dx[2];
      ldot_eh += F->X1[U3][j][iEH]*dx[2];
    }
}

/**
 * @brief Main diagnostic routine: compute integrals, write log, trigger dumps.
 *
 * @details Behaviour depends on call_code:
 * - **DIAG_INIT**: Opens `dumps/log.out`, performs initial integrals, writes
 *   first log row, and calls dump() (unless restarting).
 * - **DIAG_LOG**: Recomputes integrals and writes a row to the log file.
 * - **DIAG_DUMP**: Calls dump() and increments dump_cnt.
 * - **DIAG_ABORT**: Calls dump_backend() with IO_ABORT.
 * - **DIAG_FINAL**: Computes integrals, writes final log row, and calls dump().
 *
 * Computed quantities (for LOG and FINAL call codes):
 *   - Grid-integrated conserved RHO (rmed), U3 (pp), UU (e).
 *   - Max |div B| across all active zones (via flux_ct_divb()).
 *   - Total mass (mass) and gas energy (egas).
 *   - Magnetic flux Φ = 0.5 * ∫|B1| gdet dx2 at i = NG+5 (EH proxy).
 *   - EHT 2019 synchrotron luminosity proxy j_eht.
 *
 * @param G          Grid geometry.
 * @param S          Fluid state.
 * @param call_code  One of DIAG_INIT, DIAG_LOG, DIAG_DUMP, DIAG_ABORT, DIAG_FINAL.
 */
// Additional diagnostics: divB, conserved quantities: rho*u^t, T^_t, T^t_k,
// jet luminosity (see EHT CC'19), total integrated conserved mass and conserved energy
// Calls dump_backend to write out dump
// Prints diagnostics to log file
void diag(struct GridGeom *G, struct FluidState *S, int call_code)
{
  static FILE *ener_file;

  if (call_code == DIAG_INIT)
  {
    // Set things up
      ener_file = fopen("dumps/log.out", "a");
      if (ener_file == NULL)
      {
        fprintf(stderr, "Error opening log file!\n");
        exit(1);
      }
  }

  double pp = 0.;
  double divbmax = 0.; 
  int imax = 0; int jmax = 0;
  double rmed = 0.;
  double e = 0.;
  // Calculate conserved quantities
  if (call_code == DIAG_INIT || call_code == DIAG_LOG || call_code == DIAG_FINAL)
  {
    get_state_vec(G, S, CENT, 0, N2 - 1, 0, N1 - 1);
    prim_to_flux_vec(G, S, 0, CENT, 0, N2 - 1, 0, N1 - 1, S->U);
#pragma omp parallel for \
  reduction(+:rmed) reduction(+:pp) reduction(+:e) \
  reduction(max:divbmax) collapse(2)
    ZLOOP
    {
      rmed += S->U[RHO][j][i]*dV;
      pp += S->U[U3][j][i]*dV;
      e += S->U[UU][j][i]*dV;
      if (i > 0+NG && j > 0+NG)
      {
        double divb = flux_ct_divb(G, S, i, j);
        if (divb > divbmax)
        {
          divbmax = divb;
          imax = i;
          jmax = j;
        }
      }
    }
  }

  double mass = 0.;
  double egas = 0.;
  double Phi = 0.;
  double jet_EM_flux = 0.;
  double lum_eht = 0.;
#pragma omp parallel for reduction(+:mass) reduction(+:egas) reduction(+:Phi) reduction(+:jet_EM_flux) reduction(+:lum_eht) collapse(2)
  ZLOOP
  {
    mass += S->U[RHO][j][i]*dV;
    egas += S->U[UU][j][i]*dV;
    double rho = S->P[RHO][j][i];
    double Pg = (gam - 1.)*S->P[UU][j][i];
    double bsq = bsq_calc(S, i, j);
    double Bmag = sqrt(bsq);
    double C_eht = 0.2;
    double j_eht = pow(rho,3.)*pow(Pg,-2.)*exp(-C_eht*pow(rho*rho/(Bmag*Pg*Pg),1./3.));
    lum_eht += j_eht*dV*G->gdet[CENT][j][i];
    if (i == 5+NG)
      Phi += 0.5*fabs(sqrt(4*M_PI)*S->P[B1][j][i])*dx[2]*G->gdet[CENT][j][i];
  }
  
  if ((call_code == DIAG_INIT && !is_restart) || call_code == DIAG_DUMP || call_code == DIAG_FINAL)
  {
    dump(G, S);
    dump_cnt++;
  }

  if (call_code == DIAG_ABORT)
  {
    dump_backend(G, S, IO_ABORT);
  }

  if (call_code == DIAG_INIT || call_code == DIAG_LOG || call_code == DIAG_FINAL)
  {
    //mdot will be negative w/scheme above
    double phi = Phi/sqrt(fabs(mdot) + SMALL);
    fprintf(stdout, "LOG      t=%g \t divbmax: %g at %d %d\n", t, divbmax, imax, jmax);
    fprintf(ener_file, "%10.5g %10.5g %10.5g %10.5g %15.8g %15.8g ", t, rmed, pp, e, S->P[UU][N2/2][N1/2]*pow(S->P[RHO][N2/2][N1/2], -gam), S->P[UU][N2/2][N1/2]);
    fprintf(ener_file, "%15.8g %15.8g %15.8g ", mdot, edot, ldot);
    fprintf(ener_file, "%15.8g %15.8g ", mass, egas);
    fprintf(ener_file, "%15.8g %15.8g %15.8g ", Phi, phi, jet_EM_flux);
    fprintf(ener_file, "%15.8g ", divbmax);
    fprintf(ener_file, "%15.8g ", lum_eht);
    fprintf(ener_file, "%15.8g %15.8g %15.8g ", mdot_eh, edot_eh, ldot_eh);
    fprintf(ener_file, "\n");
    fflush(ener_file);
  }
}

/**
 * @brief Compute the cell-centred constrained-transport estimate of |div B| at zone (i, j).
 *
 * @details Uses the standard four-zone stencil for the CT magnetic flux divergence:
 * @f[
 *   (\nabla \cdot B)_{ij} = \frac{1}{2\Delta X^1}
 *     \left[ B^1_{ij}\,g_{ij} + B^1_{i,j-1}\,g_{i,j-1}
 *          - B^1_{i-1,j}\,g_{i-1,j} - B^1_{i-1,j-1}\,g_{i-1,j-1} \right]
 *   + \frac{1}{2\Delta X^2}
 *     \left[ B^2_{ij}\,g_{ij} + B^2_{i-1,j}\,g_{i-1,j}
 *          - B^2_{i,j-1}\,g_{i,j-1} - B^2_{i-1,j-1}\,g_{i-1,j-1} \right]
 * @f]
 * where @f$g = \sqrt{-g}@f$ (gdet).  In Constrained Transport, this measure should
 * remain at machine precision if the induction equation is evolved exactly.
 *
 * @param G  Grid geometry (gdet at CENT is read at four neighbouring zones).
 * @param S  Fluid state (P[B1] and P[B2] are read at the same four zones).
 * @param i  X1 zone index (must satisfy i > NG).
 * @param j  X2 zone index (must satisfy j > NG).
 * @return   |∇·B| at zone (i, j).
 */
// Calculate divB
double flux_ct_divb(struct GridGeom *G, struct FluidState *S, int i, int j)
{
  return fabs(0.5*(S->P[B1][j][i]*G->gdet[CENT][j][i] + S->P[B1][j-1][i]*G->gdet[CENT][j-1][i] 
  - S->P[B1][j][i-1]*G->gdet[CENT][j][i-1] - S->P[B1][j-1][i-1]*G->gdet[CENT][j-1][i-1])/dx[1] 
  + 0.5*(S->P[B2][j][i]*G->gdet[CENT][j][i] + S->P[B2][j][i-1]*G->gdet[CENT][j][i-1] 
  - S->P[B2][j-1][i]*G->gdet[CENT][j-1][i] - S->P[B2][j-1][i-1]*G->gdet[CENT][j-1][i-1])/dx[2]);
}
