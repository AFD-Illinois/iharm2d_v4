/*---------------------------------------------------------------------------------

  DIAG.C

  -Diagnostic output. Also used for writing dumps
  -Computes mass, magnetic, angular momentum and energy flux at first physical
   zone and at EH
  -Computes total mass, angular momentum and energy
  -Computes divB
  -Computes pseudo-emissivity for thermal synchrotron radiation and luminosity
  -Write out diagnostics to dumps/log.out

---------------------------------------------------------------------------------*/

#include "decs.h"

// Evaluate flux based diagnostics; put results in global variables
// Note this is still per-process
void diag_flux(struct FluidFlux *F)
{
  mdot = edot = ldot = 0.;
  mdot_eh = edot_eh = ldot_eh = 0.;
  int iEH = NG + 5;
#if !INTEL_WORKAROUND
#pragma omp parallel for \
  reduction(+:mdot) reduction(+:edot) reduction(+:ldot) \
  reduction(+:mdot_eh) reduction(+:edot_eh) reduction(+:ldot_eh)
#endif
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
#if !INTEL_WORKAROUND
#pragma omp parallel for \
  reduction(+:rmed) reduction(+:pp) reduction(+:e) \
  reduction(max:divbmax) collapse(2)
#endif
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
#if !INTEL_WORKAROUND
#pragma omp parallel for reduction(+:mass) reduction(+:egas) reduction(+:Phi) reduction(+:jet_EM_flux) reduction(+:lum_eht) collapse(2)
#endif
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

// Calculate divB
double flux_ct_divb(struct GridGeom *G, struct FluidState *S, int i, int j)
{
  return fabs(0.5*(S->P[B1][j][i]*G->gdet[CENT][j][i] + S->P[B1][j-1][i]*G->gdet[CENT][j-1][i] 
  - S->P[B1][j][i-1]*G->gdet[CENT][j][i-1] - S->P[B1][j-1][i-1]*G->gdet[CENT][j-1][i-1])/dx[1] 
  + 0.5*(S->P[B2][j][i]*G->gdet[CENT][j][i] + S->P[B2][j][i-1]*G->gdet[CENT][j][i-1] 
  - S->P[B2][j-1][i]*G->gdet[CENT][j-1][i] - S->P[B2][j-1][i-1]*G->gdet[CENT][j-1][i-1])/dx[2]);
}
