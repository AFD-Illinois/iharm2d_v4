/*---------------------------------------------------------------------------------

  IO.C

  -WRITE DUMP FILE
  -WRITE GRID FILE

---------------------------------------------------------------------------------*/

#include "decs.h"

#include <sys/stat.h>
#include <ctype.h>

#define FML_DBL_OUT "%28.18e"
#define FML_INT_OUT "%10d"
#define STRING_OUT "%15s"

// Write REGULAR dump file by calling dump_write_backend with type=IO_REGULAR
void dump(struct GridGeom *G, struct FluidState *S)
{
  dump_backend(G, S, IO_REGULAR);
}

// The function that actually writes the dump file
void dump_backend(struct GridGeom *G, struct FluidState *S, int type)
{
  // Start timer for dump write
  timer_start(TIMER_IO);

  static GridDouble *data;
  char fname[80];
  
  // Set prim names
  // Not writing it out though
  #if ELECTRONS
  #if ALLMODELS
  const char varNames[NVAR][STRLEN] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3", "KTOT", "KEL0", "KEL1", "KEL2", "KEL3"};
  #else
  const char varNames[NVAR][STRLEN] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3", "KTOT", "KEL0"};
  #endif
  #else
  const char varNames[NVAR][STRLEN] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3"};
  #endif

  static int firstc = 1;
  if (firstc)
  {
    data = calloc(1, sizeof(GridDouble));
    firstc = 0;
  }

  // Dont re-dump the grid after a restart
  if (dump_cnt == 0)
    dump_grid(G);

  if (type == IO_REGULAR)
    sprintf(fname, "dumps/dump_%08d", dump_cnt);
  else if(type == IO_ABORT)
    sprintf(fname, "dumps/dump_abort");
  
  fprintf(stdout, "DUMP %s\n", fname);

  FILE *fp;
  fp = fopen(fname, "wt");

  // Write dump header
  save_problem_data(fp); 
 
  fprintf(fp, " %s", VERSION);
  int has_electrons = ELECTRONS;
  fprintf(fp, FML_INT_OUT , has_electrons);
  char gridfile[15] = "grid";
  fprintf(fp, STRING_OUT, gridfile);

  #if METRIC == MINKOWSKI
  fprintf(fp, STRING_OUT, "MINKOWSKI");
  #elif METRIC == MKS
  #if DEREFINE_POLES
  fprintf(fp, STRING_OUT, "FMKS");
  #else
  fprintf(fp, STRING_OUT, "MKS");
  #endif
  #endif
  
  #if RECONSTRUCTION == LINEAR
  fprintf(fp, STRING_OUT, "LINEAR");
  #elif RECONSTRUCTION == PPM
  fprintf(fp, STRING_OUT, "PPM");
  #elif RECONSTRUCTION == WENO
  fprintf(fp, STRING_OUT, "WENO");
  #elif RECONSTRUCTION == MP%
  fprintf(fp, STRING_OUT, "MP5");
  #endif

  fprintf(fp, FML_INT_OUT, N1);
  fprintf(fp, FML_INT_OUT, N2);

  int n_prims = NVAR;
  fprintf(fp, FML_INT_OUT, n_prims);
  int n_prims_passive = 0;
  fprintf(fp, FML_INT_OUT, n_prims_passive);
  //for (i=0; i<NVAR; i++)
  //  fprintf(fp, STRING_OUT, varNames[i]);

  #if ELECTRONS
  fprintf(fp, FML_DBL_OUT, game);
  fprintf(fp, FML_DBL_OUT, gamp);
  fprintf(fp, FML_DBL_OUT, fel0);
  fprintf(fp, FML_DBL_OUT, tptemin);
  fprintf(fp, FML_DBL_OUT, tptemax);
  #endif
  
  fprintf(fp, FML_DBL_OUT, gam);
  fprintf(fp, FML_DBL_OUT, cour);
  fprintf(fp, FML_DBL_OUT, tf);
   
  fprintf(fp, FML_DBL_OUT, startx[1]);
  fprintf(fp, FML_DBL_OUT, startx[2]);
  fprintf(fp, FML_DBL_OUT, dx[1]);
  fprintf(fp, FML_DBL_OUT, dx[2]);
  int n_dim = NDIM;
  fprintf(fp, FML_INT_OUT, n_dim);

  #if METRIC == MKS
  #if DEREFINE_POLES
  fprintf(fp, FML_DBL_OUT, poly_xt);
  fprintf(fp, FML_DBL_OUT, poly_alpha);
  fprintf(fp, FML_DBL_OUT, mks_smooth);
  #endif
  fprintf(fp, FML_DBL_OUT, Rin);
  fprintf(fp, FML_DBL_OUT, Rout);
  fprintf(fp, FML_DBL_OUT, Rhor);
  double z1 = 1 + pow(1 - a*a, 1./3.)*(pow(1 + a, 1./3.) + pow(1 - a, 1./3.));
  double z2 = sqrt(3*a*a + z1*z1);
  double Risco = 3 + z2 - sqrt((3 - z1)*(3 + z1 + 2*z2));
  fprintf(fp, FML_DBL_OUT, Risco);
  fprintf(fp, FML_DBL_OUT, hslope);
  fprintf(fp, FML_DBL_OUT, a);
  #endif

  fprintf(fp, FML_DBL_OUT, t);
  fprintf(fp, FML_DBL_OUT, dt);
  fprintf(fp, FML_INT_OUT, nstep);
  fprintf(fp, FML_INT_OUT, dump_cnt);
  fprintf(fp, FML_DBL_OUT, DTd);
  fprintf(fp, FML_DBL_OUT, DTf);

  fprintf(fp, "\n");

  // Note that the arrays must be reversed along (X1,X2)
  // Write prims
  ZLOOP_OUT
  {
    PLOOP
      fprintf(fp, FML_DBL_OUT, S->P[ip][j][i]);

  // Write current
    DLOOP1
      fprintf(fp, FML_DBL_OUT, S->jcon[mu][j][i]);
  
  // Write Lorentz factor
    fprintf(fp, FML_DBL_OUT, mhd_gamma_calc(G, S, i, j, CENT));
  
  // Write divB
    if (i > 0+NG && j > 0+NG)
      fprintf(fp, FML_DBL_OUT, flux_ct_divb(G, S, i, j));
    else
      fprintf(fp, FML_DBL_OUT, 0.);


  // Write u_to_p fails
    fprintf(fp, FML_INT_OUT, fail_save[j][i]);
    
  // Write fixup fails
    fprintf(fp, FML_INT_OUT, fflag[j][i]);

    fprintf(fp, "\n");
  }

  // Reset u_to_p failure flags to zero
  ZLOOP
    fail_save[j][i] = 0;

  // Close file pointer
  fclose(fp);

  // Stop timer for dump write
  timer_stop(TIMER_IO);
}

#define NGRIDVARS 6
// The function that actually writes the grid file
void dump_grid(struct GridGeom *G)
{
  GridDouble *x[NGRIDVARS];
  for (int d = 0; d < NGRIDVARS; d++) x[d] = calloc(1, sizeof(GridDouble));
  const char *coordNames[] = {"X", "Y", "Z", "r", "th", "phi", "X1", "X2", "X3"};

  char *fname = "dumps/grid";
  fprintf(stdout, "GRID %s\n", fname);

  FILE *fp;
  fp = fopen(fname, "wt");

  ZLOOP
  {
    double xp[NDIM];
    coord(i, j, CENT, xp);
    #if METRIC == MINKOWSKI
    (*x[0])[j][i] = xp[1];
    (*x[1])[j][i] = xp[2];
    (*x[2])[j][i] = 0;
    (*x[3])[j][i] = 0;
    #elif METRIC == MKS
    double r, th;
    bl_coord(xp, &r, &th);
    (*x[0])[j][i] = r*sin(th);
    (*x[1])[j][i] = r*cos(th);
    (*x[2])[j][i] = r;
    (*x[3])[j][i] = th;
    #endif
    (*x[4])[j][i] = xp[1];
    (*x[5])[j][i] = xp[2];
  }
  
  ZLOOP_OUT
  {
    for (int d = 0; d < NGRIDVARS; d++)
    {
      fprintf(fp, FML_DBL_OUT, (*x[d])[j][i]);
    }
    fprintf(fp, FML_DBL_OUT, G->gdet[CENT][j][i]);    
    fprintf(fp, FML_DBL_OUT, G->lapse[CENT][j][i]);
    DLOOP2
      fprintf(fp, FML_DBL_OUT, G->gcon[CENT][mu][nu][j][i]);
    DLOOP2
      fprintf(fp, FML_DBL_OUT, G->gcov[CENT][mu][nu][j][i]);

    fprintf(fp, "\n");
  }
  // Close file pointer
  fclose(fp);

  // Free memory
  for (int d = 0; d < NGRIDVARS; d++) free(x[d]);
}
