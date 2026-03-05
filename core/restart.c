/**
 * @file restart.c
 * @brief ASCII restart file writing and reading.
 *
 * @details Provides four functions for checkpointing:
 *
 * - **restart_write()**: Writes an indexed restart file
 *   `restarts/restart_NNNNNNNN` and updates the symlink `restarts/restart.last`
 *   to point to the newest file.  Calls restart_write_backend() internally.
 *
 * - **restart_write_backend()**: Core writer.  Outputs a compact ASCII header
 *   (VERSION, grid dimensions, simulation time, runtime parameters, metric
 *   parameters) followed by all NVAR primitives for every zone including
 *   ghost zones (ZLOOPALL).  Can be called with type=IO_ABORT to write
 *   `restarts/restart_abort` on a fatal error.
 *
 * - **restart_read()**: Reads a restart file by file name, parses the header,
 *   and fills in all global simulation parameters and S->P[][j][i].
 *
 * - **restart_init()**: High-level restart entry point called from main().
 *   Looks for `restarts/restart.last`; if found, calls restart_read() then
 *   re-initializes the grid, 4-vectors, conserved variables, and boundary
 *   conditions so the simulation is ready to continue.
 *   Returns 1 on success, 0 if no restart file is present.
 *
 * @note For MKS runs, tf (final time) is intentionally not restored from the
 * restart file; it is kept from the current parameter file so the user can
 * extend a run without editing the restart file.
 */

/*---------------------------------------------------------------------------------

  RESTART.C

  -Write restart file
  -Read restart file
  -Restart a simulation by reading an existing restart file

---------------------------------------------------------------------------------*/

#include "decs.h"

#include <ctype.h>

// Set output precision
#define FML_DBL_OUT "%28.18e"
#define FML_INT_OUT "%10d"
#define STRING_OUT "%25s"

static int restart_id = 0;

// Declare known sizes for outputting primitives
//static size_t fdims[] = {NVAR, N3, N2, N1};
//static size_t mdims = {NVAR, N3+2*NG, N2+2*NG, N1+2*NG};
//static size_t mstart = {0,NG, NG, NG};

/**
 * @brief Write a regular restart checkpoint file.
 *
 * @details Delegates to restart_write_backend() with type=IO_REGULAR, then
 * updates `restarts/restart.last` to symlink to the newly written file.
 *
 * @param S  Fluid state; P[][j][i] is serialized to disk.
 */
// Write REGULAR restart file by calling restart_write_backend with type=IO_REGULAR
void restart_write(struct FluidState *S)
{
  restart_write_backend(S, IO_REGULAR);
}

/**
 * @brief Core restart writer: serializes the full simulation state to an ASCII file.
 *
 * @details Writes:
 *  1. A header line: VERSION, N1, N2, simulation time t, step counter nstep,
 *     tf, gam, (game, gamp, fel0 if ELECTRONS), cour, DTd/DTf/DTl/DTr/DTp,
 *     restart_id, dump_cnt, dt, and metric parameters.
 *  2. Zone data: all NVAR primitives for every zone (including ghost zones),
 *     one zone per line (ZLOOPALL).
 * On completion, `restarts/restart.last` is re-linked to the new file.
 *
 * @param S     Fluid state.
 * @param type  IO_REGULAR → `restarts/restart_NNNNNNNN`; IO_ABORT → `restarts/restart_abort`.
 */
// The function that actually writes the restart file
void restart_write_backend(struct FluidState *S, int type)
{
  // Start timer for restart write
  timer_start(TIMER_RESTART);

  // Keep track of our own index
  restart_id++;
  
  char fname[STRLEN];

  // Set restart file name
  if (type == IO_REGULAR)
    sprintf(fname, "restarts/restart_%08d", restart_id);
  else if (type == IO_ABORT)
    sprintf(fname, "restarts/restart_abort");

  // Open file
  FILE *fp;
  fp = fopen(fname, "wt");

  // Write restart header data
  fprintf(fp, STRING_OUT, VERSION);
  fprintf(fp, FML_INT_OUT, N1);
  fprintf(fp, FML_INT_OUT, N2);
  fprintf(fp, FML_DBL_OUT, t);
  fprintf(fp, FML_INT_OUT, nstep);
  fprintf(fp, FML_DBL_OUT, tf);
  fprintf(fp, FML_DBL_OUT, gam);
  #if ELECTRONS
  fprintf(fp, FML_DBL_OUT, game);
  fprintf(fp, FML_DBL_OUT, gamp);
  fprintf(fp, FML_DBL_OUT, fel0);
  #endif
  fprintf(fp, FML_DBL_OUT, cour);
  fprintf(fp, FML_DBL_OUT, DTd);
  fprintf(fp, FML_DBL_OUT, DTf);
  fprintf(fp, FML_DBL_OUT, DTl);
  fprintf(fp, FML_INT_OUT, DTr);
  fprintf(fp, FML_INT_OUT, DTp);
  fprintf(fp, FML_INT_OUT, restart_id);
  fprintf(fp, FML_INT_OUT, dump_cnt);
  fprintf(fp, FML_DBL_OUT, dt);
  #if METRIC==MKS
  fprintf(fp, FML_DBL_OUT, Rin);
  fprintf(fp, FML_DBL_OUT, Rout);
  fprintf(fp, FML_DBL_OUT, a);
  fprintf(fp, FML_DBL_OUT, hslope);
  fprintf(fp, FML_DBL_OUT, Rhor);
  #else
  fprintf(fp, FML_DBL_OUT, x1Min);
  fprintf(fp, FML_DBL_OUT, x1Max);
  fprintf(fp, FML_DBL_OUT, x2Min);
  fprintf(fp, FML_DBL_OUT, x2Max);
  #endif
   
  fprintf(fp, "\n");

  // Write prims
  ZLOOPALL
  {
    PLOOP
      fprintf(fp, FML_DBL_OUT, S->P[ip][j][i]);
    fprintf(fp, "\n");
  }

  // Close file pointer
  fclose(fp);

  fprintf(stdout, "RESTART %s\n", fname);

  // Symlink when we're done writing (link to last good file)
  char fname_nofolder[80];
  sprintf(fname_nofolder, "restart_%08d", restart_id);

  // Chained OS functions: switch to restart directory,
  // remove current link, link last file, switch back
  int errcode;
  errcode = chdir("restarts");
  if (access("restart.last", F_OK) != -1)
    errcode = errcode || remove("restart.last");
  errcode = errcode || symlink(fname_nofolder, "restart.last");
  errcode = errcode || chdir("..");
  if (errcode != 0)
  {
    printf("Symlink failed: errno %d\n", errno);
    exit(-1);
  }

  // Stop timer for restart write
  timer_stop(TIMER_RESTART);
}

/**
 * @brief Read a restart file and restore all global simulation state.
 *
 * @details Parses the ASCII header written by restart_write_backend(), restoring:
 * simulation time t, step counter nstep, gam, cour, DTd/DTf/DTl/DTr/DTp,
 * restart_id, dump_cnt, dt, and metric parameters (a, hslope, Rin, Rout,
 * Rhor for MKS; x-domain extents for MINKOWSKI).
 * Then reads all NVAR primitives into S->P[][j][i] for every zone (ZLOOPALL).
 * Aborts if N1 or N2 in the file differs from the compiled values.
 *
 * @param fname  Path to the restart file.
 * @param S      Fluid state; P[][] is overwritten.
 */
// Read restart file pointed to by 'fname'
void restart_read(char *fname, struct FluidState *S)
{
  // Set file pointer
  FILE *fp;
  fp = fopen(fname, "r");
  
  int buffer_int;
  double buffer_double;
  
  // Read header info
  char version[20];
  fscanf(fp, "%s", version);
  fprintf(stderr, "Restarting from %s, file version %s\n\n", fname, version);
  
  int n1, n2;
  fscanf(fp, "%d", &n1);
  fscanf(fp, "%d", &n2);
  if (n1 != N1 || n2 != N2)
  {
    fprintf(stderr, "Restart file is wrong size!\n");
    exit(-1);
  }
  fscanf(fp, "%lf", &t);  
  fscanf(fp, "%d", &nstep);  
  if (METRIC != MKS)
    fscanf(fp, "%lf", &tf);
  else
    fscanf(fp, "%lf", &buffer_double);  
  fscanf(fp, "%lf", &gam);
  #if ELECTRONS
  fscanf(fp, "%lf", &game);
  fscanf(fp, "%lf", &gamp);
  fscanf(fp, "%lf", &fel0);
  #endif
  if (METRIC != MKS)
  {
    fscanf(fp, "%lf", &cour);
    fscanf(fp, "%lf", &DTd);
    fscanf(fp, "%lf", &DTf);
    fscanf(fp, "%lf", &DTl);
    fscanf(fp, "%d", &DTr);
    fscanf(fp, "%d", &DTp);
  } 
  else
  {
    fscanf(fp, "%lf", &buffer_double);
    fscanf(fp, "%lf", &buffer_double);
    fscanf(fp, "%lf", &buffer_double);
    fscanf(fp, "%lf", &buffer_double);
    fscanf(fp, "%d", &buffer_int);
    fscanf(fp, "%d", &buffer_int);
  }
  fscanf(fp, "%d", &restart_id);
  fscanf(fp, "%d", &dump_cnt);
  fscanf(fp, "%lf", &dt);
  #if METRIC==MKS
  fscanf(fp, "%lf", &Rin);
  fscanf(fp, "%lf", &Rout);
  fscanf(fp, "%lf", &a);
  fscanf(fp, "%lf", &hslope);
  fscanf(fp, "%lf", &Rhor);
  #else
  fscanf(fp, "%lf", &x1Min);
  fscanf(fp, "%lf", &x1Max);
  fscanf(fp, "%lf", &x2Min);
  fscanf(fp, "%lf", &x2Max);
  #endif

  // Read prims
  ZLOOPALL
  {
    PLOOP
      fscanf(fp, "%lf", &(S->P[ip][j][i]));
  }
}

/**
 * @brief Restart the simulation from the most recent checkpoint.
 *
 * @details Looks for `restarts/restart.last` (a symlink maintained by
 * restart_write()).  If not found, returns 0 and the caller initialises the
 * simulation from scratch.  If found:
 * 1. Calls zero_arrays() to clear pflag / fail_save.
 * 2. Calls restart_read() to fill S->P and global parameters.
 * 3. Calls set_grid() to (re-)compute metric quantities.
 * 4. Calls get_state_vec() to compute all 4-vectors from P.
 * 5. Calls prim_to_flux_vec() to fill S->U (conserved variables).
 * 6. Calls set_bounds() to apply boundary conditions.
 *
 * @param G  Grid geometry (overwritten by set_grid()).
 * @param S  Fluid state (overwritten by restart_read() and get_state_vec()).
 * @return   1 if a restart file was found and loaded; 0 otherwise.
 */
// Restarting a run by reading restart, setting grid, computing conserved quantities and applying boundary conditions
int restart_init(struct GridGeom *G, struct FluidState *S)
{
  char fname[STRLEN];
  sprintf(fname, "restarts/restart.last");

  FILE *fp = fopen(fname, "r");
  if (fp == NULL)
  {
    fprintf(stdout, "No restart file: error %d\n\n", errno);
    return 0;
  }
  fclose(fp);

  fprintf(stdout, "Loading restart file %s\n\n", fname);

  // Initialize flags and fails to zero
  zero_arrays();
  // Read the restart file
  restart_read(fname, S);
  // Set grid
  set_grid(G);
  // Compute four vectors
  get_state_vec(G, S, CENT, 0, N2-1, 0, N1-1);
  // Compute conserved quantities
  prim_to_flux_vec(G, S, 0, CENT, 0, N2-1, 0, N1-1, S->U);
  // Set boundary conditions
  set_bounds(G, S);

  return 1;
}
