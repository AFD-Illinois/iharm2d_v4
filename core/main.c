/*---------------------------------------------------------------------------------

  MAIN.C (Driver routine)
 
  -Read command line arguments
  -Call PARAMETERS.C and PROBLEM.C to read and set parameters
  -Memory allocation for grid and fluid state objects
  -Call RESTART.C to read restart file if it exists
  -Call PROBLEM.C to initialize problem if no restart file
  -Call DIAGNOSTICS.C to print initial diagnostics and write 0th dump
  -Loop in time, perform time-stepping by calling STEP.C

-----------------------------------------------------------------------------------*/

#include "decs.h"
#include "defs.h"
#include <time.h>
#include <sys/stat.h>

int main(int argc, char *argv[])
{
  fprintf(stdout, "\n          ************************************************************\n");
  fprintf(stdout, "          *                                                          *\n");
  fprintf(stdout, "          *                          IHARM2D_V4                      *\n");
  fprintf(stdout, "          *                                                          *\n");
  fprintf(stdout, "          *          Prather et al. JOSS, 2021 [UNDER REVIEW]        *\n");
  fprintf(stdout, "          *          Gammie, McKinney & Toth ApJ 509:444, 2003       *\n");
  fprintf(stdout, "          *                                                          *\n");
  fprintf(stdout, "          *  V K Dhruv                                               *\n");
  fprintf(stdout, "          *  B S Prather                                             *\n");
  fprintf(stdout, "          *  G N Wong                                                *\n");
  fprintf(stdout, "          *  C F Gammie                                              *\n");
  fprintf(stdout, "          *                                                          *\n");
  fprintf(stdout, "          *                          SYNTAX                          *\n");
  fprintf(stdout, "          *                                                          *\n");
  fprintf(stdout, "          *    -p /path/to/param.dat                                 *\n");
  fprintf(stdout, "          *    -o /path/to/output/dir                                *\n");
  fprintf(stdout, "          *                                                          *\n");
  fprintf(stdout, "          ************************************************************\n\n");

  // Read command line arguments
  char pfname[STRLEN] = "param.dat"; // Assumed parameter file name
  char outputdir[STRLEN] = ".";
  for (int n = 0; n < argc; n++)
  {
    if (*argv[n] == '-' && *(argv[n]+1) != '\0' && *(argv[n]+2) == '\0'&& n < argc-1)
    {
      if (*(argv[n]+1)=='o')
      {
        strcpy(outputdir, argv[++n]); // Set ouput directory path
      }
      if (*(argv[n]+1)=='p')
      {
        strcpy(pfname, argv[++n]); // Set parameter file path
      }
    }
  }

  // Read parameter file and set parameters
  set_core_params(); // initialize core params in paramtable dict
  set_problem_params(); // initialize problem-specific params in paramstable dict
  read_params(pfname); // read parameters from parameter file

  // Remove 'abort' file if it exists
  char abort_fname[256] = "abort";
  remove(abort_fname);

  // Change to output directory and make 'dumps' and 'restarts' subdirectories
  if (chdir(outputdir) != 0)
  {
    fprintf(stderr, "Output directory does not exist!\n");
    exit(2);
  }
  int is_error = mkdir("dumps/", 0777) || mkdir("restarts/", 0777); 
  if (is_error == -1 && errno != EEXIST)
  {
    fprintf(stderr, "Could not make dumps/restarts directory. Is output directory writeable?\n");
    exit(1);
  }

  // Set number of threads  
  #pragma omp parallel
  {
    #pragma omp master
    {
      nthreads = omp_get_num_threads();
    }
  }

  // Declare grid and fluid struct
  struct GridGeom *G = calloc(1,sizeof(struct GridGeom));
  struct FluidState *S = calloc(1,sizeof(struct FluidState));
  
  // Read restart file if it exists, else initialize problem
  is_restart = restart_init(G, S);
  if (!is_restart)
  {
    init_rng(1); // Standard seed for reproducibility
    init(G, S);
    nstep = 0;
    t = 0;
    dump_cnt = 0;
    zero_arrays();
    fprintf(stdout, "Initial conditions generated\n\n");
  }

  // In case we're restarting and these changed
  tdump = t + DTd;
  tlog = t + DTl;

  // Initial diagnostics
  diag(G, S, DIAG_INIT);
  if (!is_restart) restart_write(S);
  
  fprintf(stdout, "t= %e tf = %e\n", t, tf);

/*-----------------------------------------------------------------------------------
                                    Main loop
------------------------------------------------------------------------------------*/
  
  fprintf(stdout, "\nEntering main loop\n");

  // Initialize timers for processes to zero
  time_init();
  
  int dumpThisStep = 0;
  while (t < tf)
  {
    if (access(abort_fname, F_OK) != -1)
    {
      fprintf(stdout, "\nFound 'abort' file. Quitting now.\n\n");
			diag(G, S, DIAG_ABORT);
			restart_write_backend(S, IO_ABORT);
			return 0;
		}
  
    dumpThisStep = 0;
    // Start timer for time-step
    timer_start(TIMER_ALL);

    // Step forward in time
    step(G, S);
    nstep++;

    // Not evolving beyond tf
    if (t + dt > tf)
    {
      dt = tf-t;
    }

    fprintf(stdout, "t = %10.5g dt = %10.5g n = %d\n", t, dt, nstep);
  
    // Write dump, log or restart based on the respective cadence
    if (t < tf)
    {
      if (t >= tdump)
      {
        dumpThisStep = 1;
        diag(G, S, DIAG_DUMP);
        tdump += DTd;
      }
      if (t>=tlog)
      {
        diag(G, S, DIAG_LOG);
        tlog += DTl;
      }
      if (nstep % DTr == 0)
        restart_write(S);
    }
    
    // Stop timer for time-step
    timer_stop(TIMER_ALL);

    // Print code performance 
    if (nstep % DTp == 0)
      report_performance();
  }
/*-----------------------------------------------------------------------------------
                                    End main loop
------------------------------------------------------------------------------------*/

  if (dumpThisStep == 0)
    diag(G, S, DIAG_FINAL);
  
  return 0;
}
