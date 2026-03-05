/**
 * @file timing.c
 * @brief Wall-clock performance timing and reporting.
 *
 * @details Maintains two static arrays:
 * - `timers[NUM_TIMERS]`: stores the wall time recorded at the last timer_start() call.
 * - `times[NUM_TIMERS]`: accumulates total elapsed wall time per timer.
 *
 * Usage pattern:
 * ```c
 * timer_start(TIMER_RECON);
 * // ... work to be timed ...
 * timer_stop(TIMER_RECON);
 * ```
 *
 * Only the OpenMP master thread updates the timing arrays; all calls inside
 * parallel regions are safe due to the `#pragma omp master` guard.
 * When the TIMERS flag is disabled, only the TIMER_ALL entry is still updated.
 *
 * report_performance() prints a per-timer breakdown (wall seconds per step +
 * percentage of TIMER_ALL) and the zone-cycles-per-second throughput.
 *
 * @note The file header incorrectly says "METRIC.C"; the file actually
 * implements timing functionality for the whole simulation.
 */

/*---------------------------------------------------------------------------------

  METRIC.C

  -Performance timing and reporting

---------------------------------------------------------------------------------*/

#include "decs.h"

static double timers[NUM_TIMERS];
static double times[NUM_TIMERS];

static int nstep_start = 0;

/**
 * @brief Initialize (reset) all accumulated timing counters to zero.
 *
 * @details Sets all entries in `times[]` to 0.0 and records the current
 * nstep so that report_performance() can compute a per-step average.
 * Must be called once after the first step is complete and the warm-up
 * overhead should be excluded.
 */
void time_init()
{
  for (int n = 0; n < NUM_TIMERS; n++) {
    times[n] = 0.;
  }
  nstep_start = nstep;
}

/**
 * @brief Record the start time for a named timer.
 *
 * @details Stores omp_get_wtime() in `timers[timerCode]`.
 * Only executed by the OpenMP master thread (harmless to call inside or outside
 * parallel regions).  No-op unless TIMERS is set or timerCode == TIMER_ALL.
 *
 * @param timerCode  Timer identifier (one of the TIMER_* constants in decs.h).
 */
inline void timer_start(int timerCode)
{
  if (TIMERS || timerCode == TIMER_ALL) {
#pragma omp master
    {
      timers[timerCode] = omp_get_wtime();
    }
  }
}

/**
 * @brief Accumulate elapsed time for a named timer.
 *
 * @details Adds `omp_get_wtime() - timers[timerCode]` to `times[timerCode]`.
 * Must be called after the corresponding timer_start().
 * Only executed by the OpenMP master thread.
 *
 * @param timerCode  Timer identifier matching a prior timer_start() call.
 */
inline void timer_stop(int timerCode)
{
  if (TIMERS || timerCode == TIMER_ALL) {
#pragma omp master
    {
      times[timerCode] += (omp_get_wtime() - timers[timerCode]);
    }
  }
}

/**
 * @brief Print a summary of timing statistics to stdout.
 *
 * @details Computes the number of steps elapsed since time_init(), then for
 * each timer prints the per-step wall time and its percentage of TIMER_ALL.
 * Also reports:
 * - Zone-cycles per core-second: N1*N2 / (T_ALL * nthreads / nsteps).
 * - Zone-cycles per second:      N1*N2 / (T_ALL / nsteps).
 * Only the timer breakdowns are gated by the TIMERS compile-time flag;
 * the final two lines are always printed.
 */
// Report a running average of performance data
void report_performance()
{
  int steps = nstep - nstep_start;
#if TIMERS
  fprintf(stdout, "\n********** PERFORMANCE **********\n");
  fprintf(stdout, "   RECON:    %8.4g s (%.4g %%)\n",
    times[TIMER_RECON]/steps, 100.*times[TIMER_RECON]/times[TIMER_ALL]);
  fprintf(stdout, "   LR_TO_F:  %8.4g s (%.4g %%)\n",
    times[TIMER_LR_TO_F]/steps, 100.*times[TIMER_LR_TO_F]/times[TIMER_ALL]);
  fprintf(stdout, "   CMAX:     %8.4g s (%.4g %%)\n",
    times[TIMER_CMAX]/steps, 100.*times[TIMER_CMAX]/times[TIMER_ALL]);
  fprintf(stdout, "   FLUX_CT:  %8.4g s (%.4g %%)\n",
    times[TIMER_FLUX_CT]/steps, 100.*times[TIMER_FLUX_CT]/times[TIMER_ALL]);
  fprintf(stdout, "   UPDATE_U: %8.4g s (%.4g %%)\n",
    times[TIMER_UPDATE_U]/steps, 100.*times[TIMER_UPDATE_U]/times[TIMER_ALL]);
  fprintf(stdout, "   U_TO_P:   %8.4g s (%.4g %%)\n",
    times[TIMER_U_TO_P]/steps, 100.*times[TIMER_U_TO_P]/times[TIMER_ALL]);
  fprintf(stdout, "   FIXUP:    %8.4g s (%.4g %%)\n",
    times[TIMER_FIXUP]/steps, 100.*times[TIMER_FIXUP]/times[TIMER_ALL]);
  fprintf(stdout, "   BOUND:    %8.4g s (%.4g %%)\n",
    times[TIMER_BOUND]/steps, 100.*times[TIMER_BOUND]/times[TIMER_ALL]);
  fprintf(stdout, "   BOUND_COMMS:    %8.4g s (%.4g %%)\n",
    times[TIMER_BOUND_COMMS]/steps, 100.*times[TIMER_BOUND_COMMS]/times[TIMER_ALL]);
  fprintf(stdout, "   DIAG:     %8.4g s (%.4g %%)\n",
    times[TIMER_DIAG]/steps, 100.*times[TIMER_DIAG]/times[TIMER_ALL]);
  fprintf(stdout, "   IO:     %8.4g s (%.4g %%)\n",
    times[TIMER_IO]/steps, 100.*times[TIMER_IO]/times[TIMER_ALL]);
  fprintf(stdout, "   RESTART:     %8.4g s (%.4g %%)\n",
    times[TIMER_RESTART]/steps, 100.*times[TIMER_RESTART]/times[TIMER_ALL]);
  fprintf(stdout, "   CURRENT:     %8.4g s (%.4g %%)\n",
    times[TIMER_CURRENT]/steps, 100.*times[TIMER_CURRENT]/times[TIMER_ALL]);
  fprintf(stdout, "   LR_STATE:     %8.4g s (%.4g %%)\n",
    times[TIMER_LR_STATE]/steps, 100.*times[TIMER_LR_STATE]/times[TIMER_ALL]);
  fprintf(stdout, "   LR_PTOF:     %8.4g s (%.4g %%)\n",
    times[TIMER_LR_PTOF]/steps, 100.*times[TIMER_LR_PTOF]/times[TIMER_ALL]);
  fprintf(stdout, "   LR_VCHAR:     %8.4g s (%.4g %%)\n",
    times[TIMER_LR_VCHAR]/steps, 100.*times[TIMER_LR_VCHAR]/times[TIMER_ALL]);
  fprintf(stdout, "   LR_CMAX:     %8.4g s (%.4g %%)\n",
    times[TIMER_LR_CMAX]/steps, 100.*times[TIMER_LR_CMAX]/times[TIMER_ALL]);
  fprintf(stdout, "   LR_FLUX:     %8.4g s (%.4g %%)\n",
    times[TIMER_LR_FLUX]/steps, 100.*times[TIMER_LR_FLUX]/times[TIMER_ALL]);
#if ELECTRONS
  fprintf(stdout, "   E_HEAT:   %8.4g s (%.4g %%)\n",
    times[TIMER_ELECTRON_HEAT]/steps,
    100.*times[TIMER_ELECTRON_HEAT]/times[TIMER_ALL]);
  fprintf(stdout, "   E_FIXUP:  %8.4g s (%.4g %%)\n",
    times[TIMER_ELECTRON_FIXUP]/steps,
    100.*times[TIMER_ELECTRON_FIXUP]/times[TIMER_ALL]);
#endif
#endif

  fprintf(stdout, "   ALL:      %8.4g s\n", times[TIMER_ALL]/steps);
  fprintf(stdout, "   ZONE CYCLES PER\n");
  fprintf(stdout, "     CORE-SECOND: %e\n",
    N1*N2/(times[TIMER_ALL]*nthreads/steps));
  fprintf(stdout, "          SECOND: %e\n",
        N1*N2/(times[TIMER_ALL]/steps));
}
