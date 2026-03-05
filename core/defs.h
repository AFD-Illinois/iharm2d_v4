/**
 * @file defs.h
 * @brief Global variable definitions (storage allocation for all shared simulation state).
 *
 * @details This file is included by exactly one source file (main.c) via the
 * standard C pattern of declaring variables @c extern in the header (decs.h) and
 * defining them (allocating storage) here.
 *
 * The file defines:
 * - Per-zone integer flag arrays (pflag, fail_save, fflag).
 * - Physical parameters (a, gam, Rhor, tp_over_te).
 * - Geometry and grid parameters (Rin, Rout, hslope, dx, startx, ...).
 * - Timing and output cadence variables (DTd, DTl, DTr, DTp, ...).
 * - Diagnostic accumulators (mdot, edot, ldot, ...).
 * - Electron-heating parameters (game, gamp, fel0, ...) when ELECTRONS is defined.
 *
 * @note Do NOT include this file more than once; it allocates global storage.
 */

/*---------------------------------------------------------------------------------

  DEFS.H

  -GLOBAL VARIABLE DEFINITIONS

---------------------------------------------------------------------------------*/

#pragma once

// Zone flags
GridInt pflag;
GridInt fail_save;
GridInt fflag;

// Parameters
// physical
double a;
double gam;
double Rhor;
double tp_over_te;

// geometry
double Rin, Rout, hslope;
double poly_norm, poly_xt, poly_alpha, mks_smooth;
double cour;
double dV, dx[NDIM], startx[NDIM];
double x1Min, x1Max, x2Min, x2Max;
double dt, dt_light;
double t, tf;
double rcurr, hcurr;
int istart, istop, jstart, jstop;
int nstep;
int is_restart;

// fluid dumps
double DTd;
double DTf;
double DTl;
int DTr;
int DTp;
int dump_cnt;
double tdump, tlog;

// derived logged output
double mdot, edot, ldot;
double mdot_eh, edot_eh, ldot_eh;
int icurr, jcurr;

// Number of OpenMP threads
int nthreads;

// Electron-heating related globals
#if ELECTRONS
double game, gamp;
double fel0;
double tptemin, tptemax;
#endif
