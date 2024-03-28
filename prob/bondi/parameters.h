/*---------------------------------------------------------------------------------

  PARAMETERS.H

  -COMPILE-TIME PARAMETERS

---------------------------------------------------------------------------------*/

// Global resolution
// Since there's no MPI, NiTOT==Ni
#define N1TOT 64
#define N2TOT 64


#define METRIC MKS
#define DEREFINE_POLES 0

// Floors:
// Wind term is a small source for torii only
// Maximum magnetization parameters should be set high for most problems
 
#define WIND_TERM 0
#define BSQORHOMAX (100.)
#define UORHOMAX (100.)

// Electrons and options:
// ALLMODELS - Flag for enabling all models
// SUPPRESS_MAG_HEAT - (0,1) No electron heating when sigma > 1
// BETA_HEAT         - (0,1) Beta-dependent electron heating
 
#define ELECTRONS           0
#define ALLMODELS           0

// Reconstruction algorithm: 'LINEAR' or 'WENO'
#define RECONSTRUCTION WENO

// Boundary conditions:
// 'OUTFLOW' for X1 inner boundary
// 'USER' for outer X1 boundary (analytic solution)
// 'POLAR' for X2 
#define X1L_BOUND OUTFLOW
#define X1R_BOUND USER
#define X2L_BOUND POLAR
#define X2R_BOUND POLAR

// Ensure no inflow at inner boundary
#define X1L_INFLOW 0
#define X1R_INFLOW 1
