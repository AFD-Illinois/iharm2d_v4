
/*---------------------------------------------------------------------------------

  PARAMETERS.H

  -COMPILE-TIME PARAMETERS

---------------------------------------------------------------------------------*/

// Global resolution
// Since there's no MPI, NiTOT==Ni
#define N1TOT 128
#define N2TOT 128


// Metric: 'MINKOWSKI' or 'MKS'
#define METRIC MKS
// Set 'DEREFINE_POLES' to 1 for FMKS
#define DEREFINE_POLES 1


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
 
#define ELECTRONS           1
#define ALLMODELS           1
#define SUPPRESS_HIGHB_HEAT 1
#define BETA_HEAT           1

// Reconstruction algorithm: 'LINEAR' or 'PPM' or 'WENO' or 'MP5'

#define RECONSTRUCTION WENO

// Boundary conditions:
// 'OUTFLOW' for X1
// 'POLAR' for X2 
#define X1L_BOUND OUTFLOW
#define X1R_BOUND OUTFLOW
#define X2L_BOUND POLAR
#define X2R_BOUND POLAR

#define X1L_INFLOW 0
#define X1R_INFLOW 0
