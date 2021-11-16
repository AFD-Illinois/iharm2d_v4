/*---------------------------------------------------------------------------------

  PARAMETERS.H

  -COMPILE-TIME PARAMETERS

---------------------------------------------------------------------------------*/

// Global resolution
// Since there's no MPI, NiTOT==Ni
#define N1TOT 64
#define N2TOT 64


// Metric: 'MINKOWSKI''
#define METRIC MINKOWSKI

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
#define SUPPRESS_HIGHB_HEAT 1
#define BETA_HEAT           1

// Reconstruction algorithm: 'LINEAR' or 'PPM' or 'WENO' or 'MP5'

#define RECONSTRUCTION LINEAR

// Boundary conditions:
// 'OUTFLOW' for X1
// 'POLAR' for X2 
#define X1L_BOUND PERIODIC
#define X1R_BOUND PERIODIC
#define X2L_BOUND PERIODIC
#define X2R_BOUND PERIODIC
