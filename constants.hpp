/* vim: ft=cpp: */

#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP


//////////////////////////////////////////////
// Main switches to change the model in use //
//////////////////////////////////////////////

// Re-generates the results of the BG model from Tsirogiannis et al, 2010
//#define TSIROGIANNIS_2010 

// Simulates two distinct sources of cortical neurons: IT and PT
#define ITPTCTX 

// Simulates D1 and D2 MSN separatly
//#define MSN_SEPARATION 

// Simulates multi-channel nuclei
#define MULTICHANNELSMODEL




//////////////////////////////////////////////
// Minor options (convergence criteria etc) //
//////////////////////////////////////////////

// Change the integration timestep of the second order differential equation (default: 1e-4 seconds)
//#define ISSMALLDT /* dt = 1e-3s => one simulation takes less than 1 second */
//#define ISBIGDT /* dt = 1e-5s => one simulation takes about 20 seconds */
//#define ISHUGEDT /* dt = 1e-6s => one simulation takes 30+ minutes */

// Should we check twice for convergence, with the second check using a lower timestep?
#define CHECK_CONV_TWICE

// Should we bypass all convergence checks?
//#define BYPASS_CONV_CHECK

// When should we say that the convergence is attained?
#define TESTCONV
//#define LIGHTCONV
//#define SOUNDCONV
//#define OBJECTCONV

// Enables the contiguous check for convergence (only for single-channel nuclei)
#define SMALLECHCONV

// When should we measure the activity of nuclei? Set to lower values for quicker and less precise simulations
#define MINRECORDT 19
#define MAXRECORDT 20

// Disables the output of nuclei on std::cerr
//#define NOOUTPUT




/////////////////////////////////
// Definition of array indices //
// (change at your own risks!) //
/////////////////////////////////

enum BCBG_NUCLEI_INDICES {
  CTX_N = 0,
  CMPf_N,
  MSN_N,
// Set the indices for the simulation of special nuclei
#define MSND1_N MSN_N
#define FIRST_NUCLEUS_SIM MSN_N
  FSI_N,
  STN_N,
  GPe_N,
  GPi_N,

#ifdef MSN_SEPARATION
  // adds a nucleus to represent the D2 MSN neurons in striatum */
  MSND2_N,
#define MSND2_N MSND2
#define LAST_NUCLEUS_SIM (MSND2_N+1)
#else
#define MSND2_N MSND1_N
#define LAST_NUCLEUS_SIM (GPi_N+1)
#endif

#ifdef ITPTCTX
  // adds a nucleus to represent the PT population in cortex */
  CTXPT_N,
#endif

  NUCLEUS_NUMBER
};

enum TSIROGIANNIS_PARAMETERS {
  TSIROGIANNIS_2010_CTXe_D1_D2 = 0,
  TSIROGIANNIS_2010_CTXe_STN,
  TSIROGIANNIS_2010_STN_STN,
  TSIROGIANNIS_2010_GPe_STN,
  TSIROGIANNIS_2010_STN_GPe,
  TSIROGIANNIS_2010_D2_GPe,
  TSIROGIANNIS_2010_GPe_GPe,
  TSIROGIANNIS_2010_STN_GPi,
  TSIROGIANNIS_2010_D1_GPi,
  TSIROGIANNIS_2010_GPe_GPi,
  TSIROGIANNIS_2010_THETA_D1_D2,
  TSIROGIANNIS_2010_THETA_STN,
  TSIROGIANNIS_2010_THETA_GPe,
  TSIROGIANNIS_2010_THETA_GPi,
  TSIROGIANNIS_2010_PARAMS_COUNT
};

enum BCBG_PARAMS {
  CTX_MSN = 0, CTX_FSI, CTX_STN,
  MSN_GPe, MSN_GPi, STN_GPe, STN_GPi, STN_MSN, STN_FSI,
  GPe_STN, GPe_GPi, GPe_MSN, GPe_FSI, GPe_GPe,
  FSI_MSN, FSI_FSI,
  MSN_MSN,
  CMPf_MSN, CMPf_FSI, CMPf_STN, CMPf_GPe, CMPf_GPi,
#ifdef ITPTCTX
  CTXPT_MSN, CTXPT_FSI,
#endif
  DIST_CTX_MSN, 
#define CONNECT_NUMBER DIST_CTX_MSN
  DIST_CTX_FSI, DIST_CTX_STN,
  DIST_MSN_GPe, DIST_MSN_GPi,
  DIST_STN_GPe, DIST_STN_GPi, DIST_STN_MSN, DIST_STN_FSI,
  DIST_GPe_STN, DIST_GPe_GPi, DIST_GPe_MSN, DIST_GPe_FSI, DIST_GPe_GPe,
  DIST_FSI_MSN, DIST_FSI_FSI,
  DIST_MSN_MSN,
  DIST_CMPf_MSN, DIST_CMPf_FSI, DIST_CMPf_STN, DIST_CMPf_GPe, DIST_CMPf_GPi,
#ifdef ITPTCTX
  DIST_CTXPT_MSN, DIST_CTXPT_FSI,
#endif
  THETA_MSN, THETA_FSI, THETA_STN, THETA_GPe, THETA_GPi,
  FSI_SMAX,
  PARAMS_NUMBER
};

#define ARGS_NUMBER PARAMS_NUMBER

enum BCBG_DESACT {
  CTX_MSN_AMPA, CTX_FSI_AMPA, CTX_STN_AMPA,
  MSN_GPe_GABAA, MSN_GPi_GABAA,
  STN_GPe_AMPA, STN_GPi_AMPA, STN_MSN_AMPA, STN_FSI_AMPA,
  GPe_STN_GABAA, GPe_GPi_GABAA, GPe_MSN_GABAA, GPe_FSI_GABAA, GPe_GPe_GABAA,
  FSI_MSN_GABAA, FSI_FSI_GABAA,
  MSN_MSN_GABAA,
  CTX_MSN_NMDA, CTX_FSI_NMDA, CTX_STN_NMDA,
  STN_GPe_NMDA, STN_GPi_NMDA, STN_MSN_NMDA, STN_FSI_NMDA,
  MSN_GPe_GABAB, MSN_GPi_GABAB,
  GPe_STN_GABAB, GPe_GPi_GABAB, GPe_MSN_GABAB, GPe_FSI_GABAB, GPe_GPe_GABAB,
  FSI_MSN_GABAB, FSI_FSI_GABAB,
  MSN_MSN_GABAB,
  CMPf_MSN_AMPA, CMPf_FSI_AMPA, CMPf_STN_AMPA, CMPf_GPe_AMPA, CMPf_GPi_AMPA, CMPf_MSN_NMDA, CMPf_FSI_NMDA, CMPf_STN_NMDA, CMPf_GPe_NMDA, CMPf_GPi_NMDA,
  CTXPT_MSN_AMPA, CTXPT_FSI_AMPA, CTXPT_MSN_NMDA, CTXPT_FSI_NMDA,
  DESACT_NUMBER 
};


#endif
