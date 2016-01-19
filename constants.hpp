#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! //
// !! DO NOT FORGET TO UNCOMMENT ABOVE IF UNSET ITPTCTX !! //
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! //
#define ITPTCTX




#define TSIROGIANNIS_2010_CTXe_D1_D2 0
#define TSIROGIANNIS_2010_CTXe_STN 1
#define TSIROGIANNIS_2010_STN_STN 2
#define TSIROGIANNIS_2010_GPe_STN 3
#define TSIROGIANNIS_2010_STN_GPe 4
#define TSIROGIANNIS_2010_D2_GPe 5
#define TSIROGIANNIS_2010_GPe_GPe 6
#define TSIROGIANNIS_2010_STN_GPi 7
#define TSIROGIANNIS_2010_D1_GPi 8
#define TSIROGIANNIS_2010_GPe_GPi 9

#define TSIROGIANNIS_2010_THETA_D1_D2 10
#define TSIROGIANNIS_2010_THETA_STN 11
#define TSIROGIANNIS_2010_THETA_GPe 12
#define TSIROGIANNIS_2010_THETA_GPi 13

#define CTX_MSN 0
#define CTX_FSI 1
#define CTX_STN 2
#define MSN_GPe 3
#define MSN_GPi 4
#define STN_GPe 5
#define STN_GPi 6
#define STN_MSN 7
#define STN_FSI 8
#define GPe_STN 9
#define GPe_GPi 10
#define GPe_MSN 11
#define GPe_FSI 12
#define GPe_GPe 13
#define FSI_MSN 14
#define FSI_FSI 15
#define MSN_MSN 16
#define CMPf_MSN 17
#define CMPf_FSI 18
#define CMPf_STN 19
#define CMPf_GPe 20
#define CMPf_GPi 21
#define CTXPT_MSN 22
#define CTXPT_FSI 23


#ifdef ITPTCTX

#define DIST_CTX_MSN 24
#define DIST_CTX_FSI 25
#define DIST_CTX_STN 26
#define DIST_MSN_GPe 27
#define DIST_MSN_GPi 28
#define DIST_STN_GPe 29
#define DIST_STN_GPi 30
#define DIST_STN_MSN 31
#define DIST_STN_FSI 32
#define DIST_GPe_STN 33
#define DIST_GPe_GPi 34
#define DIST_GPe_MSN 35
#define DIST_GPe_FSI 36
#define DIST_GPe_GPe 37
#define DIST_FSI_MSN 38
#define DIST_FSI_FSI 39
#define DIST_MSN_MSN 40
#define DIST_CMPf_MSN 41
#define DIST_CMPf_FSI 42
#define DIST_CMPf_STN 43
#define DIST_CMPf_GPe 44
#define DIST_CMPf_GPi 45
#define DIST_CTXPT_MSN 46
#define DIST_CTXPT_FSI 47
#define THETA_MSN 48
#define THETA_FSI 49
#define THETA_STN 50
#define THETA_GPe 51
#define THETA_GPi 52
#define FSI_SMAX 53

#else

#define DIST_CTX_MSN 22
#define DIST_CTX_FSI 23
#define DIST_CTX_STN 24
#define DIST_MSN_GPe 25
#define DIST_MSN_GPi 26
#define DIST_STN_GPe 27
#define DIST_STN_GPi 28
#define DIST_STN_MSN 29
#define DIST_STN_FSI 30
#define DIST_GPe_STN 31
#define DIST_GPe_GPi 32
#define DIST_GPe_MSN 33
#define DIST_GPe_FSI 34
#define DIST_GPe_GPe 35
#define DIST_FSI_MSN 36
#define DIST_FSI_FSI 37
#define DIST_MSN_MSN 38
#define DIST_CMPf_MSN 39
#define DIST_CMPf_FSI 40
#define DIST_CMPf_STN 41
#define DIST_CMPf_GPe 42
#define DIST_CMPf_GPi 43
#define THETA_MSN 44
#define THETA_FSI 45
#define THETA_STN 46
#define THETA_GPe 47
#define THETA_GPi 48
#define FSI_SMAX 49

#endif



#define CTX_MSN_AMPA 0
#define CTX_FSI_AMPA 1
#define CTX_STN_AMPA 2
#define MSN_GPe_GABAA 3
#define MSN_GPi_GABAA 4
#define STN_GPe_AMPA 5
#define STN_GPi_AMPA 6
#define STN_MSN_AMPA 7
#define STN_FSI_AMPA 8
#define GPe_STN_GABAA 9
#define GPe_GPi_GABAA 10
#define GPe_MSN_GABAA 11
#define GPe_FSI_GABAA 12
#define GPe_GPe_GABAA 13
#define FSI_MSN_GABAA 14
#define FSI_FSI_GABAA 15
#define MSN_MSN_GABAA 16
#define CTX_MSN_NMDA 17
#define CTX_FSI_NMDA 18
#define CTX_STN_NMDA 19
#define STN_GPe_NMDA 20
#define STN_GPi_NMDA 21
#define STN_MSN_NMDA 22
#define STN_FSI_NMDA 23
#define MSN_GPe_GABAB 24
#define MSN_GPi_GABAB 25
#define GPe_STN_GABAB 26
#define GPe_GPi_GABAB 27
#define GPe_MSN_GABAB 28
#define GPe_FSI_GABAB 29
#define GPe_GPe_GABAB 30
#define FSI_MSN_GABAB 31
#define FSI_FSI_GABAB 32
#define MSN_MSN_GABAB 33
#define CMPf_MSN_AMPA 34
#define CMPf_FSI_AMPA 35
#define CMPf_STN_AMPA 36
#define CMPf_GPe_AMPA 37
#define CMPf_GPi_AMPA 38
#define CMPf_MSN_NMDA 39
#define CMPf_FSI_NMDA 40
#define CMPf_STN_NMDA 41
#define CMPf_GPe_NMDA 42
#define CMPf_GPi_NMDA 43
#define CTXPT_MSN_AMPA 44
#define CTXPT_FSI_AMPA 45
#define CTXPT_MSN_NMDA 46
#define CTXPT_FSI_NMDA 47


#define CTX_N 0
#define CMPf_N 1
#define MSN_N 2
#define FSI_N 3
#define STN_N 4
#define GPe_N 5
#define GPi_N 6

#define MSND1_N 2
#define MSND2_N 7

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! //
// !! DO NOT FORGET TO UNCOMMENT IF UNSET MSN_SEPARATION !! //
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! //
//#define CTXPT_N 8
#define CTXPT_N 7 
//// to use with msn not separated but itptctx

//#define TSIROGIANNIS_2010

//#define MSN_SEPARATION
//#define REDUCEDCMPF

#ifdef TSIROGIANNIS_2010

#define PARAMS_NUMBER 10
#define CONNECT_NUMBER 0
#define NUCLEUS_NUMBER 6
#define ARGS_NUMBER 14
#define DESACT_NUMBER -1

#elif defined(MSN_SEPARATION)

#ifdef ITPTCTX
// MSN_SEPARATION & ITPTCTX
#define PARAMS_NUMBER 54
#define CONNECT_NUMBER 24
#define NUCLEUS_NUMBER 9
#define ARGS_NUMBER 54
#define DESACT_NUMBER 48

#else
// MSN_SEPARATION
#define PARAMS_NUMBER 50
#define CONNECT_NUMBER 22
#define NUCLEUS_NUMBER 8
#define ARGS_NUMBER 50
#define DESACT_NUMBER 44

#endif

#else

#ifdef ITPTCTX
// ITPTCTX
#define PARAMS_NUMBER 54
#define CONNECT_NUMBER 24
#define NUCLEUS_NUMBER 8
#define ARGS_NUMBER 54
#define DESACT_NUMBER 48

#else
//
#define PARAMS_NUMBER 50
#define CONNECT_NUMBER 22
#define NUCLEUS_NUMBER 7
#define ARGS_NUMBER 50
#define DESACT_NUMBER 44

#endif

#endif


#define MULTICHANNELSMODEL

// various options regarding the integration timestep
//#define ISSMALLDT
//#define ISBIGDT
//#define ISHUGEDT
#define MIXEDDT
//#define MIXEDFULLDT
#define MIXEDMEDIUMDT

#define TRONQGALL

// when should we say that the convergence is attained?
#define TESTCONV
//#define LIGHTCONV
//#define SOUNDCONV
//#define OBJECTCONV
//#define FASTESTCONV
//#define FASTERCONV
//#define FASTCONV
#define SMALLECHCONV

#define NOTMUCHDELETEDCON

#define REALISTICDELAYS


#define GABAA025

#define RECTIF_NB_NEURONS

//#define MINRECORDT 5
//#define MAXRECORDT 5.09999
//#define ONLYMEANOUTPUT
//#define NOOUTPUT

#define MINRECORDT 19
//#define MINRECORDT 5
#define MAXRECORDT 20

#endif
