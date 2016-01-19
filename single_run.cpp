#include "bcbg2.hpp"

#include "constants.hpp"
#include "run_sim.hpp"

#include "helper_fct.hpp"
#include "stdlib.h"
#include <numeric>

void printout(std::vector <float> means, int ch_n)
{
  if (ch_n == 1) {
    std::cout << "  CTX = " << means[CTX_N] << std::endl;
    std::cout << "  CMPf = " << means[CMPf_N] << std::endl;
#ifdef MSN_SEPARATION
    std::cout << "  MSN D1 = " << means[MSND1_N] << std::endl;
    std::cout << "  MSN D2 = " << means[MSND2_N] << std::endl;
#else
    std::cout << "  MSN  = " << means[MSN_N] << std::endl;
#endif
    std::cout << "  FSI  = " << means[FSI_N] << std::endl;
    std::cout << "  STN  = " << means[STN_N] << std::endl;
    std::cout << "  GPe  = " << means[GPe_N] << std::endl;
    std::cout << "  GPi  = " << means[GPi_N] << std::endl<< std::endl;
  } else {
    std::cout << "  CTX = "; 
    for (int ch_i=0; ch_i<ch_n; ch_i++) {  std::cout << means[CTX_N*ch_n+ch_i] << " ; "; }
    std::cout << std::endl;

    std::cout << "  CMPf = "; 
    for (int ch_i=0; ch_i<ch_n; ch_i++) {  std::cout << means[CMPf_N*ch_n+ch_i] << " ; "; }
    std::cout << std::endl;

    std::cout << "  MSN = "; 
    for (int ch_i=0; ch_i<ch_n; ch_i++) {  std::cout << means[MSN_N*ch_n+ch_i] << " ; "; }
    std::cout << std::endl;

    std::cout << "  FSI = "; 
    for (int ch_i=0; ch_i<ch_n; ch_i++) {  std::cout << means[FSI_N*ch_n+ch_i] << " ; "; }
    std::cout << std::endl;

    std::cout << "  STN = "; 
    for (int ch_i=0; ch_i<ch_n; ch_i++) {  std::cout << means[STN_N*ch_n+ch_i] << " ; "; }
    std::cout << std::endl;

    std::cout << "  GPe = "; 
    for (int ch_i=0; ch_i<ch_n; ch_i++) {  std::cout << means[GPe_N*ch_n+ch_i] << " ; "; }
    std::cout << std::endl;

    std::cout << "  GPi = "; 
    for (int ch_i=0; ch_i<ch_n; ch_i++) {  std::cout << means[GPi_N*ch_n+ch_i] << " ; "; }
    std::cout << std::endl;
  }

}

int _is_near(
    float reference,
    float mean,
    float absolute_radius)
{
  if ((mean >= reference - absolute_radius) and (mean <= reference + absolute_radius)) {
    return 1;
  }
  return 0;
}

int _has_changed_near(
    float reference,
    float mean,
    float proportional_change,
    float proportional_radius)
{
  if ((mean >= reference * ( 1 + (proportional_change - proportional_radius)))
      && (mean <= reference * ( 1 + (proportional_change + proportional_radius)))) {
    return 1;
  }
  return 0;
}

float theta(float value)
{
  return ((float)(value/4.));
}

float pseudolog(float value)
{
  if (value < 100) {
    return (value); // 0, 1, ..., 99
  } else if (value < 200) {
    return ((100+(value-100)*10)); // 100, 110, ..., 990
  } else if (value < 300) {
    return ((1000+(value-200)*100)); // 1000, 1100, ..., 9900
  } else {
    std::cerr << "invalid parameter : " << value << std::endl;
    return (-1.); // throw an error
  }
}

void calc_score_1(std::vector <float> & params)
{
std::cout << "to be implemented" << std::endl;
std::cerr << "to be implemented" << std::endl;
}


int main(int argc, char** argv)
{
  int i,j;
#ifdef ISSMALLDT
  float dt = 1e-3;
#elif defined(ISBIGDT)
  float dt = 1e-5;
#elif defined(ISHUGEDT)
  float dt = 1e-6;
#else
#define CUSTOMDT 1.0f
  float dt = (1.0f/CUSTOMDT)*1e-4;
#endif
#ifdef MULTICHANNELSMODEL
    int ch_n = 8;
#elif defined(CELLSMODEL)
    int ch_n = 0;
#else
    int ch_n = 1;
#endif

  std::vector<float> params_synaptic;
  std::vector<int> params_delay;
  std::vector<float> means;
  std::vector<float> modulators_synaptic;
  std::vector<float> params_cs;


    int c = (int)(1./(dt*1000.0)+0.5);
    int im=0;
    params_delay.assign(ARGS_NUMBER,c);

    //////// delays given in Van Albada et al 2009
    //int delay_model_CTX_Str = 2;
    //int delay_model_CTX_STN = 1;
    //int delay_model_Str_GPe = 1;
    //int delay_model_Str_GPi = 1;
    //int delay_model_STN_GPe = 1;
    //int delay_model_STN_GPi = 1;
    //int delay_model_GPe_STN = 1;
    //int delay_model_GPe_GPi = 1;
    //int striatal_afferents = 1;
    
    ////// delays given in Tsirogiannis et al 2010
    int delay_model_CTX_Str = 4;
    int delay_model_CTX_STN = 1;
    int delay_model_Str_GPe = 3;
    int delay_model_Str_GPi = 3;
    int delay_model_STN_GPe = 1;
    int delay_model_STN_GPi = 1;
    int delay_model_GPe_STN = 1;
    int delay_model_GPe_GPi = 1;
    int delay_model_GPe_Str = 3;
    int delay_model_STN_Str = 3;
    int delay_model_CMPf    = 1;

    params_delay[CTX_MSN]  = c * delay_model_CTX_Str;
    params_delay[CTX_FSI]  = c * delay_model_CTX_Str;
    params_delay[CTX_STN]  = c * delay_model_CTX_STN;
    params_delay[CMPf_MSN] = c * 1;
    params_delay[CMPf_FSI] = c * 1;
    params_delay[CMPf_STN] = c * 1;
    params_delay[CMPf_GPe] = c * 1;
    params_delay[CMPf_GPi] = c * 1;
    params_delay[MSN_GPe]  = c * delay_model_Str_GPe;
    params_delay[MSN_GPi]  = c * delay_model_Str_GPi;
    params_delay[MSN_MSN]  = c * 1;
    params_delay[FSI_MSN]  = c * 1;
    params_delay[FSI_FSI]  = c * 1;
    params_delay[STN_GPe]  = c * delay_model_STN_GPe;
    params_delay[STN_GPi]  = c * delay_model_STN_GPi;
    params_delay[STN_MSN]  = c * delay_model_STN_Str;
    params_delay[STN_FSI]  = c * delay_model_STN_Str;
    params_delay[GPe_STN]  = c * delay_model_GPe_STN;
    params_delay[GPe_GPi]  = c * delay_model_GPe_GPi;
    params_delay[GPe_MSN]  = c * delay_model_GPe_Str;
    params_delay[GPe_FSI]  = c * delay_model_GPe_Str;
    params_delay[GPe_GPe]  = c * 1;
    params_delay[CTXPT_MSN]=   c * delay_model_CTX_STN;
    params_delay[CTXPT_FSI]=   c * delay_model_CTX_STN;

    std::vector<float> raw_params;
    raw_params.assign(ARGS_NUMBER,0.);


    float manual_input[ARGS_NUMBER];
    if (ch_n == 1) {
      std::cerr << std::endl;
    }
    for (int i=0; i<ARGS_NUMBER; i++) {
			manual_input[i] = atof(argv[i+1+3]); // +3 : the first three runs correspond to the run n°, score bio previously computed, score electro previously computed
      if (manual_input[i] < 0.0) {
        manual_input[i] = 0.0;
      } else if (manual_input[i] > 1.0) {
        manual_input[i] = 1.0;
      }
      if (ch_n == 1) {
        std::cerr << " " << manual_input[i] ;
      }
    }
    if (ch_n == 1) {
      std::cerr << std::endl;
    }

    for (int i=0; i<ARGS_NUMBER; i++) {raw_params[i] = manual_input[i];}

       raw_params[CTX_MSN]  = raw_params[CTX_MSN]*5995.0f+5.0f;//
       raw_params[CTX_FSI]  = raw_params[CTX_FSI]*5995.0f+5.0f;//Warning: these are not axonal boutons count, but synapse number
       raw_params[CTX_STN]  = raw_params[CTX_STN]*5995.0f+5.0f;//
       raw_params[CMPf_MSN] = param2boutons(raw_params[CMPf_MSN]);
       raw_params[CMPf_FSI] = param2boutons(raw_params[CMPf_FSI]);
       raw_params[CMPf_STN] = param2boutons(raw_params[CMPf_STN]);
       raw_params[CMPf_GPe] = param2boutons(raw_params[CMPf_GPe]);
       raw_params[CMPf_GPi] = param2boutons(raw_params[CMPf_GPi]);
       raw_params[MSN_GPe]  = param2boutons(raw_params[MSN_GPe]);
       raw_params[MSN_GPi]  = param2boutons(raw_params[MSN_GPi]);
       raw_params[MSN_MSN]  = param2boutons(raw_params[MSN_MSN]);
       raw_params[FSI_MSN]  = param2boutons(raw_params[FSI_MSN]);
       raw_params[FSI_FSI]  = param2boutons(raw_params[FSI_FSI]);
       raw_params[STN_GPe]  = param2boutons(raw_params[STN_GPe]);
       raw_params[STN_GPi]  = param2boutons(raw_params[STN_GPi]);
       raw_params[STN_MSN]  = param2boutons2(raw_params[STN_MSN]);
       raw_params[STN_FSI]  = param2boutons2(raw_params[STN_FSI]);
       raw_params[GPe_STN]  = param2boutons(raw_params[GPe_STN]);
       raw_params[GPe_GPi]  = param2boutons(raw_params[GPe_GPi]);
       raw_params[GPe_MSN]  = param2boutons2(raw_params[GPe_MSN]);
       raw_params[GPe_FSI]  = param2boutons2(raw_params[GPe_FSI]);
       raw_params[GPe_GPe]  = param2boutons(raw_params[GPe_GPe]);
       raw_params[CTXPT_MSN]= param2boutons(raw_params[CTXPT_MSN]);
       raw_params[CTXPT_FSI]= param2boutons(raw_params[CTXPT_FSI]);

    float score_0 = calc_score_selective_axons(raw_params,false,-1);

    // (factor 10³ for the numbers of neurons)
    float neurons_nb_CTX = 1400000; // total cortex (Christensen07, Collins10)
    float neurons_nb_CMPf = 86; // From Hunt91 and stereotaxic altases, this concerns only the non-gabaergic neurons of the CM/Pf. See count.odt for details of calculation
#ifdef MSN_SEPARATION
    float neurons_nb_MSN = 15200; // Yelnik91
#else
    float neurons_nb_MSN = 15200*2; // Yelnik91
#endif
    float neurons_nb_FSI = 611; // 2% (cf Deng10 or Yelnik91) of the total striatal count 30 400 (Yelnik91)
    float neurons_nb_STN = 77; // Hardman02: (STN)/2
    float neurons_nb_GPe = 251; // Hardman02: (GPe)/2
    float neurons_nb_GPi = 143; // Hardman02: (GPi + SN Non Dopaminergic)/2

#ifdef RECTIF_NB_NEURONS
    neurons_nb_MSN *= 0.87;
    neurons_nb_FSI *= 0.87;
#endif

    params_synaptic.assign(ARGS_NUMBER,0.);

    params_synaptic[CTX_MSN] = raw_params[CTX_MSN];
    params_synaptic[CTX_FSI] = raw_params[CTX_FSI];
    params_synaptic[CTX_STN] = raw_params[CTX_STN];
    params_synaptic[CMPf_MSN] = (1.  * raw_params[CMPf_MSN] * neurons_nb_CMPf) / (neurons_nb_MSN);
    params_synaptic[CMPf_FSI] = (1.  * raw_params[CMPf_FSI] * neurons_nb_CMPf) / (neurons_nb_FSI);
    params_synaptic[CMPf_STN] = (1.  * raw_params[CMPf_STN] * neurons_nb_CMPf) / (neurons_nb_STN);
    params_synaptic[CMPf_GPe] = (1.  * raw_params[CMPf_GPe] * neurons_nb_CMPf) / (neurons_nb_GPe);
    params_synaptic[CMPf_GPi] = (1.  * raw_params[CMPf_GPi] * neurons_nb_CMPf) / (neurons_nb_GPi);
    params_synaptic[MSN_GPe] = (1.   * raw_params[MSN_GPe] * neurons_nb_MSN) / (neurons_nb_GPe);
    params_synaptic[MSN_GPi] = (0.82 * raw_params[MSN_GPi] * neurons_nb_MSN) / (neurons_nb_GPi); // based on Levesque 2005
    params_synaptic[MSN_MSN] = (1.   * raw_params[MSN_MSN] * neurons_nb_MSN) / (neurons_nb_MSN);
    params_synaptic[FSI_MSN] = (1.   * raw_params[FSI_MSN] * neurons_nb_FSI) / (neurons_nb_MSN);
    params_synaptic[FSI_FSI] = (1.   * raw_params[FSI_FSI] * neurons_nb_FSI) / (neurons_nb_FSI);
    params_synaptic[STN_GPe] =1.0* (0.83 * raw_params[STN_GPe] * neurons_nb_STN) / (neurons_nb_GPe);
    params_synaptic[STN_GPi] =1.0* (0.72 * raw_params[STN_GPi] * neurons_nb_STN) / (neurons_nb_GPi);
    params_synaptic[STN_MSN] =1.0* (0.17 * raw_params[STN_MSN] * neurons_nb_STN) / (neurons_nb_MSN);
    params_synaptic[STN_FSI] =1.0* (0.17 * raw_params[STN_FSI] * neurons_nb_STN) / (neurons_nb_FSI);
    params_synaptic[GPe_STN] = (0.84 * raw_params[GPe_STN] * neurons_nb_GPe) / (neurons_nb_STN);
    params_synaptic[GPe_GPi] = 1.0*(0.84 * raw_params[GPe_GPi] * neurons_nb_GPe) / (neurons_nb_GPi);
    params_synaptic[GPe_MSN] = (0.16 * raw_params[GPe_MSN] * neurons_nb_GPe) / (neurons_nb_MSN);
    params_synaptic[GPe_FSI] = (0.16 * raw_params[GPe_FSI] * neurons_nb_GPe) / (neurons_nb_FSI);
    params_synaptic[GPe_GPe] = (1.   * raw_params[GPe_GPe] * neurons_nb_GPe) / (neurons_nb_GPe);
    params_synaptic[CTXPT_MSN] = raw_params[CTXPT_MSN];
    params_synaptic[CTXPT_FSI] = raw_params[CTXPT_FSI];

    params_synaptic[DIST_CTX_MSN]  = raw_params[DIST_CTX_MSN];
    params_synaptic[DIST_CTX_FSI]  = raw_params[DIST_CTX_FSI];
    params_synaptic[DIST_CTX_STN]  = raw_params[DIST_CTX_STN];
    params_synaptic[DIST_CMPf_MSN] = raw_params[DIST_CMPf_MSN];
    params_synaptic[DIST_CMPf_FSI] = raw_params[DIST_CMPf_FSI];
    params_synaptic[DIST_CMPf_STN] = raw_params[DIST_CMPf_STN];
    params_synaptic[DIST_CMPf_GPe] = raw_params[DIST_CMPf_GPe];
    params_synaptic[DIST_CMPf_GPi] = raw_params[DIST_CMPf_GPi];
    params_synaptic[DIST_MSN_GPe]  = raw_params[DIST_MSN_GPe];
    params_synaptic[DIST_MSN_GPi]  = raw_params[DIST_MSN_GPi];
    params_synaptic[DIST_MSN_MSN]  = raw_params[DIST_MSN_MSN];
    params_synaptic[DIST_FSI_MSN]  = raw_params[DIST_FSI_MSN];
    params_synaptic[DIST_FSI_FSI]  = raw_params[DIST_FSI_FSI];
    params_synaptic[DIST_STN_GPe]  = raw_params[DIST_STN_GPe];
    params_synaptic[DIST_STN_GPi]  = raw_params[DIST_STN_GPi];
    params_synaptic[DIST_STN_MSN]  = raw_params[DIST_STN_MSN];
    params_synaptic[DIST_STN_FSI]  = raw_params[DIST_STN_FSI];
    params_synaptic[DIST_GPe_STN]  = raw_params[DIST_GPe_STN];
    params_synaptic[DIST_GPe_GPi]  = raw_params[DIST_GPe_GPi];
    params_synaptic[DIST_GPe_MSN]  = raw_params[DIST_GPe_MSN];
    params_synaptic[DIST_GPe_FSI]  = raw_params[DIST_GPe_FSI];
    params_synaptic[DIST_GPe_GPe]  = raw_params[DIST_GPe_GPe];
    params_synaptic[DIST_CTXPT_MSN]  = raw_params[DIST_CTXPT_MSN];
    params_synaptic[DIST_CTXPT_FSI]  = raw_params[DIST_CTXPT_FSI];


    params_synaptic[THETA_MSN]  = raw_params[THETA_MSN];
    params_synaptic[THETA_FSI]  = raw_params[THETA_FSI];
    params_synaptic[THETA_STN]  = raw_params[THETA_STN];
    params_synaptic[THETA_GPe]  = raw_params[THETA_GPe];
    params_synaptic[THETA_GPi]  = raw_params[THETA_GPi];

    params_synaptic[FSI_SMAX]  = param2hz(raw_params[FSI_SMAX]);

    params_cs.assign(ARGS_NUMBER,0.);
    std::vector <float> cs;
    cs.assign(10,0.);

    // 0. => one-to-one
    // 1. => one-to-all
    cs[8] = 1.; // D* -> D*
    cs[0] = 0.; // CTX -> STN // no influence here
    cs[6] = 0.; // GPe -> D*  // according to the general consensus, but see Spooren et al 1996
    cs[3] = 1.; // STN -> D* // cf. Smith et al. 1990 (see most figures, and in particular figures 4 & 5)

    cs[7] = 1.; // GPe -> GPe // could be justified with Sato el al 200a
    cs[5] = 1.; // GPe -> GPi // !!!! Toggle this variable to check the focuse/diffused inhibition

    cs[4] = 0.; // GPe -> STN // Cf general consensus in rat & Sato et al 2000a 

    cs[2] = 1.; // STN -> GPi  // Cf general consensus in rat & Sato et al 2000a 
    cs[1] = 1.; // STN -> GPe // Cf general consensus in rat & Sato et al 2000a 

    params_cs[CTX_MSN] = 0.;
    params_cs[CTX_FSI] = 0.; // no influence
    params_cs[CTX_STN] = cs[0]; // no influence
    params_cs[MSN_GPe] = 0.; // for example, see Parent et al 1995c..
    params_cs[MSN_GPi] = 0.; // same
    params_cs[STN_GPe] = cs[1];
    params_cs[STN_GPi] = cs[2];
    params_cs[STN_MSN] = cs[3];
    params_cs[STN_FSI] = 1.; // cf Smith et al 1990 (see most figures, and in particular figures 4 & 5). Plus, there is no influence here
    params_cs[GPe_STN] = cs[4];
    params_cs[GPe_GPi] = cs[5];
    params_cs[GPe_MSN] = cs[6];
    params_cs[GPe_FSI] = 1.; // Spooren96
    params_cs[GPe_GPe] = cs[7];
    params_cs[FSI_MSN] = 1.; // Cf general consensus
    params_cs[FSI_FSI] = 1.; // Cf general consensus
    params_cs[MSN_MSN] = cs[8];
    params_cs[CMPf_MSN] = 1.; // Cf Parent04, but no influence here
    params_cs[CMPf_FSI] = 1.; // Cf Parent04, but no influence here
    params_cs[CMPf_STN] = 1.; // |
    params_cs[CMPf_GPe] = 1.; // | => cf. Sadikot et al 1992 (but note that there is no influence here)
    params_cs[CMPf_GPi] = 1.; // | 
    params_cs[CTXPT_MSN] = 0.; // it is consistent with Parent et al. 2006
    params_cs[CTXPT_FSI] = 0.; // no influence here

#ifdef MSN_SEPARATION
    int msn_separation = 1;
#else
    int msn_separation = 0;
#endif

    modulators_synaptic.assign(DESACT_NUMBER,1.0f);

    int result;

    means.assign(NUCLEUS_NUMBER*ch_n,0.); // note that we never store the value for the cortex (but we can display it)

    MemoryBCBG2 mem = {-1};

    // different convergence options are available as #define:
#ifdef FASTCONV
    float sim_time = 50;
#elif defined(FASTERCONV)
    float sim_time = 25;
#elif defined(LIGHTCONV)
    float sim_time = 0.5;
#elif defined(TESTCONV)
    float sim_time = 10;
#elif defined(SOUNDCONV)
    float sim_time = 5;
#elif defined(OBJECTCONV)
    float sim_time = 6;
#else
    float sim_time = 100;
#endif

    // // //
    // for the test of deactivations
    // // //
  float score_desact = 0;
  float score_desact_other = 0;

    // to get the logs
    // last and bef-bef-last do not matter in the case of multichannels nuclei
    // verbose does not matter if >4 and multichannels
    // 
    result = _run_sim(sim_time,0.001,dt,modulators_synaptic,params_cs,params_synaptic,params_delay,means,4,ch_n,msn_separation,0,mem,0); // verbose version

        for (int ch_i=0; ch_i < ch_n; ch_i++) {
          std::cout << means[MSN_N*ch_n+ch_i] << " ";
        }
        for (int ch_i=0; ch_i < ch_n; ch_i++) {
          std::cout << means[FSI_N*ch_n+ch_i] << " ";
        }
        for (int ch_i=0; ch_i < ch_n; ch_i++) {
          std::cout << means[STN_N*ch_n+ch_i] << " ";
        }
        for (int ch_i=0; ch_i < ch_n; ch_i++) {
          std::cout << means[GPe_N*ch_n+ch_i] << " ";
        }
        for (int ch_i=0; ch_i < ch_n; ch_i++) {
          std::cout << means[GPi_N*ch_n+ch_i] << " ";
        }


  score_desact_other = calc_score_desactivation(means, params_synaptic, params_delay, 0.0f, sim_time, mem,true);

    result = _run_sim(sim_time,0.001,dt,modulators_synaptic,params_cs,params_synaptic,params_delay,means,5,ch_n,msn_separation,0,mem,0); // verbose version
        for (int ch_i=0; ch_i < ch_n; ch_i++) {
          std::cout << means[MSN_N*ch_n+ch_i] << " ";
        }
        for (int ch_i=0; ch_i < ch_n; ch_i++) {
          std::cout << means[FSI_N*ch_n+ch_i] << " ";
        }
        for (int ch_i=0; ch_i < ch_n; ch_i++) {
          std::cout << means[STN_N*ch_n+ch_i] << " ";
        }
        for (int ch_i=0; ch_i < ch_n; ch_i++) {
          std::cout << means[GPe_N*ch_n+ch_i] << " ";
        }
        for (int ch_i=0; ch_i < ch_n; ch_i++) {
          std::cout << means[GPi_N*ch_n+ch_i] << " ";
        }

        std::cout << std::endl;
    return 0;

}



int main_tsirogiannis(int argc, char** argv)
{
  // re-implementation of Tsirogiannis et al 2010
  int i,j;
  float dt = 1e-5;
  int nb_channels = 1;

  std::vector<float> params_synaptic;
  std::vector<int> params_delay;
  std::vector<float> means;
  std::vector<float> modulators_synaptic;

  means.assign(NUCLEUS_NUMBER,0.);

  int c = (int)(1./(dt*1000.0)+0.5); // c corresponds to 1ms in timesteps
  int im=0;

  params_delay.assign(ARGS_NUMBER,c);
  params_synaptic.assign(ARGS_NUMBER,1);
  modulators_synaptic.assign(ARGS_NUMBER,1);

  params_synaptic[TSIROGIANNIS_2010_CTXe_D1_D2] = 50; params_delay[TSIROGIANNIS_2010_CTXe_D1_D2] = 4*c;
  params_synaptic[TSIROGIANNIS_2010_CTXe_STN] = 2; params_delay[TSIROGIANNIS_2010_CTXe_STN] = c;
  params_synaptic[TSIROGIANNIS_2010_STN_STN] = 2; params_delay[TSIROGIANNIS_2010_STN_STN] = c;
  params_synaptic[TSIROGIANNIS_2010_GPe_STN] = 10; params_delay[TSIROGIANNIS_2010_GPe_STN] = c;
  params_synaptic[TSIROGIANNIS_2010_STN_GPe] = 3; params_delay[TSIROGIANNIS_2010_STN_GPe] = c;
  params_synaptic[TSIROGIANNIS_2010_D2_GPe] = 33; params_delay[TSIROGIANNIS_2010_D2_GPe] = 3*c;
  params_synaptic[TSIROGIANNIS_2010_GPe_GPe] = 1; params_delay[TSIROGIANNIS_2010_GPe_GPe] = c;
  params_synaptic[TSIROGIANNIS_2010_STN_GPi] = 3; params_delay[TSIROGIANNIS_2010_STN_GPi] = c;
  params_synaptic[TSIROGIANNIS_2010_D1_GPi] = 22; params_delay[TSIROGIANNIS_2010_D1_GPi] = 3*c;
  params_synaptic[TSIROGIANNIS_2010_GPe_GPi] = 5; params_delay[TSIROGIANNIS_2010_GPe_GPi] = c;

  params_synaptic[TSIROGIANNIS_2010_THETA_D1_D2]  = 27.;
  params_synaptic[TSIROGIANNIS_2010_THETA_STN]  = 18.5;
  params_synaptic[TSIROGIANNIS_2010_THETA_GPe]  = 14.;
  params_synaptic[TSIROGIANNIS_2010_THETA_GPi]  = 12.;

  int result = _run_sim_tsirogiannis_2010(1000,0.001,dt,modulators_synaptic,params_synaptic,params_delay,means,2); //output

  std::cout << "At rest" << std::endl;
  std::cout << "  D1  = " << means[0] << std::endl;
  std::cout << "  D2  = " << means[1] << std::endl;
  std::cout << "  STN = " << means[2] << std::endl;
  std::cout << "  GPe = " << means[3] << std::endl;
  std::cout << "  GPi = " << means[4] << std::endl;
  std::cout << "  CTXe = " << means[5] << std::endl;

  if (result) {
    std::cout << "CONVERGENCE OK" << std::endl;
  } else {
    std::cout << "PAS DE CONVERGENCE" << std::endl;
  }

  return 0;
}

