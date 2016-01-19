#include "constants.hpp"
#include "bcbg2.hpp"
#include "helper_fct.hpp"
#include "run_sim.hpp"


int is_in(std::vector <int>& c, int i)
{
  for (int j=0; j<c.size(); j++) {
    if (i==c[j]) {
      return 1;
    }
  }
  return 0;
}

float param2hz(float value) {
  return 200.0 + value * 300.0f;
}

float param2boutons(float value) {
#ifdef NODELETEDCON
  return value * 5995.0f + 5.0f;
#elif defined(NOTMUCHDELETEDCON)
  return value * 5995.0f + 5.0f;
#else
  return value * 6000.0f;
#endif
}

float param2boutons2(float value) {
#ifdef NOTMUCHDELETEDCON
  return value * 6000.0f;
#elif defined(NODELETEDCON)
  return value * 5995.0f + 5.0f;
#else
  return value * 6000.0f;
#endif
}

float _do_trial(
    std::vector<float>& means,
    std::vector<float>& ref,
    std::vector<float>& params,
    std::vector<int>& delays,
    std::vector<float>& activations,
    int nucleus,
    float proportional_change,
    float proportional_radius,
    float sim_time,
    MemoryBCBG2& mem,
    bool checkedtrial)
{
  // this functions performs one trial of a deactivation
  // there are many options with #define implemented here to allow various speed-up of the optimization
#ifdef MSN_SEPARATION
  int msn_separation = 1;
#else
  int msn_separation = 0;
#endif

#ifdef MULTICHANNELSMODEL
    int ch_n = 3;
#elif defined(CELLSMODEL)
    int ch_n = 0;
#else
    int ch_n = 1;
#endif

#ifdef ISSMALLDT
  float dt = 1e-3;
#elif defined(ISBIGDT)
  float dt = 1e-5;
#elif defined(ISHUGEDT)
  float dt = 1e-6;
#else
  float dt = 1e-4;
#endif

  std::vector<float> cs;
  cs.assign(ARGS_NUMBER, 0.);

#ifdef LIGHTESTCONV
  float sim_step = 0.01;
#else
  float sim_step = 0.1;
#endif

#if defined(MIXEDMEDIUMDT)
  int converged;
  if (checkedtrial) {
    converged = _run_sim(sim_time,sim_step,dt,activations,cs,params,delays,means,0,ch_n,msn_separation,0,mem,1);
  } else {
    converged = _run_sim(sim_time,sim_step,dt,activations,cs,params,delays,means,0,ch_n,msn_separation,0,mem,0);
    if (converged == 1) {
      std::vector <float> smalldt_means(means);
      converged = _run_sim(sim_time,sim_step,dt,activations,cs,params,delays,means,0,ch_n,msn_separation,0,mem,1);
      if (converged != -1) {
        for (int i=0;i<NUCLEUS_NUMBER;i++) {
          if (abs(means[i*ch_n] - smalldt_means[i*ch_n]) > 1) {
            converged = -1;
          }
        }
      }
    }
  }
#elif defined(MIXEDFULLDT)
  delays.assign(ARGS_NUMBER,1); // 1ms everywhere
  int converged = _run_sim(sim_time,0.1,1e-3,activations,cs,params,delays,means,0,ch_n,msn_separation,0,mem,0);
  if (converged != -1) {
    std::vector <float> smalldt_means(means);
    delays.assign(ARGS_NUMBER,10); // 1ms everywhere
    converged = _run_sim(sim_time,0.1,1e-4,activations,cs,params,delays,means,0,ch_n,msn_separation,0,mem,0);
    if (converged != -1) {
      float diff = 0.;
      for (int i=0;i<NUCLEUS_NUMBER;i++) {
        diff += abs(means[i*ch_n] - smalldt_means[i*ch_n]);
      }
      if (diff > 1) {
        converged = 0;
      }
    }
  }
#else
  int converged = _run_sim(sim_time,sim_step,dt,activations,cs,params,delays,means,0,ch_n,msn_separation,0,mem,0);
#endif

#ifdef TRONQGALL
  float score = _has_changed_near_tronqgaussian(ref[nucleus*ch_n],means[nucleus*ch_n],proportional_change,proportional_radius);
#else
  float score = _has_changed_near_gaussian(ref[nucleus*ch_n],means[nucleus*ch_n],proportional_change,proportional_radius);
#endif

  if ( converged == -1) {
    return -10000;
  } else {
    return score;
  }

}



float calc_score_desactivation(
    std::vector <float>& means,
    std::vector <float>& params,
    std::vector <int>& delays,
    float desactivation_level,
    float sim_time,
    MemoryBCBG2& mem,
    bool verbose)
{
  // computes part of the electrophysiological plausibility objective
  // these are the deactivation studies taken into account in the paper Lienard and Girard 2013.
#ifdef MULTICHANNELSMODEL
    int ch_n = 3;
#elif defined(CELLSMODEL)
    int ch_n = 0;
#else
    int ch_n = 1;
#endif

  float scoreg = 0.;
  float score = 0.;

  std::vector <float> activations;
  activations.assign(DESACT_NUMBER,1.0f);

  std::vector <float> ref(means);
  std::vector <float> ref_tmp(ref);

  // when NBQX is applied in GPe, we expect a 56.7% decrease in GPe
  // (cf Kita et al. 2004 : "Role of Ionotropic Glutamatergic and GABAergic Inputs
  // on the Firing Activity of Neurons in the External Pallidum in Awake Monkeys")
  activations.assign(DESACT_NUMBER,1);
  activations[STN_GPe_AMPA] = desactivation_level; // to block the AMPA channel of GPe
  activations[CMPf_GPe_AMPA] = desactivation_level; // to block the AMPA channel of GPe
  score = _do_trial(means,ref,params,delays,activations,GPe_N,-0.567,0.356,sim_time,mem,false);
  scoreg += score;
  if (verbose) {
    std::cout << means[GPe_N*ch_n] << " ";
  }

  // WARNING : this desactivation has to follow the NBQX application in GPe.
  ref_tmp.assign(NUCLEUS_NUMBER,0.);
  ref_tmp[GPe_N*ch_n] = means[GPe_N*ch_n];
  // when NBQX then gabazine are applied in GPe, we expect a 116.5% further increase in the GPe
  // (cf Kita et al. 2004 : "Role of Ionotropic Glutamatergic and GABAergic Inputs
  // on the Firing Activity of Neurons in the External Pallidum in Awake Monkeys")
  activations.assign(DESACT_NUMBER,1);
  activations[STN_GPe_AMPA] = desactivation_level; // to block the AMPA channel of GPe
  activations[CMPf_GPe_AMPA] = desactivation_level; // to block the AMPA channel of GPe
  activations[MSN_GPe_GABAA] = desactivation_level; // to block the GABA channel of GPe (from D2)
  activations[GPe_GPe_GABAA] = desactivation_level; // to block the GABA channel of GPe (from GPe)
  if (scoreg >= 0) {
    score = _do_trial(means,ref_tmp,params,delays,activations,GPe_N,+1.165,0.167,sim_time,mem,false);
    scoreg += score;
  } else if (verbose) {
    score = _do_trial(means,ref_tmp,params,delays,activations,GPe_N,+1.165,0.167,sim_time,mem,false);
  }
  if (verbose) {
    std::cout << means[GPe_N*ch_n] << " ";
  }

  // when CPP is applied in GPe, we expect a 32.4% decrease in GPe
  // (cf Kita et al. 2004 : "Role of Ionotropic Glutamatergic and GABAergic Inputs
  // on the Firing Activity of Neurons in the External Pallidum in Awake Monkeys")
  activations.assign(DESACT_NUMBER,1);
  activations[STN_GPe_NMDA] = desactivation_level; // to block the NMDA channel of GPe
  activations[CMPf_GPe_NMDA] = desactivation_level; // to block the NMDA channel of GPe
  if (scoreg >= 0) {
    score = _do_trial(means,ref,params,delays,activations,GPe_N,-0.324,0.145,sim_time,mem,false);
    scoreg += score;
  } else if (verbose) {
    score = _do_trial(means,ref,params,delays,activations,GPe_N,-0.324,0.145,sim_time,mem,false);
  }
  if (verbose) {
    std::cout << means[GPe_N*ch_n] << " ";
  }

  // when gabazine is applied in GPe, we expect a 115.8% increase in GPe
  // (cf Kita et al. 2004 : "Role of Ionotropic Glutamatergic and GABAergic Inputs
  // on the Firing Activity of Neurons in the External Pallidum in Awake Monkeys")
  activations.assign(DESACT_NUMBER,1);
  activations[MSN_GPe_GABAA] = desactivation_level; // to block the GABA channel of GPe (from D2)
  activations[GPe_GPe_GABAA] = desactivation_level; // to block the GABA channel of GPe (from GPe)
  if (scoreg >= 0) {
    score = _do_trial(means,ref,params,delays,activations,GPe_N,+1.158,0.815,sim_time,mem,false);
    scoreg += score;
  } else if (verbose) {
    score = _do_trial(means,ref,params,delays,activations,GPe_N,+1.158,0.815,sim_time,mem,false);
  }
  if (verbose) {
    std::cout << means[GPe_N*ch_n] << " ";
  }



  // when CPP is applied in GPi, we expect a 27.5% decrease in GPi
  // (cf Tachibana et al. 2008 : "Motor cortical control of internal pallidal
  // activity through glutamatergic and GABAergic inputs in awake monkeys")
  activations.assign(DESACT_NUMBER,1);
  activations[STN_GPi_NMDA] = desactivation_level; // to block the NMDA channel of GPi
  activations[CMPf_GPi_NMDA] = desactivation_level; // to block the NMDA channel of GPi
  if (scoreg >= 0) {
    score = _do_trial(means,ref,params,delays,activations,GPi_N,-0.275,0.264,sim_time,mem,false);
    scoreg += score;
  } else if (verbose) {
    score = _do_trial(means,ref,params,delays,activations,GPi_N,-0.275,0.264,sim_time,mem,false);
  }
  if (verbose) {
    std::cout << means[GPi_N*ch_n] << " ";
  }

  // WARNING : this desactivation has to follow the CPP application in GPi.
  ref_tmp.assign(NUCLEUS_NUMBER,0.);
  ref_tmp[GPi_N*ch_n] = means[GPi_N*ch_n];
  // when CPP and then NBQX are applied in GPi, we expect a 54.2% further decrease in GPi
  // (cf Tachibana et al. 2008 : "Motor cortical control of internal pallidal
  // activity through glutamatergic and GABAergic inputs in awake monkeys")
  activations.assign(DESACT_NUMBER,1);
  activations[STN_GPi_AMPA] = desactivation_level; // to block the AMPA channel of GPi
  activations[CMPf_GPi_AMPA] = desactivation_level; // to block the AMPA channel of GPi
  activations[STN_GPi_NMDA] = desactivation_level; // to block the NMDA channel of GPi
  activations[CMPf_GPi_NMDA] = desactivation_level; // to block the NMDA channel of GPi
  if (scoreg >= 0) {
    score = _do_trial(means,ref_tmp,params,delays,activations,GPi_N,-0.542,0.208,sim_time,mem,false);
    scoreg += score;
  } else if (verbose) {
    score = _do_trial(means,ref_tmp,params,delays,activations,GPi_N,-0.542,0.208,sim_time,mem,false);
  }
  if (verbose) {
    std::cout << means[GPi_N*ch_n] << " ";
  }

  // when NBQX is applied in GPi, we expect a 53.6% decrease in GPi
  // (cf Tachibana et al. 2008 : "Motor cortical control of internal pallidal
  // activity through glutamatergic and GABAergic inputs in awake monkeys")
  activations.assign(DESACT_NUMBER,1);
  activations[STN_GPi_AMPA] = desactivation_level; // to block the AMPA channel of GPi
  activations[CMPf_GPi_AMPA] = desactivation_level; // to block the AMPA channel of GPi
  if (scoreg >= 0) {
    score = _do_trial(means,ref,params,delays,activations,GPi_N,-0.536,0.367,sim_time,mem,false);
    scoreg += score;
  } else if (verbose) {
    score = _do_trial(means,ref,params,delays,activations,GPi_N,-0.536,0.367,sim_time,mem,false);
  }
  if (verbose) {
    std::cout << means[GPi_N*ch_n] << " ";
  }

  // when gabazine is applied in GPi, we expect a 92% increase in GPi
  // (cf Tachibana et al. 2008 : "Motor cortical control of internal pallidal
  // activity through glutamatergic and GABAergic inputs in awake monkeys")
  activations.assign(DESACT_NUMBER,1);
  activations[MSN_GPi_GABAA] = desactivation_level; // to block the GABA channel of GPi (from D1)
  activations[GPe_GPi_GABAA] = desactivation_level; // to block the GABA channel of GPi (from GPe)
  if (scoreg >= 0) {
    score = _do_trial(means,ref,params,delays,activations,GPi_N,+0.92,1.173,sim_time,mem,false);
    scoreg += score;
  } else if (verbose) {
    score = _do_trial(means,ref,params,delays,activations,GPi_N,+0.92,1.173,sim_time,mem,false);
  }
  if (verbose) {
    std::cout << means[GPi_N*ch_n] << " ";
  }

  // this desactivation was unset before. Reason ? Because the data is not supplied with changes, but as an absolute value. The choice is to integrate it with the absolute value
  // when CPP, NBQX and gabazine are applied in GPi, we expect a normal firing rate
  // (ie : near 75.1Hz)
  // (cf Tachibana et al. 2008 : "Motor cortical control of internal pallidal
  // activity through glutamatergic and GABAergic inputs in awake monkeys")
  activations.assign(DESACT_NUMBER,1);
  activations[STN_GPi_AMPA] = desactivation_level; // to block the AMPA channel of GPi (from STN)
  activations[CMPf_GPi_AMPA] = desactivation_level; // to block the AMPA channel of GPi (from CM/Pf)
  activations[STN_GPi_NMDA] = desactivation_level; // to block the NMDA channel of GPi (from STN)
  activations[CMPf_GPi_NMDA] = desactivation_level; // to block the NMDA channel of GPi (from CM/Pf)
  activations[MSN_GPi_GABAA] = desactivation_level; // to block the GABA channel of GPi (from D1)
  activations[GPe_GPi_GABAA] = desactivation_level; // to block the GABA channel of GPi (from GPe)
  ref_tmp[GPi_N*ch_n] = 75.1;
  if (scoreg >= 0) {
    score = _do_trial(means,ref_tmp,params,delays,activations,GPi_N,0.0f,0.217,sim_time,mem,false);
    scoreg += score;
  } else if (verbose) {
    score = _do_trial(means,ref_tmp,params,delays,activations,GPi_N,0.0f,0.217,sim_time,mem,false);
  }
  if (verbose) {
    std::cout << means[GPi_N*ch_n] << " ";
  }

  if (scoreg < 0) {
    scoreg = 0;
  }

  return scoreg;
}


float calc_score_desactivation_other(
    std::vector <float>& means,
    std::vector <float>& params,
    std::vector <int>& delays,
    float desactivation_level,
    float sim_time,
    MemoryBCBG2& mem,
    bool verbose)
{
  // these deactivation studies were included later in the process of model making, and are not part of the initial constraints for the optimization
#ifdef MULTICHANNELSMODEL
    int ch_n = 3;
#elif defined(CELLSMODEL)
    int ch_n = 0;
#else
    int ch_n = 1;
#endif

  float scoreg = 0.;
  float score = 0.;

  std::vector <float> activations;
  activations.assign(DESACT_NUMBER,1.0f);

  std::vector <float> ref(means);
  std::vector <float> ref_tmp(ref);

  // NB for GPe: "The STN blockade greatly decreased the firing rate, to complete silence in some neurons. However, 5–10 min after the muscimol injection, the activity began to increase with repeated occurrences of short grouped spike discharges. As time progressed, the activity further increased and developed into repeated occurrences of 2–12 s of a very high-frequency active phase and then 2–12 s of a completely silent period"

  // when NBQX is applied in GPe, we expect a 51.2% decrease in GPe
  // (cf Kita et al. 2004 : "Role of Ionotropic Glutamatergic and GABAergic Inputs
  // on the Firing Activity of Neurons in the External Pallidum in Awake Monkeys")
  activations.assign(DESACT_NUMBER,1);
  activations[STN_GPe_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_GPi_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_MSN_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_FSI_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_GPe_NMDA] = 0.0f; // to block the NMDA channel from STN
  activations[STN_GPi_NMDA] = 0.0f; // to block the NMDA channel from STN
  activations[STN_MSN_NMDA] = 0.0f; // to block the NMDA channel from STN
  activations[STN_FSI_NMDA] = 0.0f; // to block the NMDA channel from STN
  activations[CMPf_GPe_AMPA] = desactivation_level; // to block the AMPA channel of GPe
  score = _do_trial(means,ref,params,delays,activations,GPe_N,-0.512,0.42,sim_time,mem,false);
  scoreg += score;
  if (verbose) {
    std::cerr << score << " ";
    for (int i=0; i<NUCLEUS_NUMBER; i++) {
      std::cout << means[i*ch_n] << " ";
    }
  }

  // when gabazine is applied in GPe, we expect a 198% increase in GPe
  // (cf Kita et al. 2004 : "Role of Ionotropic Glutamatergic and GABAergic Inputs
  // on the Firing Activity of Neurons in the External Pallidum in Awake Monkeys")
  activations.assign(DESACT_NUMBER,1);
  activations[STN_GPe_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_GPi_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_MSN_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_FSI_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_GPe_NMDA] = 0.0f; // to block the NMDA channel from STN
  activations[STN_GPi_NMDA] = 0.0f; // to block the NMDA channel from STN
  activations[STN_MSN_NMDA] = 0.0f; // to block the NMDA channel from STN
  activations[STN_FSI_NMDA] = 0.0f; // to block the NMDA channel from STN
  activations[MSN_GPe_GABAA] = desactivation_level; // to block the GABA channel of GPe (from D2)
  activations[GPe_GPe_GABAA] = desactivation_level; // to block the GABA channel of GPe (from GPe)
  score = _do_trial(means,ref,params,delays,activations,GPe_N,+1.98,1.38,sim_time,mem,false);
  scoreg += score;
  if (verbose) {
    std::cerr << score << " ";
    for (int i=0; i<NUCLEUS_NUMBER; i++) {
      std::cout << means[i*ch_n] << " ";
    }
  }


  // the data is not supplied with changes, but as an absolute value. The choice is to integrate it with the absolute value
  // when no drug is employed, the firing rate is close to normal (ie near 67.8 +/- 30.2 Hz)
  // (cf Tachibana et al. 2008 : "Motor cortical control of internal pallidal
  // activity through glutamatergic and GABAergic inputs in awake monkeys")
  activations.assign(DESACT_NUMBER,1);
  activations[STN_GPe_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_GPi_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_MSN_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_FSI_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_GPe_NMDA] = 0.0f; // to block the NMDA channel from STN
  activations[STN_GPi_NMDA] = 0.0f; // to block the NMDA channel from STN
  activations[STN_MSN_NMDA] = 0.0f; // to block the NMDA channel from STN
  activations[STN_FSI_NMDA] = 0.0f; // to block the NMDA channel from STN
  ref_tmp[GPi_N*ch_n] = 67.8;
  score = _do_trial(means,ref_tmp,params,delays,activations,GPi_N,0.0f,0.302,sim_time,mem,false);
  scoreg += score;
  if (verbose) {
    std::cerr << score << " ";
    for (int i=0; i<NUCLEUS_NUMBER; i++) {
      std::cout << means[i*ch_n] << " ";
    }
  }

  // when gabazine is applied in GPi, we expect a 28.6% increase in GPi
  // (cf Tachibana et al. 2008 : "Motor cortical control of internal pallidal
  // activity through glutamatergic and GABAergic inputs in awake monkeys")
  activations.assign(DESACT_NUMBER,1);
  activations[STN_GPe_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_GPi_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_MSN_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_FSI_AMPA] = 0.0f; // to block the AMPA channel from STN
  activations[STN_GPe_NMDA] = 0.0f; // to block the NMDA channel from STN
  activations[STN_GPi_NMDA] = 0.0f; // to block the NMDA channel from STN
  activations[STN_MSN_NMDA] = 0.0f; // to block the NMDA channel from STN
  activations[STN_FSI_NMDA] = 0.0f; // to block the NMDA channel from STN
  activations[MSN_GPi_GABAA] = desactivation_level; // to block the GABA channel of GPi (from D1)
  activations[GPe_GPi_GABAA] = desactivation_level; // to block the GABA channel of GPi (from GPe)
  score = _do_trial(means,ref,params,delays,activations,GPi_N,+0.286,0.128,sim_time,mem,false);
  scoreg += score;
  if (verbose) {
    std::cerr << score << " ";
    for (int i=0; i<NUCLEUS_NUMBER; i++) {
      std::cout << means[i*ch_n] << " ";
    }
  }
  // NB "very regular and oscillatory firing"
  
  // when muscimol is applied in GPe, shutting it, we expect a 38.3% increase in GPi
  // (cf Tachibana et al. 2008 : "Motor cortical control of internal pallidal
  // activity through glutamatergic and GABAergic inputs in awake monkeys")
  activations.assign(DESACT_NUMBER,1);
  activations[GPe_MSN_GABAA] = 0.0f; // to block the GABAA channel from GPe
  activations[GPe_FSI_GABAA] = 0.0f; // to block the GABAA channel from GPe
  activations[GPe_STN_GABAA] = 0.0f; // to block the GABAA channel from GPe
  activations[GPe_GPi_GABAA] = 0.0f; // to block the GABAA channel from GPe
  score = _do_trial(means,ref,params,delays,activations,GPi_N,+0.383,0.331,sim_time,mem,false);
  scoreg += score;
  if (verbose) {
    std::cerr << score << " ";
    for (int i=0; i<NUCLEUS_NUMBER; i++) {
      std::cout << means[i*ch_n] << " ";
    }
  }

  // when gabazine is applied in GPe, we expect a 30.3% decrease in GPi
  // (cf Tachibana et al. 2008 : "Motor cortical control of internal pallidal
  // activity through glutamatergic and GABAergic inputs in awake monkeys")
  activations.assign(DESACT_NUMBER,1);
  activations[MSN_GPe_GABAA] = desactivation_level; // to block the GABAA channel of GPe (from MSN)
  activations[GPe_GPe_GABAA] = desactivation_level; // to block the GABAA channel of GPe (from GPe)
  score = _do_trial(means,ref,params,delays,activations,GPi_N,-0.303,0.229,sim_time,mem,false);
  scoreg += score;
  if (verbose) {
    std::cerr << score << " ";
    for (int i=0; i<NUCLEUS_NUMBER; i++) {
      std::cout << means[i*ch_n] << " ";
    }
  }

  return scoreg;
}


float calc_score_selective_boutons(std::vector <float>& params, bool verbose, int studied_nucleus)
{
  // computes part of the anatomical plausibility objective
  float weak[] = {15, 149}; 
  float high[] = {150, 1499};
  float massive[] = {1500, 4999};
  float ctxpt_stn[] = {0, 1499};
  float not_nonexistent[] = {15, 4999};
  float ctx_msn_range[] = {250, 4999};
  float ctxpt_msn_fsi_range[] = {100, 2499};
  std::vector <int> c_weak; c_weak.resize(0);
  std::vector <int> c_high; c_high.resize(0);
  std::vector <int> c_massive; c_massive.resize(0);
  std::vector <int> c_ctxpt_stn; c_ctxpt_stn.resize(0);
  std::vector <int> c_not_nonexistent; c_not_nonexistent.resize(0);
  std::vector <int> c_ctx_msn_range; c_ctx_msn_range.resize(0);
  std::vector <int> c_ctxpt_msn_fsi_range; c_ctxpt_msn_fsi_range.resize(0);

  if (studied_nucleus == MSN_N || studied_nucleus == -1) {
    //axons
    ///CTX -> MSN+FSI treated separately after this function (for now)
    c_weak.push_back(STN_MSN);
    //c_not_nonexistent.push_back(CTX_MSN);
    //c_not_nonexistent.push_back(GPe_MSN); // data exists only for the rat (Bevan98) 
    c_massive.push_back(FSI_MSN); // estimation from Humphries10b is ~3000
    c_high.push_back(MSN_MSN); // estimation from Humphries10b is ~700; Estimation from Wilson07 & Wickens07 is ~450; Wickens07 wth Lee 1997, estimates at ~550
    c_massive.push_back(CMPf_MSN); // estimation from Parent05
    c_ctx_msn_range.push_back(CTX_MSN); //TODO
  }

  if (studied_nucleus == FSI_N || studied_nucleus == -1) {
    //axons
    //c_not_nonexistent.push_back(CTX_FSI);
    c_weak.push_back(STN_FSI);
    //c_not_nonexistent.push_back(GPe_FSI); // data exists only for the rat (Bevan98)
    c_weak.push_back(FSI_FSI); // estimation from Humphries10b is < 100
    c_high.push_back(CMPf_FSI); // estimation from Sidibe99: 1/3 of the contacts happen on FSI
    c_ctxpt_msn_fsi_range.push_back(CTX_FSI); //TODO
  }

  if (studied_nucleus == STN_N || studied_nucleus == -1) {
    //axons
    ///CTX -> STN treated separately after this function (for now)
    //c_not_nonexistent.push_back(CTX_STN);
    c_high.push_back(GPe_STN);
    c_weak.push_back(CMPf_STN); // estimation from Sadikot92b and Tande06: weak
    c_ctxpt_msn_fsi_range.push_back(CTX_STN); //TODO
  }

  if (studied_nucleus == GPe_N || studied_nucleus == -1) {
    //axons
    c_high.push_back(MSN_GPe); // Wickens07: ~749
    c_high.push_back(STN_GPe);
    c_high.push_back(GPe_GPe); // (hazardous ?) estimate: the Shink96 count of synapses (2.5 times more synapses from STN than from GPe) and we can compare the number of neurons (~3 more neurons in the GPe than in the STN), leading to the estimate that their should be slighty more boutons in the GPe -> GPe connection than in the STN -> GPe connection. Hence the value is expected to be around 500-1000, in the "high" range as the STN -> GPe connection.
    c_weak.push_back(CMPf_GPe); // estimation from Sadikot92b, Parent05 and Tande06: weak
  }

  if (studied_nucleus == GPi_N || studied_nucleus == -1) {
    //axons
    c_high.push_back(MSN_GPi); // Wickens07: ~749
    c_high.push_back(STN_GPi);
    c_high.push_back(GPe_GPi);
    c_weak.push_back(CMPf_GPi); // estimation from Parent05 and Tande06: weak
  }

  float borneinf = 0.;
  float bornesup = 5000.;

  float score = 0.;
  float score_c = 1.; float score_d = 1.;

  float *mini = NULL; float *maxi = NULL;

  int nb_elem_treated = 0;
  bool elem_already_treated = false;

  for (int i=0; i<PARAMS_NUMBER; i++) {
    if (is_in(c_weak,i)) {
      mini = weak; maxi = weak+1;
    } else if (is_in(c_high,i)) {
      mini = high; maxi = high+1;
    } else if (is_in(c_massive,i)) {
      mini = massive; maxi = massive+1;
    } else if (is_in(c_ctxpt_stn,i)) {
      mini = ctxpt_stn; maxi = ctxpt_stn+1;
    } else if (is_in(c_not_nonexistent,i)) {
      mini = not_nonexistent; maxi = not_nonexistent+1;
    } else if (is_in(c_ctx_msn_range,i)) {
      mini = ctx_msn_range; maxi = ctx_msn_range+1;
    } else if (is_in(c_ctxpt_msn_fsi_range,i)) {
      mini = ctxpt_msn_fsi_range; maxi = ctxpt_msn_fsi_range+1;
    } else {
      mini = NULL;
    }
    if (mini) {
      nb_elem_treated++;
      elem_already_treated = true;
      score_c = _is_near_tronqgaussian((*mini+*maxi)/2,params[i],((*mini+*maxi)/2)- *mini);
    } else {
      score_c = 0.;
    }
    if (verbose) {
      if (score_c < 1 && mini) {
        std::cout << "!!!! boutons count n° " << i << " out of range (" << params[i] << " not in [" << *mini << "," << *maxi << "])" << std::endl;
      } else if (mini) {
        std::cout << "     boutons count n° " << i << " inside range (" << params[i] << " in [" << *mini << "," << *maxi << "])" << std::endl;
      }
    }
    if (studied_nucleus == -1 || elem_already_treated) { 
      score += score_c;
    }
    elem_already_treated = false;
  }

  //	c_prox.push_back(DIST_GPe_MSN); //no data for monkey (for the rat, it would have been "middle" according to Bevan98)
  //	c_prox.push_back(DIST_GPe_FSI); //no data for monkey (for the rat, it would have been "middle" according to Bevan98)
  if (studied_nucleus == -1) {
    nb_elem_treated++;
#ifdef LINEARDIST
    score_c = _linear_dist(params[GPe_MSN]+params[GPe_FSI],high[0],high[1],borneinf,bornesup);
#else
    score_c = _is_near_tronqgaussian((high[0]+high[1])/2,params[GPe_MSN]+params[GPe_FSI],((high[0]+high[1])/2)- high[0]);
#endif
    if (score_c < 1 && verbose) {
      std::cout << "pallido-striatal inputs out of range (score = " << score_c << ")" << std::endl;
    }
    score += score_c;
  }

  if (verbose) {
    std::cout << "Summary: " << score << "/" << nb_elem_treated << std::endl;
  }


  return score;
}

float calc_score_selective_axons(std::vector <float>& params, bool verbose, int studied_nucleus)
{
  // computes part of the anatomical plausibility objective
  float weak[] = {15, 149}; // this constant and following ranges are for the bouton counts
  float high[] = {150, 749};
  float pretty_massive[] = {500, 1999};
  float massive[] = {1500, 4999};
  float not_nonexistent[] = {15, 749};
  float ctx_msn_range[] = {250, 4999};
  float ctxpt_msn_fsi_range[] = {1, 999};
  float ctxpt_stn[] = {25, 4999};
  std::vector <int> c_weak; c_weak.resize(0);
  std::vector <int> c_high; c_high.resize(0);
  std::vector <int> c_massive; c_massive.resize(0);
  std::vector <int> c_pretty_massive; c_pretty_massive.resize(0);
  std::vector <int> c_not_nonexistent; c_not_nonexistent.resize(0);
  std::vector <int> c_ctx_msn_range; c_ctx_msn_range.resize(0);
  std::vector <int> c_ctxpt_msn_fsi_range; c_ctxpt_msn_fsi_range.resize(0);
  std::vector <int> c_ctxpt_stn; c_ctxpt_stn.resize(0);

  // range for the synaptic localization
  float prox[] = {0., 0.2}; 
  float middle[] = {0.2, 0.6};
  float far[] = {0.6, 1.};
  std::vector <int> c_prox; c_prox.resize(0);
  std::vector <int> c_middle; c_middle.resize(0);
  std::vector <int> c_far; c_far.resize(0);

  if (studied_nucleus == MSN_N || studied_nucleus == -1) {
    //axons
    ///CTX -> MSN+FSI treated separately after this function (for now)
//    c_weak.push_back(STN_MSN); // see below
    //c_not_nonexistent.push_back(CTX_MSN);
    //c_not_nonexistent.push_back(GPe_MSN); // data exists only for the rat (Bevan98) and is difficult to interpret (Mallet12)
    c_massive.push_back(FSI_MSN); // estimation from Humphries10b is ~3000
    c_high.push_back(MSN_MSN); // estimation from Humphries10b is ~700; Estimation from Wilson07 & Wickens07 is ~450; Wickens07 wth Lee 1997, estimates at ~550
    c_massive.push_back(CMPf_MSN); // estimation from Parent05 (~1500 to ~3500 varicosities) and Sidibe99 (2/3 of these should target MSN neurons)
    c_ctx_msn_range.push_back(CTX_MSN); //TODO
    c_ctxpt_msn_fsi_range.push_back(CTXPT_MSN); //TODO

    // synapses
    c_far.push_back(DIST_CTX_MSN); // [rat] Wilson07 : "The synapses formed by cortical axons end almost exclusively (about 95%) on dendritic spines (Kemp and Powell, 1971b; Somogyi et al., 1981; Xu et al., 1989)"
    c_far.push_back(DIST_CTXPT_MSN); // [rat] Wilson07 : "The synapses formed by cortical axons end almost exclusively (about 95%) on dendritic spines (Kemp and Powell, 1971b; Somogyi et al., 1981; Xu et al., 1989)"
    //	c_prox.push_back(DIST_STN_MSN); // not known
    c_prox.push_back(DIST_FSI_MSN); // [rat] Tepper08 (p. 274); Wilson07 (p. 95)
    c_far.push_back(DIST_MSN_MSN); // [rat] Tepper08 (p. 276); Wilson07 (p. 95)
    c_middle.push_back(DIST_CMPf_MSN); // Sidibe96 + Sidibe99
  }

  if (studied_nucleus == FSI_N || studied_nucleus == -1) {
    //axons
    //c_not_nonexistent.push_back(CTX_FSI);
//    c_weak.push_back(STN_FSI); // see below
    //c_not_nonexistent.push_back(GPe_FSI); // data exists only for the rat (Bevan98)
    c_weak.push_back(FSI_FSI); // estimation from Humphries10b is < 100
    c_pretty_massive.push_back(CMPf_FSI); // estimation from Sidibe99: 1/3 of the contacts happen on FSI, and Parent05 : (~1500 to ~3500 varicosities for these axons)
    c_ctx_msn_range.push_back(CTX_FSI); //TODO
    c_ctxpt_msn_fsi_range.push_back(CTXPT_FSI); //TODO

    // synapses
    c_far.push_back(DIST_CTX_FSI); // Lapper92 (table 1, p. 219)
    c_far.push_back(DIST_CTXPT_FSI); // Lapper92 (table 1, p. 219)
    //	c_prox.push_back(DIST_STN_FSI); // not known
    //	c_prox.push_back(DIST_FSI_FSI);
    c_prox.push_back(DIST_CMPf_FSI); // Sidibe99
  }

  if (studied_nucleus == STN_N || studied_nucleus == -1) {
    //axons
    ///CTX -> STN treated separately after this function (for now)
    //c_not_nonexistent.push_back(CTX_STN);

    //c_high.push_back(GPe_STN); // Following Martin Parent advice
    c_weak.push_back(GPe_STN); // cf Karachi 05
    //c_not_nonexistent.push_back(GPe_STN); // we can't really know better

    c_weak.push_back(CMPf_STN); // estimation from Sadikot92b and Tande06: weak
    c_ctxpt_stn.push_back(CTX_STN); //TODO

    // synapses
    c_far.push_back(DIST_CTX_STN); // Marani 08 (ch 2)
    // c_prox.push_back(DIST_GPe_STN); // cf Shink 96 : 31% axosomatic and 69% axodendritic; and cf Parent 95b : 30% soma, 40% proximal dendrites and 30% distant dendrite. Hassler82 + Marani08 => labelled it as "prox"
    c_middle.push_back(DIST_GPe_STN); // there was an ambiguity, as this is the term always got wrong after evolution, we labelled it as "middle". Coherent with the rat (Bevan97) which show a repartion : 31% on perikarya, 39% on large dendrites, 30% on small dendrites
  }

  if (studied_nucleus == GPe_N || studied_nucleus == -1) {
    //axons
    c_high.push_back(MSN_GPe); // Wickens07: ~749
    c_high.push_back(STN_GPe);
    //c_high.push_back(GPe_GPe); // (hazardous ?) estimate: the Shink96 count of synapses (2.5 times more synapses from STN than from GPe) and we can compare the number of neurons (~3 more neurons in the GPe than in the STN), leading to the estimate that their should be slighty more boutons in the GPe -> GPe connection than in the STN -> GPe connection. Hence the value is expected to be around 500-1000, in the "high" range as the STN -> GPe connection.
#ifdef TRYGPEGPE
    c_high.push_back(GPe_GPe);
#else
    c_weak.push_back(GPe_GPe); // Following Martin Parent advice
#endif
    c_weak.push_back(CMPf_GPe); // estimation from Sadikot92b, Parent05 and Tande06: weak

    // synapses
    c_middle.push_back(DIST_MSN_GPe); // cf Shink 95: but we believe it should be more 0.2 than 0.6 here !
    c_middle.push_back(DIST_STN_GPe); // cf Shink 95 and Shink 96
    c_prox.push_back(DIST_GPe_GPe); // cf Shink 95
  }

  if (studied_nucleus == GPi_N || studied_nucleus == -1) {
    //axons
    c_high.push_back(MSN_GPi); // Wickens07: ~749
    c_high.push_back(STN_GPi);
#ifdef TRYGPEGPI
    c_high.push_back(GPe_GPi); // following Martin Parent advice
#else
    c_not_nonexistent.push_back(GPe_GPi); // we can't really know better
#endif
    c_weak.push_back(CMPf_GPi); // estimation from Parent05 and Tande06: weak

    // synapses
    c_middle.push_back(DIST_MSN_GPi); // cf Shink 95: as well axosynaptic as axodendritic
    c_middle.push_back(DIST_STN_GPi); // cf Shink 95
    c_prox.push_back(DIST_GPe_GPi); // cf Shink 95. More or less coherent with the rat (Bevan97)
  }

  float borneinf = 0.;
  float bornesup = 5000.;

  float score = 0.;
  float score_c; float score_d;

  float *mini = NULL; float *maxi = NULL;

  int nb_elem_treated = 0;
  bool elem_already_treated = false;

  for (int i=0; i<CONNECT_NUMBER; i++) {
    if (is_in(c_weak,i)) {
      mini = weak; maxi = weak+1;
    } else if (is_in(c_high,i)) {
      mini = high; maxi = high+1;
    } else if (is_in(c_massive,i)) {
      mini = massive; maxi = massive+1;
    } else if (is_in(c_pretty_massive,i)) {
      mini = pretty_massive; maxi = pretty_massive+1;
    } else if (is_in(c_ctxpt_stn,i)) {
      mini = ctxpt_stn; maxi = ctxpt_stn+1;
    } else if (is_in(c_not_nonexistent,i)) {
      mini = not_nonexistent; maxi = not_nonexistent+1;
    } else if (is_in(c_ctx_msn_range,i)) {
      mini = ctx_msn_range; maxi = ctx_msn_range+1;
    } else if (is_in(c_ctxpt_msn_fsi_range,i)) {
      mini = ctxpt_msn_fsi_range; maxi = ctxpt_msn_fsi_range+1;
    } else {
      mini = NULL;
    }
    if (mini) {
      nb_elem_treated++;
      elem_already_treated = true;
      score_c = _is_near_tronqgaussian((*mini+*maxi)/2,params[i],((*mini+*maxi)/2)- *mini);
    } else {
      score_c = 1.;
    }
    if (verbose) {
      if (score_c < 1) {
        std::cout << "!!!! boutons count n° " << i << " out of range (" << params[i] << " not in [" << *mini << "," << *maxi << "])" << std::endl;
      } else if (mini) {
        std::cout << "     boutons count n° " << i << " inside range (" << params[i] << " in [" << *mini << "," << *maxi << "])" << std::endl;
      }
    }
    if (is_in(c_prox,i+CONNECT_NUMBER)) {
      mini = prox; maxi = prox+1;
    } else if (is_in(c_middle,i+CONNECT_NUMBER)) {
      mini = middle; maxi = middle+1;
    } else if (is_in(c_far,i+CONNECT_NUMBER)) {
      mini = far; maxi = far+1;
    } else {
      mini = NULL;
    }
    if (mini) {
      elem_already_treated = true;
      if (!elem_already_treated) {
        nb_elem_treated++;
      }
      score_d = _is_near_tronqgaussian((*mini+*maxi)/2,params[i+CONNECT_NUMBER],((*mini+*maxi)/2)- *mini);
    } else {
      score_d = 1.;
    }
    if (verbose) {
      if (score_d < 1) {
        std::cout << "!!!! distance n° " << i+CONNECT_NUMBER << " out of range (" << params[i+CONNECT_NUMBER] << " not in [" << *mini << "," << *maxi << "])" << std::endl;
      } else if (mini) {
        std::cout << "     distance n° " << i+CONNECT_NUMBER << " inside range (" << params[i+CONNECT_NUMBER] << " in [" << *mini << "," << *maxi << "])" << std::endl;
      }
    }
    if (studied_nucleus == -1 || elem_already_treated) { 
      score += score_c * score_d;
    }
    elem_already_treated = false;
  }

  //	c_prox.push_back(DIST_GPe_MSN); //no data for monkey (for the rat, it would have been "middle" according to Bevan98)
  //	c_prox.push_back(DIST_GPe_FSI); //no data for monkey (for the rat, it would have been "middle" according to Bevan98)
  if (studied_nucleus == -1) {
    nb_elem_treated++;
    score_c = _is_near_tronqgaussian((high[0]+high[1])/2,params[GPe_MSN]+params[GPe_FSI],((high[0]+high[1])/2)- high[0]);
    if (score_c < 1 && verbose) {
      std::cout << "pallido-striatal inputs out of range (score = " << score_c << ")" << std::endl;
    }
    score += score_c;
  }

  if (studied_nucleus == -1) {
    nb_elem_treated++;
    score_c = _is_near_tronqgaussian((weak[0]+weak[1])/2,params[STN_MSN]+params[STN_FSI],((weak[0]+weak[1])/2)- weak[0]);
    if (score_c < 1 && verbose) {
      std::cout << "subthalamico-striatal inputs out of range (score = " << score_c << ")" << std::endl;
    }
    score += score_c;
  }

  if (verbose) {
    std::cout << "Summary: " << score << "/" << nb_elem_treated << std::endl;
  }


  return score;
}

float _is_near_gaussian(
    float reference,
    float mean,
    float absolute_radius)
{
  return exp((-0.5)*(mean-reference)*(mean-reference)/(absolute_radius*absolute_radius));
}

float _is_near_tronqgaussian(
    float reference,
    float mean,
    float absolute_radius)
{
  if ((mean >= reference - absolute_radius) and (mean <= reference + absolute_radius)) {
    return 1.;
  }
  return (exp((-0.5)*(mean-reference)*(mean-reference)/(absolute_radius*absolute_radius)))/(exp(-0.5));
}

float _has_changed_near_gaussian(
    float reference,
    float mean,
    float proportional_change,
    float proportional_radius)
{
  float new_ref = reference + reference * proportional_change;
  float absolute_radius = abs(new_ref - (reference + reference * (proportional_change+proportional_radius)));
  if (absolute_radius < 1) {
    return 0.;
  } else {
    return exp((-0.5)*(mean-new_ref)*(mean-new_ref)/(absolute_radius*absolute_radius));
  }
}

float _has_changed_near_tronqgaussian(
    float reference,
    float mean,
    float proportional_change,
    float proportional_radius)
{
  float new_ref = reference + reference * proportional_change;
  float absolute_radius = abs(new_ref - (reference + reference * (proportional_change+proportional_radius)));
  if ((mean >= new_ref - absolute_radius) and (mean <= new_ref + absolute_radius)) {
    return 1.;
  }
  return (exp((-0.5)*(mean-new_ref)*(mean-new_ref)/(absolute_radius*absolute_radius)))/(exp(-0.5));
}

float _is_near_step(
    float reference,
    float mean,
    float absolute_radius)
{
  if ((mean >= reference - absolute_radius) and (mean <= reference + absolute_radius)) {
    return 1;
  }
  return 0;
}



