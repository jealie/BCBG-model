#include "bcbg2.hpp"

float BCBG2::basic_test(int ch, float time)
{
  float debut_time = 450.;
  float unity = 0.25;
  float base_rate = 2.;
  float mid_rate = 3.;
  float max_rate = 4.;
  if (time < unity + debut_time) {
    return base_rate;
  } else if (time < 2*unity + debut_time) {
    if (ch==0) {
      return mid_rate;
    } else {
      return base_rate;
    }
  } else if (time < 3*unity + debut_time) {
    if (ch==0) {
      return mid_rate;
    } else if (ch==1) {
      return max_rate;
    } else {
      return base_rate;
    }
  } else if (time < 4*unity + debut_time) {
    if (ch==0||ch==1) {
    //if (ch==0) {
      return max_rate;
    } else {
      return base_rate;
    }
  } else if (time < 5*unity + debut_time) {
    if (ch==0) {
      return mid_rate;
    } else if (ch==1) {
    //} else if (ch==10) {
      return max_rate;
    } else {
      return base_rate;
    }
  } else {
    return base_rate;
  }
}

float BCBG2::exploexplo_testbis(int ch, float time, float debut_time, float end_time)
{
  // this function implements some selection tests for the model
  float base_rate = 0.;
  float mid_rate = 10.;
  float max_rate = 20.0;
  float transition_period = 0.001;
  if (time < debut_time || time >= end_time) {
    return base_rate;
  } else if (time >= end_time-transition_period) {
    if (ch<=0) {
      return base_rate;
    } else if (ch==1) {
      return std::max((float)(mid_rate*(0.5-0.5*((float)sin(3.1415*((time-(end_time-transition_period))/transition_period-0.5))))),base_rate);
    } else {
      return std::max(max_rate*(1-(time-(end_time-transition_period))/transition_period),base_rate);
    }
  } else if (time < debut_time+transition_period) {
    if (ch<=0) {
      return base_rate;
    } else if (ch==1) {
      return std::max((float)(mid_rate*(0.5+0.5*((float)sin(3.1415*((time-debut_time)/transition_period-0.5))))),base_rate);
    } else {
      return std::max((float)(max_rate*(0.5+0.5*((float)sin(3.1415*((time-debut_time)/transition_period-0.5))))),base_rate);
    }
  } else {
    if (ch<=0) {
      return base_rate;
    } else if (ch==1) {
      return mid_rate;
    } else {
      return max_rate;
    }
  }
}

float BCBG2::exploexplo_test(int ch, float time, int i)
{
  // simple input scheme used in the chapter 7 of the thesis
  float debut_time = 5.;
  float base_rate = 2.;
  float mid_rate = 2. + (9.0 )/((float) i);
  float max_rate = 2. + (18.0)/((float) i);
  if (time < debut_time) {
    return base_rate;
  } else {
    if (ch==0) {
      return base_rate;
    } else if (ch==1) {
      return mid_rate;
    } else {
      return max_rate;
    }
  }
}

float BCBG2::visualization_humphries_et_al_2006(int ch, float time)
{
  // mimics the selection task of Humphries et al, 2006 (see function below)
  float debut_time = 40.;
  float base_rate = 0.;
  float mid_rate = 20.;
  float max_rate = 40.;
  if (time < 1 + debut_time) {
    return base_rate;
  } else if (time < 2.5 + debut_time) {
    if (ch==0) {
      return mid_rate;
    } else {
      return base_rate;
    }
  } else {
    if (ch==0) {
      return mid_rate;
    } else if (ch==1) {
      return max_rate;
    } else {
      return base_rate;
    }
  }
}

void BCBG2::updateMultiChannelsNucleus_humphries_test(int steps, float time)
{
  // mimics the selection task of Humphries et al, 2006
  int i,s,ch_i;
  float rnd_nb;
  int S_channels_nb = MCN[0]->get_channels_number();
  // multi_channels_nucleus n°0 is always the input - we handle them differently
  for (s=0;s<steps;s++) {
    for (ch_i=0; ch_i<S_channels_nb; ch_i++) {
      rnd_nb=visualization_humphries_et_al_2006(ch_i, time);
      MCN[0]->set_S(rnd_nb,previous_tmod,ch_i);
    }
    for (i=1; i < n_nuclei; i++) {
      MCN[i]->update_multi_channels_nucleus_evo();
    }
    previous_tmod = tmod;
    tmod++;
    time+= dt;
    if (tmod == max_tau) {
      tmod = 0;
    }
  }
}

float BCBG2::scalevar(float m) {
  float r = BCBG2::_uni_one_gaussian()+m-1.0f;
  while (r < 0) {
    r = _uni_one_gaussian()+m-1.0f;
  }
  return r;
}

float BCBG2::scalevarreduc(float m, float reduc) {
  float r = BCBG2::_uni_one_gaussian()/reduc+m-1.0f;
  while (r < 0) {
    r = _uni_one_gaussian()/reduc+m-1.0f;
  }
  return r;
}

void BCBG2::updateMultiChannelsNucleus_exploexplo_test(int steps, float time, int input_f)
{
  // sinusoidal selection task (featured in Lienard and Girard 2013)
  int i,s,ch_i;
  float rnd_nb_ctx=2.;
  float rnd_nb_cmpf=4.;
#ifdef ITPTCTX
  float rnd_nb_ctxPT=15.;
#endif
  int S_channels_nb = MCN[0]->get_channels_number();
  for (s=0;s<steps;s++) {
    for (ch_i=0; ch_i<S_channels_nb; ch_i++) {
      if (input_f == 0) {
        // check only with noise
        MCN[0]->set_S(scalevar(rnd_nb_ctx),previous_tmod,ch_i);
        MCN[1]->set_S(scalevar(rnd_nb_cmpf),previous_tmod,ch_i);
#ifdef ITPTCTX
        MCN[CTXPT_N]->set_S(scalevar(rnd_nb_ctxPT),previous_tmod,ch_i);
#endif
      } else {
        // do the real selection task
        float ch_nf = S_channels_nb;
        float ch_if = (2.0 * (1.0+(int)((ch_i-1)/2)));
        float ch_if_alt = (2.0 * (1.0+(int)((ch_nf-1-(ch_i-1))/2)));
        MCN[0]->set_S(scalevar(rnd_nb_ctx)+0.2*exploexplo_testbis(1, time, 6.,11)
            *(0.5-0.5*((float)sin(3.1415*((1.0-(ch_nf/ch_if))/(ch_nf/ch_if)-0.5)))) 
            ,previous_tmod,ch_i); ///////// revised version of sinusoidal input

        //sustained input to the CMPf and CTXPT:
        MCN[CMPf_N]->set_S(scalevar(rnd_nb_cmpf),tmod,ch_i);
#ifdef ITPTCTX
        MCN[CTXPT_N]->set_S(scalevar(rnd_nb_ctxPT),tmod,ch_i);
#endif
      }
    }

#ifdef ITPTCTX
    for (i=2; i < n_nuclei-1; i++) 
#else
    for (i=2; i < n_nuclei; i++) 
#endif
    {
        MCN[i]->update_multi_channels_nucleus_evo();
    }
    previous_tmod = tmod;
    tmod++;
    time+= dt;
    if (tmod == max_tau) {
      tmod = 0;
    }
  }
}

void BCBG2::updateMultiChannelsNucleus_basic_test(int steps, float time)
{
  int i,s,ch_i;
  float rnd_nb;
  int S_channels_nb = MCN[0]->get_channels_number();
  // multi_channels_nucleus n°0 is always the input - we handle them differently
  for (s=0;s<steps;s++) {
    for (ch_i=0; ch_i<S_channels_nb; ch_i++) {
      rnd_nb=basic_test(ch_i, time);
      MCN[0]->set_S(rnd_nb,previous_tmod,ch_i);
    }
    MCN[1]->set_S(4.,previous_tmod,ch_i);
    for (i=2; i < n_nuclei; i++) {
      MCN[i]->update_multi_channels_nucleus_evo();
    }
    previous_tmod = tmod;
    tmod++;
    time+= dt;
    if (tmod == max_tau) {
      tmod = 0;
    }
  }
}

void BCBG2::updateMultiChannelsNucleus(int steps)
{
  int i,s,ch_i;
  float rnd_nb;
#ifdef ITPTCTX
  float rnd_nb_ctxPT=15.;
#endif
  int S_channels_nb = MCN[0]->get_channels_number();
  // multi_channels_nucleus n°0 is always the input - we handle them differently
  for (s=0;s<steps;s++) {
    for (ch_i=0; ch_i<S_channels_nb; ch_i++) {
      rnd_nb=_uni_one_gaussian()*2.0f;
      while (rnd_nb < 0) {
        rnd_nb=_uni_one_gaussian()*2.0f;
      }
      MCN[0]->set_S(2.,previous_tmod,ch_i);
      MCN[1]->set_S(4.,previous_tmod,ch_i);
#ifdef ITPTCTX
      MCN[CTXPT_N]->set_S(rnd_nb_ctxPT,previous_tmod,ch_i);
#endif
    }

#ifdef ITPTCTX
    for (i=2; i < n_nuclei-1; i++) 
#else
    for (i=2; i < n_nuclei; i++) 
#endif

    {
      MCN[i]->update_multi_channels_nucleus_evo();
    }
    previous_tmod = tmod;
    tmod++;
    if (tmod == max_tau) {
      tmod = 0;
    }
  }
}

void BCBG2::updateSingleChannelNucleusWithNoise(int steps, int integration_method)
{
  int i,s;
  float rnd_nb_ctx,rnd_nb_th;
  float rnd_nb_ctxPT;
  // multi_channels_nucleus n°0 is always the input - we handle them differently
  for (s=0;s<steps;s++) {
    rnd_nb_ctx=scalevar(2.0);
    SCN[CTX_N]->set_S(rnd_nb_ctx,previous_tmod);
    rnd_nb_th=4.;
    SCN[CMPf_N]->set_S(scalevar(rnd_nb_th),previous_tmod);
#ifdef ITPTCTX
    rnd_nb_ctxPT=15.;
    SCN[CTXPT_N]->set_S(rnd_nb_ctx+13.0,previous_tmod);
    for (i=2; i < n_nuclei-1; i++) 
#else
    for (i=2; i < n_nuclei; i++) 
#endif
    {
      if (integration_method == 0) {
        SCN[i]->update_single_channel_nucleus_euler();
      } else {
        SCN[i]->update_single_channel_nucleus_rk3();
      }
    }
    previous_tmod = tmod;
    tmod++;
    if (tmod == max_tau) {
      tmod = 0;
    }
  }
}

void BCBG2::updateSingleChannelNucleus(int steps, int integration_method)
{
  // this is the function that computes the next state for all nuclei (in the single channel case)
  // two integration methods are exposed here: Runge-Kutta (third order) or simple Euler method
  int i,s;
  float rnd_nb_ctx,rnd_nb_th;
#ifdef ITPTCTX
  float rnd_nb_ctxPT;
#endif
  // multi_channels_nucleus n°0 is always the input - we handle them differently
  for (s=0;s<steps;s++) {

    rnd_nb_ctx=2.;
    SCN[CTX_N]->set_S(rnd_nb_ctx,previous_tmod);

    rnd_nb_th=4.;
    SCN[CMPf_N]->set_S(rnd_nb_th,previous_tmod);

#ifdef ITPTCTX
    rnd_nb_ctxPT=15.;
    SCN[CTXPT_N]->set_S(rnd_nb_ctxPT,previous_tmod);
    for (i=2; i < n_nuclei-1; i++) 
#else
    for (i=2; i < n_nuclei; i++) 
#endif
    {
      if (integration_method == 0) {
        SCN[i]->update_single_channel_nucleus_euler();
      } else {
        SCN[i]->update_single_channel_nucleus_rk3();
      }
    }
    previous_tmod = tmod;
    tmod++;
    if (tmod == max_tau) {
      tmod = 0;
    }
  }
}

bool BCBG2::testTimingResponse()
{
  // (not currently used, this was defined originally as a test for the timings within the BG)
  int i;
  float time;
  float cortical_input;
  bool msn_passed = false;
  bool stn_passed = false;
  bool gpe_passed = false;
  bool gpi_passed = false;

  std::vector <float> refs;
  refs.assign(NUCLEUS_NUMBER,0.);
  for (i=2; i < n_nuclei; i++) {
    refs[i] = SCN[i]->get_S();
  }
  
  for (time=0.;time<1.;time+=dt) {

    // initialize the inputs
    if (time < 0.3) {
      cortical_input=40.;
    } else {
      cortical_input=2.;
    }

    // check whether MSN responded on time
    if (msn_passed == false && time > 0.75) {
      if (abs(SCN[MSN_N]->get_S() - refs[MSN_N]) > 0.25*4.) {
        msn_passed = true;
      }
    }

    // check whether STN responded on time
    if (stn_passed == false && time > 0.60) {
      if (abs(SCN[STN_N]->get_S() - refs[STN_N]) > 1.5*4.) {
        stn_passed = true;
      }
    }

    // check whether GPe responded on time
    if (gpe_passed == false && time > 0.75) {
      if (abs(SCN[GPe_N]->get_S() - refs[GPe_N]) > 2.5*4.) {
        gpe_passed = true;
      }
    }

    // check whether GPi responded on time
    if (gpi_passed == false && time > 0.75) {
      if (abs(SCN[GPi_N]->get_S() - refs[GPi_N]) > 2.5*4.) {
        gpi_passed = true;
      }
    }

    SCN[0]->set_S(cortical_input,previous_tmod);
    SCN[1]->set_S(4.,previous_tmod);
    for (i=2; i < n_nuclei; i++) {
      SCN[i]->update_single_channel_nucleus_evo();
    }
    previous_tmod = tmod;
    tmod++;
    if (tmod == max_tau) {
      tmod = 0;
    }
  }
  return (msn_passed && stn_passed && gpe_passed && gpi_passed);
}

void BCBG2::save_all()
{
  // this stores the state of all variables for the single channel nuclei
  int i;
  for (i=2; i < n_nuclei; i++) {
      if (nucleus_type == 0) {
        SCN[i]->save();
      } else {
        std::cout << "Saving memory in not implemented for this type of nucleus" << std::endl;
      }
  }
  tmod_backup = tmod;
  previous_tmod_backup = previous_tmod;
}

void BCBG2::save_all(MemoryBCBG2& mem)
{
  // more advanced version, which handles multi- and single-channels, as well as an externally defined structure (`mem`)
  int i;
  if (nucleus_type == 0) {
    mem.scns.resize(n_nuclei);
  } else if (nucleus_type == 1) {
    mem.mcns.resize(n_nuclei);
  }
  for (i=2; i < n_nuclei; i++) {
      if (nucleus_type == 0) {
        SCN[i]->save(mem.scns[i]);
      } else if (nucleus_type == 1) {
        MCN[i]->save(mem.mcns[i]);
      } else {
        std::cout << "Saving memory in not implemented for this type of nucleus" << std::endl;
      }
  }
  mem.tmod_backup = tmod;
  mem.previous_tmod_backup = previous_tmod;
}

void BCBG2::load_all()
{
  // loading counterpart of the above functions
  int i;
  for (i=2; i < n_nuclei; i++) {
      if (nucleus_type == 0) {
        SCN[i]->load();
      } else {
        std::cout << "Saving memory in not implemented for this type of nucleus" << std::endl;
      }
  }
  tmod = tmod_backup;
  previous_tmod = previous_tmod_backup;
}

void BCBG2::load_all(MemoryBCBG2& mem)
{
  // loading counterpart of the above functions
  int i;
  for (i=2; i < n_nuclei; i++) {
      if (nucleus_type == 0) {
        SCN[i]->load(mem.scns[i]);
      } else if (nucleus_type == 1) {
        MCN[i]->load(mem.mcns[i]);
      } else {
        std::cout << "Saving memory in not implemented for this type of nucleus" << std::endl;
      }
  }
  tmod = mem.tmod_backup;
  previous_tmod = mem.previous_tmod_backup;
}

void BCBG2::stabilize_all(int steps)
{
  int i,s;
    for (i=2; i < n_nuclei; i++) {
      if (nucleus_type == 0) {
        SCN[i]->update_single_channel_nucleus_stabilize(steps);
      } else if (nucleus_type == 1) {
        MCN[i]->update_multi_channels_nucleus_stabilize(steps);
      } else {
        std::cout << "stabilisation not implemented (wrong nucleus type)" << std::endl;
      }
  }
}


void BCBG2::updateSingleChannelNucleusTsiro(int steps, float value)
{
  // make the BG modeled in Tsirogiannis et al. 2010 advance of one time step
  int i,s;
  float rnd_nb;
  for (s=0;s<steps;s++) {
    if (value < 0) {
      rnd_nb=_uni_one_gaussian()*2.0f;
      while (rnd_nb < 0) {
        rnd_nb=_uni_one_gaussian()*2.0f;
      }
    } else {
      rnd_nb = value;
    }
    SCN[0]->set_S(rnd_nb,previous_tmod);
    for (i=1; i < n_nuclei; i++) {
      SCN[i]->update_single_channel_nucleus_evo();
    }
    previous_tmod = tmod;
    tmod++;
    if (tmod == max_tau) {
      tmod = 0;
    }
  }
}

float BCBG2::updateSingleChannelNucleusTsiro(int steps)
{
  // draw a random number in a similar way to Tsirogiannis et al. 2010
  float rnd_nb;
  rnd_nb=_uni_one_gaussian()*2.0f;
  while (rnd_nb < 0) {
    rnd_nb=_uni_one_gaussian()*2.0f;
  }
  this->updateSingleChannelNucleusTsiro(steps, rnd_nb);
  return rnd_nb;
}

void BCBG2::hammerizeSingleChannelNucleus(int steps)
{
  // this implements a robustness test to check if the solution can handle big deviation in its inputs
  int i,s;
  float rnd_nb;
  float hammer_nb;
  rand_t _uni_gaussian_hammer(base_generator_t(69u), boost::normal_distribution<>(0.0f, 0.1f));
  rand_t Cortical_input(base_generator_t(42u), boost::normal_distribution<>(2., 1.));
  rand_t Thalamic_input(base_generator_t(51u), boost::normal_distribution<>(4., 1.));

  for (s=0;s<steps;s++) {
    rnd_nb=Cortical_input();
    while (rnd_nb < 0) {
      rnd_nb=Cortical_input();
    }
    SCN[0]->set_S(rnd_nb,previous_tmod);
    rnd_nb=Thalamic_input();
    while (rnd_nb < 0) {
      rnd_nb=Thalamic_input();
    }
    SCN[1]->set_S(rnd_nb,previous_tmod);
    for (i=2; i < n_nuclei; i++) {
      hammer_nb=_uni_gaussian_hammer();
      if (i == MSN_N || i == FSI_N || i == MSND1_N || i == MSND2_N) {
        SCN[i]->set_S(SCN[i]->get_S(previous_tmod)+hammer_nb/10.,previous_tmod); // a variation of 0.1 Hz could be meaningful, so we divide by 10
      } else {
        SCN[i]->set_S(SCN[i]->get_S(previous_tmod)+hammer_nb,previous_tmod);
      }
    }
    for (i=2; i < n_nuclei; i++) {
      SCN[i]->update_single_channel_nucleus_evo();
    }
    previous_tmod = tmod;
    tmod++;
    if (tmod == max_tau) {
      tmod = 0;
    }
  }
}

void BCBG2::hammerizeSingleChannelNucleusTsiro(int steps)
{
// similar to above, but tailored for Tsirogiannis et al 2010
  int i,s;
  float rnd_nb;
  float hammer_nb;
  rand_t _uni_gaussian_hammer(base_generator_t(69u), boost::normal_distribution<>(0.0f, 1.0f));
  rand_t Cortical_input(base_generator_t(42u), boost::normal_distribution<>(3., 1.));

  for (s=0;s<steps;s++) {
    rnd_nb=Cortical_input();
    while (rnd_nb < 0) {
      rnd_nb=Cortical_input();
    }
    SCN[0]->set_S(rnd_nb,previous_tmod);
    for (i=1; i < n_nuclei; i++) {
      hammer_nb=_uni_gaussian_hammer();
      if (i == MSN_N || i == FSI_N || i == MSND1_N || i == MSND2_N) {
        SCN[i]->set_S(SCN[i]->get_S(previous_tmod)+hammer_nb/10.,previous_tmod); // a variation of 0.1 Hz could be meaningful, so we divide by 10
      } else {
        SCN[i]->set_S(SCN[i]->get_S(previous_tmod)+hammer_nb,previous_tmod);
      }
    }
    for (i=1; i < n_nuclei; i++) {
      SCN[i]->update_single_channel_nucleus_evo();
    }
    previous_tmod = tmod;
    tmod++;
    if (tmod == max_tau) {
      tmod = 0;
    }
  }
}

void BCBG2::updateSingleChannelNucleusDebug(float time_of_peak_in_s, float dt, float nb_of_s)
{
  // not intended for actually use
  int i,s;
  for (s=0;s<(int)(nb_of_s/dt+0.5);s++) {
    if (s<((time_of_peak_in_s+0.001)/dt) and s>((time_of_peak_in_s)/dt)) {
      SCN[0]->set_S(100);
    } else {
      SCN[0]->set_S(0);
    }
    for (i=1; i < n_nuclei; i++) {
      SCN[i]->update_single_channel_nucleus_vana();
    }
    for (i=0; i < n_nuclei; i++) {
      std::cerr << s << " : " << SCN[i]->get_name() << "," << SCN[i]->get_S() << std::endl;
    }
    std::cout << s*dt << "," << SCN[1]->get_S() << std::endl;
    tmod++;
    if (tmod == max_tau) {
      tmod = 0;
    }
  }
}

void BCBG2::updateSingleChannelNucleusVana(float input, int bunch)
{
  // make the BG modeled in Van Albada et al. 2009 advance of one time step
  int i,j,s;
  float rnd_nb=10.;
  float worse_Hp = 0;
  for (s=0;s<bunch;s++) {
    if ((1000*tmod)%((int)(1./dt+0.5)) == 0) {
      rnd_nb=_uni_one_gaussian()*2.0f;
      while (rnd_nb < 0) {
        rnd_nb=_uni_one_gaussian()*2.0f;
      }
    }
    SCN[0]->set_S(rnd_nb,(tmod+max_tau)%max_tau);
    for (i=1; i < n_nuclei; i++) {
      SCN[i]->update_single_channel_nucleus_damping_vana();
    }
    for (i=1; i < n_nuclei; i++) {
      SCN[i]->update_single_channel_nucleus_vana();
    }
    tmod++;
    if (tmod == max_tau) {
      tmod = 0;
    }
  }
}

void BCBG2::testPSPSingleChannelNucleus(int pulse, int steps)
{
  // debug function simulating a transient pulse, not intended for actual use
  int s;
  for (s=0;s<5;s++) {
		SCN[0]->set_S(0);
		SCN[1]->update_single_channel_nucleus();
    tmod++;
    if (tmod >= max_tau) {
      tmod = 0;
    }
  }
  SCN[0]->set_S(pulse);
  SCN[1]->update_single_channel_nucleus();
  tmod++;
	if (tmod >= max_tau) {
		tmod = 0;
	}
  for (s=0;s<steps;s++) {
		SCN[0]->set_S(0);
		SCN[1]->update_single_channel_nucleus();
    tmod++;
    if (tmod >= max_tau) {
      tmod = 0;
    }
  }
}

void BCBG2::update(float dt)
{
  // simple wrapper to update a single-channel nucleus
  for (int i=1; i < n_nuclei; i++) {
    if (nucleus_type == 0) {
      SCN[i]->update_single_channel_nucleus();
    } else {
      std::cerr << "ERROR UNKNOWN TYPE in updating" << std::endl;
    }
  }
  previous_tmod = tmod;
  tmod++;
  if (tmod == max_tau) {
    tmod = 0;
  }
}

void BCBG2::updateNucleusCells(int steps)
{
  // multi-channel counterpart of the function above
  int i,s;
  for (s=0;s<steps;s++) {
    for (i=0; i < n_nuclei; i++) {
      NC[i]->prepare_next_iteration();
    }
    NC[CTX_N]->set_rate(2.);
    NC[CMPf_N]->set_rate(4.);
    for (i=0; i < n_nuclei; i++) {
      NC[i]->update_cells();
    }
    previous_tmod = tmod;
    tmod++;
    treal+= dt;
    if (tmod == max_tau) {
      tmod = 0;
    }
  }
}

void BCBG2::autopiloteNucleusCells(int steps, float* values)
{
  // Nucleus modeled with cells is not finished - do not use
  int i,s;
  for (s=0;s<steps;s++) {
    for (i=0; i < n_nuclei; i++) {
      NC[i]->prepare_next_iteration();
    }
    for (i=0; i < n_nuclei; i++) {
      if (values[i] >= 0) {
        NC[i]->set_rate(values[i]);
      }
    }
    for (i=0; i < n_nuclei; i++) {
      NC[i]->update_cells();
    }
    previous_tmod = tmod;
    tmod++;
    treal += dt;
    if (tmod == max_tau) {
      tmod = 0;
    }
  }
}


void BCBG2::initialize(int nucleus_type, int tau_max, float dt)
{
  // resize the array of nuclei
  tmod = 0;
  previous_tmod = tau_max-1;
  max_tau = tau_max;
  treal = 0;
  if (nucleus_type == 0) {
    SCN.resize(0);
  } else if (nucleus_type == 1) {
    MCN.resize(0);
  } else if (nucleus_type == 2) {
    NC.resize(0);
  } else {
    std::cerr << "ERROR UNKNOWN TYPE in initialization" << std::endl;
  }
  this->nucleus_type = nucleus_type;
  n_nuclei = 0;
  this->dt = dt;
}

// the three procedures below trigger the proper initialization procedure, according to the model type
// the specific implementation is in the files singlechannelnucleus.cpp, multichannelsnucleus.cpp and cells.cpp

BCBG2::NucleusCells * BCBG2::add_nucleus_cells (int cells_number, float Smax, float Sini, float vh, float k, float V_rest, float R, float Rm, float Ri, float* dists, float* diams, int compartment_nb, const char* id)
{
  NC.resize(n_nuclei+1);
  NC[n_nuclei] = new NucleusCells(*this);
  NC[n_nuclei]->set_dt(this->dt); //only use: set the internal variable dt for the nucleus
  NC[n_nuclei]->initialize_cells(cells_number,Smax,Sini,vh,k,V_rest,R,Rm,Ri,dists,diams,compartment_nb,id,n_nuclei);
  return NC[n_nuclei++];
}

BCBG2::MultiChannelsNucleus * BCBG2::add_multi_channels_nucleus(int channels_number, float Smax, float Sini, float vh, float k, int connectivity_type, float* dists, float* diams, int compartment_nb, const char* id)
{
  MCN.resize(n_nuclei+1);
  MCN[n_nuclei] = new MultiChannelsNucleus(*this);
  MCN[n_nuclei]->initialize_multi_channels_nucleus(channels_number,Smax,Sini,vh,k,connectivity_type,dists,diams,compartment_nb,id);
  MCN[n_nuclei]->set_dt(this->dt); //only use: set the internal variable dt for the nucleus
  return MCN[n_nuclei++];
}

BCBG2::SingleChannelNucleus * BCBG2::add_single_channel_nucleus(float Smax, float Sini, float vh, float k, int damping, float* dists, float* diams, int compartment_nb, const char* id)
{
  SCN.resize(n_nuclei+1);
  SCN[n_nuclei] = new SingleChannelNucleus(*this);
  SCN[n_nuclei]->initialize_single_channel_nucleus(Smax,Sini,vh,k,damping,dists,diams,compartment_nb,id);
  SCN[n_nuclei]->set_dt(this->dt); //only use: set the internal variable dt for the nucleus
  return SCN[n_nuclei++];
}

