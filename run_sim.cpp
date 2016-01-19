#include "bcbg2.hpp"

#include "constants.hpp"
#include "run_sim.hpp"

void _add_to_fixed_l_list(int size, std::list<std::vector<float> >& l, std::vector<float>& e)
{
  if (l.size() >= size) {
    l.pop_back();
  }
  l.push_front(e);
}

bool _conv(std::list<std::vector<float> >& l)
{
  std::list<std::vector<float> >::const_iterator it1 = l.begin();
  std::list<std::vector<float> >::const_iterator it2 = ++l.begin();
  while (it2 != l.end()) {
    for (size_t i = 0; i < it1->size(); ++i) {
      if (fabs((*it1)[i] - (*it2)[i]) > 1e-3) {
        return false;
      }
    }
    ++it1;
    ++it2;
  }
  return true;
}

int _run_sim(
    float max_duration,
    float bunch_duration,
    float dt,
    std::vector<float> &activ,
    std::vector<float> &params,
    std::vector<int> &delay,
    std::vector<float> &means,
    int verbose)
{
  std::vector<float> cs;
  cs.assign(ARGS_NUMBER, 0.);
  MemoryBCBG2 mem = {-1};
  return _run_sim(max_duration,bunch_duration,dt,activ,cs,params,delay,means,verbose,1,1,1,mem,0);
}

int _run_sim(
    float max_duration,
    float bunch_duration,
    float dt,
    std::vector<float> &a,
    std::vector<float> &c,
    std::vector<float> &p,
    std::vector<int> &de,
    std::vector<float> &means,
    int verbose,
    int ch_n,
    int msn_separation,
    int do_checks,
    MemoryBCBG2& mem,
    int integration_method)
{ 
  // This is the work-horse functions which builds the whole BG model, and computes in a loop the chosen experiment (selected with the variable verbose)

  int return_value = 0;
  float minipsilon = 1e-4;
  int max_tau = 100+(*std::max_element(de.begin(),de.end()));

  float last_GPi = -10; float last_GPe = 10; float last_STN = 10000;
  float last_D1 = -10; float last_D2 = 10; float last_FS = 10000;

  std::list<std::vector<float> > last_out;
  last_out.resize(0);

  BCBG2* bg = new BCBG2();

  if (verbose==1)
    std::cout << "max_tau is " << max_tau << " & dt is " << dt << std::endl;


  float A_AMPA = 2.718281828; // 1 mV
  float D_AMPA = 0.5437; //1/5 *e
  int S_AMPA = 1;

float A_NMDA = A_AMPA / 40.; // we suppose a NMDA effect 40 times weaker. This corresponds to the ratio AMPA/NMDA when 100 synapses are concerned see Bouteiller 2008 p.190
float D_NMDA = 0.02718; // we suppose a NMDA effect 100 ms long, which is 20 times longer. See Destexhe 1998.
  // overall, the NMDA is 2 times weaker than AMPA, which is coherent with desactivations studies
  int S_NMDA = 1;

  float A_GABAA = 2.718281828; // 1mV
  A_GABAA *= 0.25;
  float D_GABAA = 0.5437; // 1/5 *e
  int S_GABAA = -1;

  float A_GABAB = A_GABAA / 100.; // we suppose a GABA_B effect 100 times weaker.
  float D_GABAB = 0.0544*2.; // we suppose a GABA_B effect 50 ms/2 long.
  // overall, the GABA_B is 10 times weaker than GABA_A.
  int S_GABAB = -1;


// the following page handles how Tsirogiannis et al 2010 model Parkinson's disease (by changing the alpha and delta values)
  float alpha_tsiro = 1.0 + ((float) integration_method)*0.01; // last argument to this function (and second on commandline)
  float delta_tsiro = 1.0 + ((float) do_checks)*0.01;;         // ante-ante-last arg to this function (and first on commandline)

  float A_AMPASTR1  = A_AMPA;
  float A_NMDASTR1  = A_NMDA;
  float A_GABAASTR1 = A_GABAA;
  float D_AMPASTR1  = D_AMPA;
  float D_NMDASTR1  = D_NMDA;
  float D_GABAASTR1 = D_GABAA;

  float A_AMPASTR2  = A_AMPA;
  float A_NMDASTR2  = A_NMDA;
  float A_GABAASTR2 = A_GABAA;
  float D_AMPASTR2  = D_AMPA;
  float D_NMDASTR2  = D_NMDA;
  float D_GABAASTR2 = D_GABAA;

  float A_AMPASTR  = A_AMPA;
  float A_NMDASTR  = A_NMDA;
  float A_GABAASTR = A_GABAA;
  float D_AMPASTR  = D_AMPA;
  float D_NMDASTR  = D_NMDA;
  float D_GABAASTR = D_GABAA;

  // A_AMPA  *= alpha_tsiro;
  // A_NMDA  *= alpha_tsiro;
  // A_GABAA *= alpha_tsiro;
  // D_AMPA  *= delta_tsiro;
  // D_NMDA  *= delta_tsiro;
  // D_GABAA *= delta_tsiro;


  float diameters_dummy[] = {};           int comp_nb_dummy=0;  // CTX & CMPf (dummy)
  //t float diameters_msn[] = {10.*1e-6, 3.3*1e-6, 1.*1e-6};  // MSN
  float diameters_msn[] = {1.*1e-6};    int comp_nb_msn=1;  // MSN (simple version)
  float diameters_fsi[] = {1.5*1e-6};    int comp_nb_fsi=1;  // FSI (simple version)
  float diameters_stn[] = {1.5*1e-6};   int comp_nb_stn=1;  // STN (simple version)
  //t float diameters_gpi[] = {20.*1e-6, 3.6*1e-6, 1.1*1e-6}; // GPi
  float diameters_gpi[] = {1.2*1e-6}; int comp_nb_gpi=1;  // GPi (simple version)
  //t float diameters_gpe[] = {20.*1e-6, 3.6*1e-6, 1.1*1e-6}; // GPe
float diameters_gpe[] = {1.7*1e-6}; int comp_nb_gpe=1;  // GPe (simple version)


float distances_dummy[] = {};                                // CTX & CMPf (dummy)
//t float distances_msn[] = {10.*1e-6, 11.*1e-6, 600.*1e-6};   // MSN
float distances_msn[] = {619.*1e-6};                       // MSN (simple version)
float distances_fsi[] = {961.*1e-6};                       // FSI (simple version)
float distances_stn[] = {750.*1e-6};                       // STN (simple version)
//t float distances_gpi[] = {10.*1e-6, 139.*1e-6, 1000.*1e-6}; // GPi
float distances_gpi[] = {1132.*1e-6};                      // GPi (simple version)
//t float distances_gpe[] = {10.*1e-6, 139.*1e-6, 700.*1e-6};  // GPe
float distances_gpe[] = {865.*1e-6};                       // GPe (simple version)

#if defined(MULTICHANNELSMODEL)
bg->initialize(1,max_tau,dt);
BCBG2::MultiChannelsNucleus* CTX  = bg->add_multi_channels_nucleus(ch_n,0.,2.,0.,0.,0,distances_dummy,diameters_dummy,comp_nb_dummy,"CTX");
BCBG2::MultiChannelsNucleus* CMPf = bg->add_multi_channels_nucleus(ch_n,0.,4.,0.,0.,0,distances_dummy,diameters_dummy,comp_nb_dummy,"CMPf");
BCBG2::MultiChannelsNucleus* MSN;
BCBG2::MultiChannelsNucleus* MSND1;
BCBG2::MultiChannelsNucleus* MSND2;
if (msn_separation) {
  MSND1  = bg->add_multi_channels_nucleus(ch_n,300. ,0.1,(p[THETA_MSN]*25.+5.),0.26,0,distances_msn,diameters_msn,comp_nb_msn,"MSND1");
} else {
  MSN  = bg->add_multi_channels_nucleus(ch_n,300. ,0.1,(p[THETA_MSN]*25.+5.),0.26,0,distances_msn,diameters_msn,comp_nb_msn,"MSN");
}
BCBG2::MultiChannelsNucleus* FSI  = bg->add_multi_channels_nucleus(ch_n,p[FSI_SMAX],10.,(p[THETA_FSI]*25.+5.),0.26,0,distances_fsi,diameters_fsi,comp_nb_fsi,"FSI");
BCBG2::MultiChannelsNucleus* STN  = bg->add_multi_channels_nucleus(ch_n,300.,20.,(p[THETA_STN]*25.+5.),0.26,0,distances_stn,diameters_stn,comp_nb_stn,"STN");
BCBG2::MultiChannelsNucleus* GPe  = bg->add_multi_channels_nucleus(ch_n,400.,65.,(p[THETA_GPe]*25.+5.),0.26,0,distances_gpe,diameters_gpe,comp_nb_gpe,"GPe");
BCBG2::MultiChannelsNucleus* GPi  = bg->add_multi_channels_nucleus(ch_n,400.,70.,(p[THETA_GPi]*25.+5.),0.26,0,distances_gpi,diameters_gpi,comp_nb_gpi,"GPi");
if (msn_separation) {
  MSND2  = bg->add_multi_channels_nucleus(ch_n,300. ,0.1,(p[THETA_MSN]*25.+5.),0.26,0,distances_msn,diameters_msn,comp_nb_msn,"MSND2");
}

#ifdef ITPTCTX
BCBG2::MultiChannelsNucleus* CTXPT  = bg->add_multi_channels_nucleus(ch_n,0.,15.,0.,0.,0,distances_dummy,diameters_dummy,comp_nb_dummy,"CTXPT");
#endif

#elif defined(CELLSMODEL)
// not used

#else

bg->initialize(0,max_tau,dt);
ch_n = 1;
BCBG2::SingleChannelNucleus* CTX  = bg->add_single_channel_nucleus(0.,2.,0.,0.,0,distances_dummy,diameters_dummy,comp_nb_dummy,"CTX");
BCBG2::SingleChannelNucleus* CMPf = bg->add_single_channel_nucleus(0.,4.,0.,0.,0,distances_dummy,diameters_dummy,comp_nb_dummy,"CMPf");
BCBG2::SingleChannelNucleus* MSN;
BCBG2::SingleChannelNucleus* MSND1;
BCBG2::SingleChannelNucleus* MSND2;
if (msn_separation) {
  MSND1  = bg->add_single_channel_nucleus(300.,0.1,(p[THETA_MSN]*25.+5.),0.26,0,distances_msn,diameters_msn,comp_nb_msn,"MSND1");
} else {
  MSN  = bg->add_single_channel_nucleus(300.,0.1,(p[THETA_MSN]*25.+5.),0.26,0,distances_msn,diameters_msn,comp_nb_msn,"MSN");
}
BCBG2::SingleChannelNucleus* FSI  = bg->add_single_channel_nucleus(p[FSI_SMAX],10.,(p[THETA_FSI]*25.+5.),0.26,0,distances_fsi,diameters_fsi,comp_nb_fsi,"FSI");
BCBG2::SingleChannelNucleus* STN  = bg->add_single_channel_nucleus(300.,20.,(p[THETA_STN]*25.+5.),0.26,0,distances_stn,diameters_stn,comp_nb_stn,"STN");
BCBG2::SingleChannelNucleus* GPe  = bg->add_single_channel_nucleus(400.,65.,(p[THETA_GPe]*25.+5.),0.26,0,distances_gpe,diameters_gpe,comp_nb_gpe,"GPe");
BCBG2::SingleChannelNucleus* GPi  = bg->add_single_channel_nucleus(400.,70.,(p[THETA_GPi]*25.+5.),0.26,0,distances_gpi,diameters_gpi,comp_nb_gpi,"GPi");
if (msn_separation) {
  MSND2  = bg->add_single_channel_nucleus(300. ,0.1,(p[THETA_MSN]*25.+5.),0.26,0,distances_msn,diameters_msn,comp_nb_msn,"MSND2");
}

#ifdef ITPTCTX
BCBG2::SingleChannelNucleus* CTXPT  = bg->add_single_channel_nucleus(0.,15.,0.,0.,0,distances_dummy,diameters_dummy,comp_nb_dummy,"CTXPT");
#endif

#endif


bool with_gabab = false;


if (msn_separation) {
#ifdef ITPTCTX
MSND1->set_afferent(A_AMPASTR1 , D_AMPASTR1 , S_AMPA , p[CTXPT_MSN]*a[CTXPT_MSN_AMPA] , de[CTXPT_MSN], c[CTXPT_MSN], p[DIST_CTXPT_MSN], CTXPT);
MSND1->set_afferent(A_NMDASTR1, D_NMDASTR1 , S_NMDA , p[CTXPT_MSN]*a[CTXPT_MSN_NMDA] , de[CTXPT_MSN], c[CTXPT_MSN], p[DIST_CTXPT_MSN], CTXPT);
#endif
  MSND1->set_afferent(A_AMPASTR1 , D_AMPASTR1 , S_AMPA , p[CTX_MSN]*a[CTX_MSN_AMPA], de[CTX_MSN], c[CTX_MSN], p[DIST_CTX_MSN], CTX);
  MSND1->set_afferent(A_NMDASTR1, D_NMDASTR1 , S_NMDA , p[CTX_MSN]*a[CTX_MSN_NMDA], de[CTX_MSN], c[CTX_MSN], p[DIST_CTX_MSN], CTX);
  MSND1->set_afferent(A_AMPASTR1 , D_AMPASTR1 , S_AMPA , p[CMPf_MSN]*a[CMPf_MSN_AMPA], de[CMPf_MSN], c[CMPf_MSN], p[DIST_CMPf_MSN], CMPf);
  MSND1->set_afferent(A_NMDASTR1, D_NMDASTR1 , S_NMDA , p[CMPf_MSN]*a[CMPf_MSN_NMDA] , de[CMPf_MSN], c[CMPf_MSN], p[DIST_CMPf_MSN], CMPf);
  MSND1->set_afferent(A_AMPASTR1 , D_AMPASTR1 , S_AMPA , p[STN_MSN]*a[STN_MSN_AMPA] , de[STN_MSN], c[STN_MSN], p[DIST_STN_MSN], STN);
  MSND1->set_afferent(A_NMDASTR1, D_NMDASTR1 , S_NMDA , p[STN_MSN]*a[STN_MSN_NMDA] , de[STN_MSN], c[STN_MSN], p[DIST_STN_MSN], STN);
  MSND1->set_afferent(A_GABAASTR1, D_GABAASTR1, S_GABAA, p[GPe_MSN]*a[GPe_MSN_GABAA], de[GPe_MSN], c[GPe_MSN], p[DIST_GPe_MSN], GPe);
  MSND1->set_afferent(A_GABAASTR1, D_GABAASTR1, S_GABAA, p[FSI_MSN]*a[FSI_MSN_GABAA], de[FSI_MSN], c[FSI_MSN], p[DIST_FSI_MSN], FSI);
  MSND1->set_afferent(A_GABAASTR1, D_GABAASTR1, S_GABAA, p[MSN_MSN]*a[MSN_MSN_GABAA], de[MSN_MSN], c[MSN_MSN], p[DIST_MSN_MSN], MSND1);
  MSND1->set_afferent(A_GABAASTR1, D_GABAASTR1, S_GABAA, p[MSN_MSN]*a[MSN_MSN_GABAA], de[MSN_MSN], c[MSN_MSN], p[DIST_MSN_MSN], MSND2);
#ifdef ITPTCTX
MSND2->set_afferent(A_AMPASTR2 , D_AMPASTR2 , S_AMPA , p[CTXPT_MSN]*a[CTXPT_MSN_AMPA] , de[CTXPT_MSN], c[CTXPT_MSN], p[DIST_CTXPT_MSN], CTXPT);
MSND2->set_afferent(A_NMDASTR2, D_NMDASTR2 , S_NMDA , p[CTXPT_MSN]*a[CTXPT_MSN_NMDA] , de[CTXPT_MSN], c[CTXPT_MSN], p[DIST_CTXPT_MSN], CTXPT);
#endif
  MSND2->set_afferent(A_AMPASTR2 , D_AMPASTR2 , S_AMPA , p[CTX_MSN]*a[CTX_MSN_AMPA], de[CTX_MSN], c[CTX_MSN], p[DIST_CTX_MSN], CTX);
  MSND2->set_afferent(A_NMDASTR2, D_NMDASTR2 , S_NMDA , p[CTX_MSN]*a[CTX_MSN_NMDA], de[CTX_MSN], c[CTX_MSN], p[DIST_CTX_MSN], CTX);
  MSND2->set_afferent(A_AMPASTR2 , D_AMPASTR2 , S_AMPA , p[CMPf_MSN]*a[CMPf_MSN_AMPA], de[CMPf_MSN], c[CMPf_MSN], p[DIST_CMPf_MSN], CMPf);
  MSND2->set_afferent(A_NMDASTR2, D_NMDASTR2 , S_NMDA , p[CMPf_MSN]*a[CMPf_MSN_NMDA], de[CMPf_MSN], c[CMPf_MSN], p[DIST_CMPf_MSN], CMPf);
  MSND2->set_afferent(A_AMPASTR2 , D_AMPASTR2 , S_AMPA , p[STN_MSN]*a[STN_MSN_AMPA], de[STN_MSN], c[STN_MSN], p[DIST_STN_MSN], STN);
  MSND2->set_afferent(A_NMDASTR2, D_NMDASTR2 , S_NMDA , p[STN_MSN]*a[STN_MSN_NMDA], de[STN_MSN], c[STN_MSN], p[DIST_STN_MSN], STN);
  MSND2->set_afferent(A_GABAASTR2, D_GABAASTR2, S_GABAA, p[GPe_MSN]*a[GPe_MSN_GABAA], de[GPe_MSN], c[GPe_MSN], p[DIST_GPe_MSN], GPe);
  MSND2->set_afferent(A_GABAASTR2, D_GABAASTR2, S_GABAA, p[FSI_MSN]*a[FSI_MSN_GABAA], de[FSI_MSN], c[FSI_MSN], p[DIST_FSI_MSN], FSI);
  MSND2->set_afferent(A_GABAASTR2, D_GABAASTR2, S_GABAA, p[MSN_MSN]*a[MSN_MSN_GABAA], de[MSN_MSN], c[MSN_MSN], p[DIST_MSN_MSN], MSND1);
  MSND2->set_afferent(A_GABAASTR2, D_GABAASTR2, S_GABAA, p[MSN_MSN]*a[MSN_MSN_GABAA], de[MSN_MSN], c[MSN_MSN], p[DIST_MSN_MSN], MSND2);
  if (with_gabab) {
    MSND1->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[GPe_MSN]*a[GPe_MSN_GABAB], de[GPe_MSN], c[GPe_MSN], p[DIST_GPe_MSN], GPe);
    MSND1->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[FSI_MSN]*a[FSI_MSN_GABAB], de[FSI_MSN], c[FSI_MSN], p[DIST_FSI_MSN], FSI);
    MSND1->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[MSN_MSN]*a[MSN_MSN_GABAB], de[MSN_MSN], c[MSN_MSN], p[DIST_MSN_MSN], MSND1);
    MSND1->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[MSN_MSN]*a[MSN_MSN_GABAB], de[MSN_MSN], c[MSN_MSN], p[DIST_MSN_MSN], MSND2);
    MSND2->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[GPe_MSN]*a[GPe_MSN_GABAB], de[GPe_MSN], c[GPe_MSN], p[DIST_GPe_MSN], GPe);
    MSND2->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[FSI_MSN]*a[FSI_MSN_GABAB], de[FSI_MSN], c[FSI_MSN], p[DIST_FSI_MSN], FSI);
    MSND2->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[MSN_MSN]*a[MSN_MSN_GABAB], de[MSN_MSN], c[MSN_MSN], p[DIST_MSN_MSN], MSND1);
    MSND2->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[MSN_MSN]*a[MSN_MSN_GABAB], de[MSN_MSN], c[MSN_MSN], p[DIST_MSN_MSN], MSND2);
  }
} else {
#ifdef ITPTCTX
MSN->set_afferent(A_AMPASTR , D_AMPASTR , S_AMPA , p[CTXPT_MSN]*a[CTXPT_MSN_AMPA] , de[CTXPT_MSN], c[CTXPT_MSN], p[DIST_CTXPT_MSN], CTXPT);
MSN->set_afferent(A_NMDASTR, D_NMDASTR , S_NMDA , p[CTXPT_MSN]*a[CTXPT_MSN_NMDA] , de[CTXPT_MSN], c[CTXPT_MSN], p[DIST_CTXPT_MSN], CTXPT);
#endif
  MSN->set_afferent(A_AMPASTR , D_AMPASTR , S_AMPA , p[CTX_MSN]*a[CTX_MSN_AMPA] , de[CTX_MSN], c[CTX_MSN], p[DIST_CTX_MSN], CTX);
  MSN->set_afferent(A_NMDASTR, D_NMDASTR , S_NMDA , p[CTX_MSN]*a[CTX_MSN_NMDA] , de[CTX_MSN], c[CTX_MSN], p[DIST_CTX_MSN], CTX);
  MSN->set_afferent(A_AMPASTR , D_AMPASTR , S_AMPA , p[CMPf_MSN]*a[CMPf_MSN_AMPA] , de[CMPf_MSN], c[CMPf_MSN], p[DIST_CMPf_MSN], CMPf);
  MSN->set_afferent(A_NMDASTR, D_NMDASTR , S_NMDA , p[CMPf_MSN]*a[CMPf_MSN_NMDA] , de[CMPf_MSN], c[CMPf_MSN], p[DIST_CMPf_MSN], CMPf);
  MSN->set_afferent(A_AMPASTR , D_AMPASTR , S_AMPA , p[STN_MSN]*a[STN_MSN_AMPA] , de[STN_MSN], c[STN_MSN], p[DIST_STN_MSN], STN);
  MSN->set_afferent(A_NMDASTR, D_NMDASTR , S_NMDA , p[STN_MSN]*a[STN_MSN_NMDA] , de[STN_MSN], c[STN_MSN], p[DIST_STN_MSN], STN);
  MSN->set_afferent(A_GABAASTR, D_GABAASTR, S_GABAA, p[GPe_MSN]*a[GPe_MSN_GABAA], de[GPe_MSN], c[GPe_MSN], p[DIST_GPe_MSN], GPe);
  MSN->set_afferent(A_GABAASTR, D_GABAASTR, S_GABAA, p[FSI_MSN]*a[FSI_MSN_GABAA], de[FSI_MSN], c[FSI_MSN], p[DIST_FSI_MSN], FSI);
  MSN->set_afferent(A_GABAASTR, D_GABAASTR, S_GABAA, p[MSN_MSN]*a[MSN_MSN_GABAA], de[MSN_MSN], c[MSN_MSN], p[DIST_MSN_MSN], MSN);
  if (with_gabab) {
    MSN->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[GPe_MSN]*a[GPe_MSN_GABAB], de[GPe_MSN], c[GPe_MSN], p[DIST_GPe_MSN], GPe);
    MSN->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[FSI_MSN]*a[FSI_MSN_GABAB], de[FSI_MSN], c[FSI_MSN], p[DIST_FSI_MSN], FSI);
    MSN->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[MSN_MSN]*a[MSN_MSN_GABAB], de[MSN_MSN], c[MSN_MSN], p[DIST_MSN_MSN], MSN);
  }
}

#ifdef ITPTCTX
FSI->set_afferent(A_AMPASTR , D_AMPASTR , S_AMPA , p[CTXPT_FSI]*a[CTXPT_FSI_AMPA] , de[CTXPT_FSI], c[CTXPT_FSI], p[DIST_CTXPT_FSI], CTXPT);
FSI->set_afferent(A_NMDASTR, D_NMDASTR , S_NMDA , p[CTXPT_FSI]*a[CTXPT_FSI_NMDA] , de[CTXPT_FSI], c[CTXPT_FSI], p[DIST_CTXPT_FSI], CTXPT);
#endif
FSI->set_afferent(A_AMPASTR , D_AMPASTR , S_AMPA , p[CTX_FSI]*a[CTX_FSI_AMPA] , de[CTX_FSI], c[CTX_FSI], p[DIST_CTX_FSI], CTX);
FSI->set_afferent(A_NMDASTR, D_NMDASTR , S_NMDA , p[CTX_FSI]*a[CTX_FSI_NMDA] , de[CTX_FSI], c[CTX_FSI], p[DIST_CTX_FSI], CTX);
FSI->set_afferent(A_AMPASTR , D_AMPASTR , S_AMPA , p[CMPf_FSI]*a[CMPf_FSI_AMPA] , de[CMPf_FSI], c[CMPf_FSI], p[DIST_CMPf_FSI], CMPf);
FSI->set_afferent(A_NMDASTR, D_NMDASTR , S_NMDA , p[CMPf_FSI]*a[CMPf_FSI_NMDA] , de[CMPf_FSI], c[CMPf_FSI], p[DIST_CMPf_FSI], CMPf);
FSI->set_afferent(A_AMPASTR , D_AMPASTR , S_AMPA , p[STN_FSI]*a[STN_FSI_AMPA] , de[STN_FSI], c[STN_FSI], p[DIST_STN_FSI], STN);
FSI->set_afferent(A_NMDASTR, D_NMDASTR , S_NMDA , p[STN_FSI]*a[STN_FSI_NMDA] , de[STN_FSI], c[STN_FSI], p[DIST_STN_FSI], STN);
FSI->set_afferent(A_GABAASTR, D_GABAASTR, S_GABAA, p[GPe_FSI]*a[GPe_FSI_GABAA], de[GPe_FSI], c[GPe_FSI], p[DIST_GPe_FSI], GPe);
FSI->set_afferent(A_GABAASTR, D_GABAASTR, S_GABAA, p[FSI_FSI]*a[FSI_FSI_GABAA], de[FSI_FSI], c[FSI_FSI], p[DIST_FSI_FSI], FSI);
if (with_gabab) {
  FSI->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[GPe_FSI]*a[GPe_FSI_GABAB], de[GPe_FSI], c[GPe_FSI], p[DIST_GPe_FSI], GPe);
  FSI->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[FSI_FSI]*a[FSI_FSI_GABAB], de[FSI_FSI], c[FSI_FSI], p[DIST_FSI_FSI], FSI);
}

STN->set_afferent(A_AMPA, D_AMPA , S_AMPA , p[CMPf_STN]*a[CMPf_STN_AMPA] , de[CMPf_STN], c[CMPf_STN], p[DIST_CMPf_STN], CMPf);
STN->set_afferent(A_NMDA, D_NMDA , S_NMDA , p[CMPf_STN]*a[CMPf_STN_NMDA] , de[CMPf_STN], c[CMPf_STN], p[DIST_CMPf_STN], CMPf);
#ifdef ITPTCTX
STN->set_afferent(A_AMPA , D_AMPA , S_AMPA , p[CTX_STN]*a[CTX_STN_AMPA] , de[CTX_STN], c[CTX_STN], p[DIST_CTX_STN], CTXPT);
STN->set_afferent(A_NMDA , D_NMDA , S_NMDA , p[CTX_STN]*a[CTX_STN_NMDA] , de[CTX_STN], c[CTX_STN], p[DIST_CTX_STN], CTXPT);
#else
STN->set_afferent(A_AMPA , D_AMPA , S_AMPA , p[CTX_STN]*a[CTX_STN_AMPA] , de[CTX_STN], c[CTX_STN], p[DIST_CTX_STN], CTX);
STN->set_afferent(A_NMDA , D_NMDA , S_NMDA , p[CTX_STN]*a[CTX_STN_NMDA] , de[CTX_STN], c[CTX_STN], p[DIST_CTX_STN], CTX);
#endif
STN->set_afferent(A_GABAA, D_GABAA, S_GABAA, p[GPe_STN]*a[GPe_STN_GABAA], de[GPe_STN], c[GPe_STN], p[DIST_GPe_STN], GPe);
if (with_gabab) {
  STN->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[GPe_STN]*a[GPe_STN_GABAB], de[GPe_STN], c[GPe_STN], p[DIST_GPe_STN], GPe);
}

GPe->set_afferent(A_AMPA , D_AMPA , S_AMPA , p[CMPf_GPe]*a[CMPf_GPe_AMPA] , de[CMPf_GPe], c[CMPf_GPe], p[DIST_CMPf_GPe], CMPf);
GPe->set_afferent(A_NMDA, D_NMDA , S_NMDA , p[CMPf_GPe]*a[CMPf_GPe_NMDA] , de[CMPf_GPe], c[CMPf_GPe], p[DIST_CMPf_GPe], CMPf);
GPe->set_afferent(A_AMPA , D_AMPA , S_AMPA , p[STN_GPe]*a[STN_GPe_AMPA] , de[STN_GPe], c[STN_GPe], p[DIST_STN_GPe], STN);
GPe->set_afferent(A_NMDA, D_NMDA , S_NMDA , p[STN_GPe]*a[STN_GPe_NMDA] , de[STN_GPe], c[STN_GPe], p[DIST_STN_GPe], STN);
if (msn_separation) {
  GPe->set_afferent(A_GABAA, D_GABAA, S_GABAA, p[MSN_GPe]*a[MSN_GPe_GABAA], de[MSN_GPe], c[MSN_GPe], p[DIST_MSN_GPe], MSND1);
  GPe->set_afferent(A_GABAA, D_GABAA, S_GABAA, p[MSN_GPe]*a[MSN_GPe_GABAA], de[MSN_GPe], c[MSN_GPe], p[DIST_MSN_GPe], MSND2);
} else {
  GPe->set_afferent(A_GABAA, D_GABAA, S_GABAA, p[MSN_GPe]*a[MSN_GPe_GABAA], de[MSN_GPe], c[MSN_GPe], p[DIST_MSN_GPe], MSN);
}
GPe->set_afferent(A_GABAA, D_GABAA, S_GABAA, p[GPe_GPe]*a[GPe_GPe_GABAA], de[GPe_GPe], c[GPe_GPe], p[DIST_GPe_GPe], GPe);
if (with_gabab) {
  if (msn_separation) {
    GPe->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[MSN_GPe]*a[MSN_GPe_GABAB], de[MSN_GPe], c[MSN_GPe], p[DIST_MSN_GPe], MSND1);
    GPe->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[MSN_GPe]*a[MSN_GPe_GABAB], de[MSN_GPe], c[MSN_GPe], p[DIST_MSN_GPe], MSND2);
  } else {
    GPe->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[MSN_GPe]*a[MSN_GPe_GABAB], de[MSN_GPe], c[MSN_GPe], p[DIST_MSN_GPe], MSN);
  }
  GPe->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[GPe_GPe]*a[GPe_GPe_GABAB], de[GPe_GPe], c[GPe_GPe], p[DIST_GPe_GPe], GPe);
}

GPi->set_afferent(A_AMPA , D_AMPA , S_AMPA , p[CMPf_GPi]*a[CMPf_GPi_AMPA] , de[CMPf_GPi], c[CMPf_GPi], p[DIST_CMPf_GPi], CMPf);
GPi->set_afferent(A_NMDA, D_NMDA , S_NMDA , p[CMPf_GPi]*a[CMPf_GPi_NMDA] , de[CMPf_GPi], c[CMPf_GPi], p[DIST_CMPf_GPi], CMPf);
GPi->set_afferent(A_AMPA , D_AMPA , S_AMPA , p[STN_GPi]*a[STN_GPi_AMPA] , de[STN_GPi], c[STN_GPi], p[DIST_STN_GPi], STN);
GPi->set_afferent(A_NMDA, D_NMDA , S_NMDA , p[STN_GPi]*a[STN_GPi_NMDA] , de[STN_GPi], c[STN_GPi], p[DIST_STN_GPi], STN);
if (msn_separation) {
  GPi->set_afferent(A_GABAA, D_GABAA, S_GABAA, p[MSN_GPi]*a[MSN_GPi_GABAA], de[MSN_GPi], c[MSN_GPi], p[DIST_MSN_GPi], MSND1);
  GPi->set_afferent(A_GABAA, D_GABAA, S_GABAA, p[MSN_GPi]*a[MSN_GPi_GABAA], de[MSN_GPi], c[MSN_GPi], p[DIST_MSN_GPi], MSND2);
} else {
  GPi->set_afferent(A_GABAA, D_GABAA, S_GABAA, p[MSN_GPi]*a[MSN_GPi_GABAA], de[MSN_GPi], c[MSN_GPi], p[DIST_MSN_GPi], MSN);
}
GPi->set_afferent(A_GABAA, D_GABAA, S_GABAA, p[GPe_GPi]*a[GPe_GPi_GABAA], de[GPe_GPi], c[GPe_GPi], p[DIST_GPe_GPi], GPe);
if (with_gabab) {
  if (msn_separation) {
    GPi->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[MSN_GPi]*a[MSN_GPi_GABAB], de[MSN_GPi], c[MSN_GPi], p[DIST_MSN_GPi], MSND1);
    GPi->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[MSN_GPi]*a[MSN_GPi_GABAB], de[MSN_GPi], c[MSN_GPi], p[DIST_MSN_GPi], MSND2);
  } else {
    GPi->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[MSN_GPi]*a[MSN_GPi_GABAB], de[MSN_GPi], c[MSN_GPi], p[DIST_MSN_GPi], MSN);
  }
  GPi->set_afferent(A_GABAB, D_GABAB, S_GABAB, p[GPe_GPi]*a[GPe_GPi_GABAB], de[GPe_GPi], c[GPe_GPi], p[DIST_GPe_GPi], GPe);
}

int bunch = (int)(bunch_duration/dt + 0.5);
int count = 0;
int old_count = 0;
float min_duration = 10.; // simulations longer than 10s
int nb_conv = (int)(2./bunch_duration+0.5); // stabilisation over at least 2 second
#ifdef LIGHTCONV
min_duration = 0.2;
nb_conv = (int)(0.1/bunch_duration+0.5); // stabilisation over at least 0.1 second
#elif defined(LIGHTESTCONV)
min_duration = 0.1;
nb_conv = (int)(0.05/bunch_duration+0.5); // stabilisation over at least 0.05 second
#elif defined(SOUNDCONV)
min_duration = 1.2;
nb_conv = (int)(1.0/bunch_duration+0.5); // stabilisation over at least 1 second
#elif defined(OBJECTCONV)
min_duration = 2.0;
nb_conv = (int)(1.0/bunch_duration+0.5); // stabilisation over at least 1 second
#endif

int contiguous_check_nb = 3;
nb_conv *= (contiguous_check_nb+1);

float hammerization_eval;


old_count += count;
count = 0;

if (verbose==1)
{
  std::cerr << bunch << " iterations (" << bunch*dt << " s) to be done between two tests of convergence, until " << max_duration/dt << " iterations (" << max_duration << "s) are done" << std::endl;
}



  if (mem.tmod_backup != -1) {
    bg->load_all(mem); // as this is not the first run, we load the previous state
  } else {
    bg->stabilize_all((int)(10./dt + 0.5)); // as this is the first run, we try to stabilize as a mean to improve convergence
  }

  float ref_rates[7]={2.,4.,0.1,10.,20.,60.,70.};

  float maxstn =0;

  float freqcount =0;
  float lastdiff  =0;
  float laststn   =0;
  float highval   =0;
  float lowval    =0;

  bool not_nan = true;

  while (
      (
      (
       (
        (
         ! _conv(last_out) || (verbose == 3) || (verbose == 4) || (verbose>4)
        )
        ||
        (count < max_tau)
        ||
        (count*dt < min_duration)
       )
       &&
       ( count*dt < max_duration )
      )
      ||
      (verbose==1 && count*dt < max_duration)
      )
      && not_nan
      )
      {
    if (verbose==2||verbose==3||verbose==4||(verbose>4))
    {
#ifdef MINRECORDT
      if ((old_count+count)*dt <= MINRECORDT) {
        for (int i=0; i<NUCLEUS_NUMBER; i++) { // displays also for the cortex & CMPf
          ref_rates[i] = 0;
          for (int ch_i=0; ch_i<ch_n; ch_i++) {
            ref_rates[i] += bg->get_multi_channels_nucleus(i)->get_S(ch_i);
          }
          ref_rates[i] /= ch_n;
        }
      }
      if ((old_count+count)*dt > MINRECORDT) {
#ifdef MAXRECORDT
      if ((old_count+count)*dt <= MAXRECORDT) {
        if (bg->get_multi_channels_nucleus(STN_N)->get_S(0) - laststn > 0 && lastdiff < 0) {
          freqcount += 1;
          highval += bg->get_multi_channels_nucleus(STN_N)->get_S(0);
        } else if (bg->get_multi_channels_nucleus(STN_N)->get_S(0) - laststn < 0 && lastdiff > 0) {
          lowval += bg->get_multi_channels_nucleus(STN_N)->get_S(0);
        }
        lastdiff   = bg->get_multi_channels_nucleus(STN_N)->get_S(0) - laststn;
        laststn    = bg->get_multi_channels_nucleus(STN_N)->get_S(0);
        if (bg->get_multi_channels_nucleus(STN_N)->get_S(0) > maxstn)
          maxstn = bg->get_multi_channels_nucleus(STN_N)->get_S(0);
#endif
#endif

#ifndef NOOUTPUT

#ifdef MINRECORDT
      std::cerr << (old_count+count)*dt - MINRECORDT;
#else
      std::cerr << (old_count+count)*dt;
#endif
      if (ch_n == 1) {
        for (int i=0; i<NUCLEUS_NUMBER; i++) { // displays also for the cortex & CMPf
          std::cerr << " , " << bg->get_single_channel_nucleus(i)->get_S();
        }
      } else if (ch_n == 0) {
        for (int i=0; i<NUCLEUS_NUMBER; i++) { // displays also for the cortex & CMPf
          std::cerr << " , " << bg->get_nucleus_cells(i)->get_meanS();
        }
      } else {
        for (int i=0; i<NUCLEUS_NUMBER; i++) { // displays also for the cortex & CMPf
          float nmeans = 0;
          for (int ch_i=0; ch_i<ch_n; ch_i++) {
          //for (int ch_i=0; ch_i<3; ch_i++) { // only the 3 first channels...
#ifndef ONLYMEANOUTPUT
#ifdef MINRECORDT
            std::cerr << " , " << bg->get_multi_channels_nucleus(i)->get_S(ch_i);
#else
            std::cerr << " , " << bg->get_multi_channels_nucleus(i)->get_S(ch_i);
#endif
#endif
            nmeans += bg->get_multi_channels_nucleus(i)->get_S(ch_i);
          }
#ifdef MINRECORDT
          std::cerr << " , " << nmeans / ((float)ch_n)/ref_rates[i];
#else
          std::cerr << " , " << nmeans / ((float)ch_n);
#endif
        }
        }
        std::cerr << std::endl;
      //nooutput end
#endif
      }


#ifdef MINRECORDT
      }
#ifdef MAXRECORDT
      }
#endif
#endif

      if (ch_n == 1) {
        if (verbose == 4) {
          bg->updateSingleChannelNucleusWithNoise(bunch,integration_method);
        } else {
          bg->updateSingleChannelNucleus(bunch,integration_method);
        }
        count+=bunch;
#ifdef ITPTCTX
        for (int i=2; i<NUCLEUS_NUMBER-1; i++)
#else
        for (int i=2; i<NUCLEUS_NUMBER; i++)
#endif
        { // do not store for the cortex & CMPf
          means[i] = bg->get_single_channel_nucleus(i)->get_S();
          if (isnan(means[i])) {
            not_nan = false;
            return_value = -1;
          }
        }
#ifdef SMALLECHCONV
        _add_to_fixed_l_list(nb_conv, last_out, means);
        for (int ijiji=0; ijiji<contiguous_check_nb; ijiji++) {
          bg->updateSingleChannelNucleus(1,integration_method);
          count+=1;
#ifdef ITPTCTX
          for (int i=2; i<NUCLEUS_NUMBER-1; i++)
#else
            for (int i=2; i<NUCLEUS_NUMBER; i++)
#endif
            { // do not store for the cortex & CMPf
              means[i] = bg->get_single_channel_nucleus(i)->get_S();
              if (isnan(means[i])) {
                not_nan = false;
                return_value = -1;
              }
            }
          _add_to_fixed_l_list(nb_conv, last_out, means);
        }
#endif
      } else if (ch_n == 0) {
        bg->updateNucleusCells(bunch);
        count+=bunch;
        for (int i=2; i<NUCLEUS_NUMBER; i++) { // do not store for the cortex & CMPf
          means[i] = bg->get_nucleus_cells(i)->get_meanS();
        }
      } else {
        if (verbose==4) {
          bg->updateMultiChannelsNucleus_exploexplo_test(bunch,(old_count+count)*dt,0);
        } else if (verbose>4) {
          bg->updateMultiChannelsNucleus_exploexplo_test(bunch,(old_count+count)*dt,verbose-4);
        } else {
          bg->updateMultiChannelsNucleus(bunch);
        }
        count+=bunch;
#ifdef ITPTCTX
        for (int i=2; i<NUCLEUS_NUMBER-1; i++)
#else
        for (int i=2; i<NUCLEUS_NUMBER; i++)
#endif
        {
          for (int ch_i=0; ch_i<ch_n; ch_i++) {
            means[i*ch_n+ch_i] = bg->get_multi_channels_nucleus(i)->get_S(ch_i);
          }
        }
      }
#ifndef SMALLECHCONV
      _add_to_fixed_l_list(nb_conv, last_out, means);
#endif
    }

    if (verbose==1) {
      std::cout << count << " iterations done." << std::endl;
    }

#ifdef MIXEDDT
    // This is intended as a bypass of validation, as it should be done by comparing the different dt runs
    if (return_value != -1 && _conv(last_out)) {
      return_value = 1;
    }
#else
    if (ch_n == 1) { // valid on for single channel nuclei
      if (_conv(last_out) && not_nan) {
        if (do_checks) {
          bg->save_all(); // If we do test, this is the "normal" run, so we have to memorize the state
          if (verbose==1) {
            std::cout << "Apparently convergence was attained. Checking resistance to noise." << std::endl;
          }
          if (verbose==2||verbose==3) {
            std::cerr << std::endl;
          }
          for (int i=0;i<(int)(5./dt+0.5);i+=(int)(0.5/dt+0.5)) { //hammerize the circuit for 5 seconds, by steps of 0.5 seconds
            bg->hammerizeSingleChannelNucleus((int)(0.5/dt+0.5));
            if (verbose==2 || verbose==3 || verbose==4) {

              std::cerr << (old_count+count+i)*dt;
              if (ch_n == 1) {
                for (int i=0; i<NUCLEUS_NUMBER; i++) { // displays also for the cortex & CMPf
                  std::cerr << " , " << bg->get_single_channel_nucleus(i)->get_S();
                }
              } else {
                for (int i=0; i<NUCLEUS_NUMBER; i++) { // displays also for the cortex & CMPf
                  for (int ch_i=0; ch_i<ch_n; ch_i++) {
                    std::cerr << " , " << bg->get_multi_channels_nucleus(i)->get_S(ch_i);
                  }
                }
              }
              std::cerr << std::endl;
            }
          }

          if (ch_n == 1) {
            hammerization_eval = 0;
            for (int i=2; i<NUCLEUS_NUMBER; i++) { // do not store for the cortex & CMPf
              hammerization_eval += fabs(means[i] - bg->get_single_channel_nucleus(i)->get_S());
            }
            if (hammerization_eval < 5.) {
              return_value = 1;
              if (verbose==1) {
                std::cout << "Resistant to noise." << std::endl;
              }
            } else {
              if (verbose==1) {
                std::cout << "WARNING !!! NOT Resistant to noise." << std::endl;
              }
            }
          }
          } else {
            return_value = 1;
          }
        }
      }
#endif

      delete CTX;
      delete CMPf;
      if (msn_separation) {
        delete MSND1;
        delete MSND2;
      } else {
        delete MSN;
      }
      delete FSI;
      delete STN;
      delete GPe;
      delete GPi;
#ifdef ITPTCTX
      delete CTXPT;
#endif
      ////delete bg;

      return maxstn; // -1: divergence; 0: no convergence; 1: convergence
    }


    int _run_sim_tsirogiannis_2010(
        float max_duration,
        float bunch_duration,
        float dt,
        std::vector<float> &activations,
        std::vector<float> &params,
        std::vector<int> &delay,
        std::vector<float> &means,
        int verbose)
    { 
      int return_value = 0;
      float minipsilon = 1e-10;
      int corr_factor = (int)(1/(dt*1000.0)+0.5);
      int max_tau = 100+(*std::max_element(delay.begin(),delay.end()));

      float last_GPi = -10; float last_GPe = 10; float last_STN = 10000;
      float last_D1 = -10; float last_D2 = 10; float last_FS = 10000;

      std::list<std::vector<float> > last_out;

      BCBG2* bg = new BCBG2();

      bg->initialize(0,max_tau,dt);
      if (verbose==1)
        std::cout << "max_tau is " << max_tau << " & dt is " << dt << std::endl;

      float A_AMPA = 5.436563657; //2 *e
      float D_AMPA = 0.679570457; //1/4 *e

      float A_NMDA = 0.271828183; // 0.1 *e
      float D_NMDA = 0.027182818; // 1/100 *e

      float A_GABA = 5.436563657; // 2 *e
      float D_GABA = 0.453046971; // 1/6 *e

      float diameters_dummy[] = {};           int comp_nb_dummy=0;  // All nuclei here (dummy)
      float distances_dummy[] = {};                                 // All nuclei here (dummy)

      BCBG2::SingleChannelNucleus* CTXe= bg->add_single_channel_nucleus(0.,0.,0.,0.,0,distances_dummy,diameters_dummy,comp_nb_dummy,"CTXe");
      BCBG2::SingleChannelNucleus* D1  = bg->add_single_channel_nucleus(300.,2., params[TSIROGIANNIS_2010_THETA_D1_D2],0.3,0,distances_dummy,diameters_dummy,comp_nb_dummy,"D1");
      BCBG2::SingleChannelNucleus* D2  = bg->add_single_channel_nucleus(300.,2., params[TSIROGIANNIS_2010_THETA_D1_D2],0.3,0,distances_dummy,diameters_dummy,comp_nb_dummy,"D2");
      BCBG2::SingleChannelNucleus* STN = bg->add_single_channel_nucleus(500.,7.5,params[TSIROGIANNIS_2010_THETA_STN],  0.2,0,distances_dummy,diameters_dummy,comp_nb_dummy,"STN");
      BCBG2::SingleChannelNucleus* GPe = bg->add_single_channel_nucleus(400.,21.,params[TSIROGIANNIS_2010_THETA_GPe],  0.2,0,distances_dummy,diameters_dummy,comp_nb_dummy,"GPe");
      BCBG2::SingleChannelNucleus* GPi = bg->add_single_channel_nucleus(300.,16.,params[TSIROGIANNIS_2010_THETA_GPi],  0.2,0,distances_dummy,diameters_dummy,comp_nb_dummy,"GPi");

      D1->set_afferent(  A_AMPA, D_AMPA, 1, params[TSIROGIANNIS_2010_CTXe_D1_D2]*activations[TSIROGIANNIS_2010_CTXe_D1_D2], delay[TSIROGIANNIS_2010_CTXe_D1_D2],0., CTXe);
      D1->set_afferent(  A_NMDA, D_NMDA, 1, params[TSIROGIANNIS_2010_CTXe_D1_D2]*activations[TSIROGIANNIS_2010_CTXe_D1_D2], delay[TSIROGIANNIS_2010_CTXe_D1_D2],0., CTXe);

      D2->set_afferent(  A_AMPA, D_AMPA, 1, params[TSIROGIANNIS_2010_CTXe_D1_D2]*activations[TSIROGIANNIS_2010_CTXe_D1_D2],  delay[TSIROGIANNIS_2010_CTXe_D1_D2],0.,  CTXe);
      D2->set_afferent(  A_NMDA, D_NMDA, 1, params[TSIROGIANNIS_2010_CTXe_D1_D2]*activations[TSIROGIANNIS_2010_CTXe_D1_D2],  delay[TSIROGIANNIS_2010_CTXe_D1_D2],0.,  CTXe);

      STN->set_afferent( A_AMPA, D_AMPA, 1, params[TSIROGIANNIS_2010_CTXe_STN]*activations[TSIROGIANNIS_2010_CTXe_STN], delay[TSIROGIANNIS_2010_CTXe_STN],0., CTXe);
      STN->set_afferent( A_NMDA, D_NMDA, 1, params[TSIROGIANNIS_2010_CTXe_STN]*activations[TSIROGIANNIS_2010_CTXe_STN], delay[TSIROGIANNIS_2010_CTXe_STN],0., CTXe);
      STN->set_afferent( A_GABA, D_GABA, -1,params[TSIROGIANNIS_2010_GPe_STN]*activations[TSIROGIANNIS_2010_GPe_STN], delay[TSIROGIANNIS_2010_GPe_STN],0., GPe);
      STN->set_afferent( A_AMPA, D_AMPA, 1, params[TSIROGIANNIS_2010_STN_STN]*activations[TSIROGIANNIS_2010_STN_STN], delay[TSIROGIANNIS_2010_STN_STN],0., STN);
      STN->set_afferent( A_NMDA, D_NMDA, 1, params[TSIROGIANNIS_2010_STN_STN]*activations[TSIROGIANNIS_2010_STN_STN], delay[TSIROGIANNIS_2010_STN_STN],0., STN);

      GPe->set_afferent( A_AMPA, D_AMPA, 1, params[TSIROGIANNIS_2010_STN_GPe]*activations[TSIROGIANNIS_2010_STN_GPe], delay[TSIROGIANNIS_2010_STN_GPe],0., STN);
      GPe->set_afferent( A_NMDA, D_NMDA, 1, params[TSIROGIANNIS_2010_STN_GPe]*activations[TSIROGIANNIS_2010_STN_GPe], delay[TSIROGIANNIS_2010_STN_GPe],0., STN);
      GPe->set_afferent( A_GABA, D_GABA, -1,params[TSIROGIANNIS_2010_D2_GPe]*activations[TSIROGIANNIS_2010_D2_GPe],  delay[TSIROGIANNIS_2010_D2_GPe],0.,  D2); 
      GPe->set_afferent( A_GABA, D_GABA, -1,params[TSIROGIANNIS_2010_GPe_GPe]*activations[TSIROGIANNIS_2010_GPe_GPe], delay[TSIROGIANNIS_2010_GPe_GPe],0., GPe);

      GPi->set_afferent( A_AMPA, D_AMPA, 1, params[TSIROGIANNIS_2010_STN_GPi]*activations[TSIROGIANNIS_2010_STN_GPi], delay[TSIROGIANNIS_2010_STN_GPi],0., STN);
      GPi->set_afferent( A_NMDA, D_NMDA, 1, params[TSIROGIANNIS_2010_STN_GPi]*activations[TSIROGIANNIS_2010_STN_GPi], delay[TSIROGIANNIS_2010_STN_GPi],0., STN);
      GPi->set_afferent( A_GABA, D_GABA, -1,params[TSIROGIANNIS_2010_D1_GPi]*activations[TSIROGIANNIS_2010_D1_GPi],  delay[TSIROGIANNIS_2010_D1_GPi],0., D1);
      GPi->set_afferent( A_GABA, D_GABA, -1,params[TSIROGIANNIS_2010_GPe_GPi]*activations[TSIROGIANNIS_2010_GPe_GPi], delay[TSIROGIANNIS_2010_GPe_GPi],0., GPe);

      int bunch = (int)(bunch_duration/dt + 0.5);
      int count = 0;
      int old_count = 0;
      float min_duration = 2.; // simulations longer than 10s
      int nb_conv = (int)(1.0/bunch_duration+0.5); // a stabilisation over at least 50 second
      float frequency_value;
      float frequency_change_period = 0.001; // here the time in s of same cortical frequency (rounded up with regards to bunch duration)
      float frequency_change_count = frequency_change_period;

      old_count += count;
      count = 0;

      if (verbose==1)
        std::cerr << bunch << " iterations (" << bunch*dt << " s) to be done between two tests of convergence, until " << max_duration/dt << " iterations (" << max_duration << "s) are done" << std::endl;

      if (verbose==2||verbose==3)
        std::cerr << "t" << " , " << "STN" << " , " << "CTXe" << " , " << "V" << std::endl;

      while (
          (
           (
            (
             (
              ! _conv(last_out) || (verbose == 3)
             )
            )
            ||
            (count < max_tau)
            ||
            (count*dt < min_duration)
           )
           &&
           ( !isnan(GPi->get_S()) )
           &&
           ( count*dt < max_duration )
          )
          ||
          (verbose==1 && count*dt < max_duration)
          ) {
        if (verbose==2||verbose==3) {
          std::cerr << (old_count+count)*dt << " , "  << STN->get_S() << " , " << CTXe->get_S() << " , " << STN->get_V() << std::endl;
        }
        if (frequency_change_count >= frequency_change_period) {
          if (frequency_change_period < 0.) {
            bg->updateSingleChannelNucleusTsiro(bunch,-1.);
          } else {
            frequency_value = bg->updateSingleChannelNucleusTsiro(bunch);
            frequency_change_count = bunch_duration;
          }
        } else {
          frequency_change_count += bunch_duration;
          bg->updateSingleChannelNucleusTsiro(bunch,frequency_value);
        }
        count+=bunch;
        means[0] = D1->get_S();
        means[1] = D2->get_S();
        means[2] = STN->get_S();
        means[3] = GPe->get_S();
        means[4] = GPi->get_S();
        _add_to_fixed_l_list(nb_conv, last_out, means);
      }

      if (verbose==1)
        std::cout << count << " iterations done." << std::endl;

      if (_conv(last_out)) {
        if (verbose==1) {
          std::cout << "Apparently convergence was attained. Checking resistance to noise." << std::endl;
        }
        if (verbose==2||verbose==3) {
          std::cerr << std::endl;
        }
        for (int i=0;i<(int)(2./dt+0.5);i+=(int)(0.1/dt+0.5)) { //hammerize the circuit for 2 seconds, by steps of 0.1 seconds
          bg->hammerizeSingleChannelNucleusTsiro((int)(0.1/dt+0.5));
          if (verbose==2)
            std::cerr << (old_count+count+i)*dt << " , " << D1->get_S() << " , " << D2->get_S() << " , " << STN->get_S() << " , " << GPe->get_S() << " , " << GPi->get_S() << " , " << CTXe->get_S() << std::endl;
        }

        if (fabs(means[0] - D1->get_S())
            + fabs(means[1] - D2->get_S())
            + fabs(means[2] - STN->get_S())
            + fabs(means[3] - GPe->get_S())
            + fabs(means[4] - GPi->get_S())
            < 1) {
          return_value = 1;
          if (verbose==1) {
            std::cout << "Resistant to noise." << std::endl;
          }
        }
      }

      // because memory leaks suck so hard
      delete CTXe;
      delete D1;
      delete D2;
      delete STN;
      delete GPe;
      delete GPi;

      return return_value;
    }
