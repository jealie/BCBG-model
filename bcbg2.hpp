#ifndef BCBG2_HPP
#define BCBG2_HPP

#include "constants.hpp"

#include <iostream>
#include "math.h"
#include "string.h"
#include <vector>
#include "stdio.h"
#include <stdlib.h>
#include <boost/random.hpp>

  template < typename T >
T **Allocate2DArray( int nRows, int nCols)
{
  T **ppi;
  T *pool;
  T *curPtr;
  //(step 1) allocate memory for array of elements of column

  ppi = new T*[nRows];

  //(step 2) allocate memory for array of elements of each row
  pool = new T [nRows * nCols];

  // Now point the pointers in the right place
  curPtr = pool;
  for( int i = 0; i < nRows; i++)
  {
    *(ppi + i) = curPtr;
    curPtr += nCols;
  }
  return ppi;
}

  template < typename T >
void Free2DArray(T** Array)
{
  delete [] *Array;
  delete [] Array;
}

typedef struct _MemorySCN
{
  std::vector <float> S_backup;
  std::vector <std::vector <float> > H_in_RK4_backup;
  std::vector <std::vector <float> > Hp_in_RK4_backup;
} MemorySCN;

typedef struct _MemoryMCN
{
  std::vector <std::vector <float> > S_backup;
  std::vector <std::vector <std::vector <float> > > H_in_RK4_backup;
  std::vector <std::vector <std::vector <float> > > Hp_in_RK4_backup;
} MemoryMCN;

typedef struct _MemoryBCBG2
{
  int tmod_backup;
  int previous_tmod_backup;
  std::vector <MemorySCN > scns;
  std::vector <MemoryMCN > mcns;
} MemoryBCBG2;

class BCBG2
{
  public:
#ifdef TSIROGIANNIS_2010
    BCBG2() : _uni_gaussian(base_generator_t(42u), boost::normal_distribution<>(3.0f, 1.0f)) {}
#else
    //BCBG2() : _uni_one_gaussian(base_generator_t(42u), boost::normal_distribution<>(1.0f, 0.5f)) {} // version with noise
    BCBG2() : _uni_one_gaussian(base_generator_t(42u), boost::normal_distribution<>(1.0f, 0.0f)) {} // no noise
#endif
    void initialize(int nucleus_type, int max_tau, float dt);
    void update(float dt);
    void updateSingleChannelNucleus(int steps);
    void updateSingleChannelNucleus(int steps, int integration_method);
    void updateSingleChannelNucleusWithNoise(int steps, int integration_method);
    void updateSingleChannelNucleusTsiro(int steps, float value);
    float updateSingleChannelNucleusTsiro(int steps); // returns the random cortical firing rate used
    void updateMultiChannelsNucleus(int steps);
    void updateMultiChannelsNucleus_basic_test(int steps, float time);
    void updateMultiChannelsNucleus_humphries_test(int steps, float time);
    float scalevar(float m);
    float scalevarreduc(float m, float reduc);
    void updateMultiChannelsNucleus_exploexplo_test(int steps, float time, int i);
    void hammerizeSingleChannelNucleus(int steps);
    void hammerizeSingleChannelNucleusTsiro(int steps);
    void updateSingleChannelNucleusVana(float input, int bunch);
    void updateSingleChannelNucleusDebug(float time_of_peak_in_s, float dt, float nb_of_s);
    void testPSPSingleChannelNucleus(int pulse, int steps);
    void stabilize_all(int steps);
    bool testTimingResponse();

    float basic_test(int ch, float time);
    float exploexplo_test(int ch, float time, int i);
    float exploexplo_testbis(int ch, float time, float debut_time, float end_time);
    float visualization_humphries_et_al_2006(int ch, float time);

    void updateNucleusCells(int steps);
    void autopiloteNucleusCells(int steps, float* values);

    void save_all();
    void load_all();
    void save_all(MemoryBCBG2& mem);
    void load_all(MemoryBCBG2& mem);

    class SingleChannelNucleus;
    friend class SingleChannelNucleus;
    class MultiChannelsNucleus;
    friend class MultiChannelsNucleus;
    class NucleusCells;
    friend class NucleusCells;

    SingleChannelNucleus* add_single_channel_nucleus(float Smax, float Sini, float vh, float k, int damping, float* dists, float* diams, int compartment_nb, const char* id);
    MultiChannelsNucleus* add_multi_channels_nucleus(int channels_number, float Smax, float Sini, float vh, float k, int connectivity_type, float* dists, float* diams, int compartment_nb, const char* id);
    NucleusCells* add_nucleus_cells(int cells_number, float Smax, float Sini, float vh, float k, float V_rest, float R, float Rm, float Ri, float* dists, float* diams, int compartment_nb, const char* id);

    SingleChannelNucleus* get_single_channel_nucleus(int i) {return SCN[i];};
    MultiChannelsNucleus* get_multi_channels_nucleus(int i) {return MCN[i];};
    NucleusCells* get_nucleus_cells(int i) {return NC[i];};

  protected:
    // general variables of the circuit
    int tmod;
    int previous_tmod;
    int max_tau;
    int nucleus_type;
    float dt;
    float treal;

    // list of nuclei
    std::vector <SingleChannelNucleus *> SCN;
    std::vector <MultiChannelsNucleus *> MCN;
    std::vector <NucleusCells *> NC;
    int n_nuclei; // counter used in the addition of nuclei

    // random generator
    typedef boost::mt19937 base_generator_t;
    typedef boost::variate_generator<base_generator_t, boost::normal_distribution<> > rand_t;
    rand_t _uni_one_gaussian;

    // for use in store/restore function
    int tmod_backup;
    int previous_tmod_backup;
};



class Constants {
  public:
    // defining some constants
    static const float A_AMPA = 0.001 * 5.43656365692;        // 2*e
    static const float D_AMPA = 0.65238763883017081;  // (1/4)*e
    static const float A_NMDA = 0.001 * 0.27182818284590454;  // 0.1*e
    static const float D_NMDA = 0.027182818284590453; // (1/100)*e
    static const float A_GABA = 0.001 * 5.4365636569180902;   // 2*e
    static const float D_GABA = 0.45304697140984085;  // (1/6)*e
};


class BCBG2::MultiChannelsNucleus
{
  public:
    friend class BCBG2;
    MultiChannelsNucleus(BCBG2& bg) : bg_(bg) {}
    void initialize_multi_channels_nucleus(int ch_n, float Smax, float Sini, float vh, float k, int connectivity_type, float* dists, float* diams, int compartment_nb, const char* id);
    void update_multi_channels_nucleus_stabilize(int steps);
    void update_multi_channels_nucleus_evo();
    void set_afferent(float A, float D, int Sign, float C, int T, float distance, MultiChannelsNucleus* N);
    void set_afferent(float A, float D, int Sign, float C, int T, float connexion_scheme, float distance, MultiChannelsNucleus* N);
    float compute_distance_factor(float distance);
    void initialize_new_afferent();
    void display(std::ostream&);
    void set_dt(float value);

    // various getter/setter
    float get_S(int tau_in, int ch_n) { return S[tau_in][ch_n]; }
    // Warning ! get_S() returns S for _last step done_ not for  _current step_
    float get_S(int ch_n) { return S[(bg_.tmod - 1 + bg_.max_tau) % bg_.max_tau][ch_n]; } 
    void set_S(float value, int tau_in, int ch_n) { S[tau_in][ch_n] = value; }
    void set_S(float value, int ch_n) { S[bg_.tmod][ch_n] = value; }
    int get_afferents_number() { return n; } 
    int get_channels_number() { return ch_n; } 
    char* get_name() { return id; }

    void save(MemoryMCN& mem);
    void load(MemoryMCN& mem);


  protected:
    // to access the constants of the circuit
    BCBG2& bg_;

    // general variables of the multi_channels_nucleus
    int ch_n; // number of concurrent channels
    std::vector <std::vector <float> > S; //output frequency
    float Smax;
    float vh;
    float k;
    float kvh;
    float dt;
    char *id;

    float Rm;                      // Membrane resistance
    float Ri;                      // Internal resistivity
    float membrane_area;           // Surface of the neuron
    std::vector <float> distances; // used to reconstruct a simple model of dendrites
    std::vector <float> diameters; // used to reconstruct a simple model of dendrites

    // afferents
    float k1,k2,k3,k4,k5,k6,k7;
    std::vector <std::vector <float> > H_in;
    std::vector <std::vector <std::vector <float> > > H_in_RK4;
    std::vector <std::vector <std::vector <float> > > Hp_in_RK4;
    std::vector <std::vector <std::vector <float> > > Hpp_in_RK4;
    std::vector <float> A_in;
    std::vector <float> D_in;
    std::vector <float> C_in;
    std::vector <float> ADC_in;
    std::vector <float> DD_in;
    std::vector <float> D2_in;
    std::vector <float> dtADC_in;
    std::vector <float> dtDD_in;
    std::vector <float> dtD2_in;
    std::vector <int> T_in;
    std::vector <int> Sign_in;
    std::vector <float> conn_scheme_in;
    std::vector <MultiChannelsNucleus *> N_in;
    int n; // counter used in the addition of afferents

    // to store/restore the state of a circuit
    std::vector <std::vector <float> > S_backup; //output frequency
    std::vector <std::vector <std::vector <float> > > H_in_RK4_backup;
    std::vector <std::vector <std::vector <float> > > Hp_in_RK4_backup;
};

class BCBG2::NucleusCells
{
  // this class implements the cell-level approximation
  // NOT TO BE USED - THIS IS NOT FINISHED
  // Use instead the classes SingleChannelNucleus (single channel allows quick simulation) and MultiChannelsNucleus (the multi channels version)
  public:

    struct Target {
      BCBG2::NucleusCells * target_nucleus;
      int cell_number;
      int ion_channel_number;
      float T; 
      float distance_factor;
    };
    
    struct IonChannel {
      float A;           // | 
      float D;           // | 
      float Vrev;        // | 
      bool is_nmda;      // | these informations are duplicated for each cell
      float time_to_add; // |
      float to_add;      // |
    };

    struct Spike {
      int target_number;
      float t;
      float t_expiration; // set to the same value as IonChannel
    };
    

    friend class BCBG2;
    NucleusCells(BCBG2& bg) : bg_(bg), is_firing(rng) {}

    void initialize_cells(int cells_n, float Smax, float Sini, float V_threshold, float k, float V_rest, float R, float Rm, float Ri, float* dists, float* diams, int compartment_nb, const char* id, int nucleus_number);
    int affect_ion_channel(float A, float D, float Vrev, bool is_nmda);
    void add_efferent_axon(NucleusCells* nc, float cs, int T, float C, float distance_factor, int ion_channel_n);
    void add_afferent_axon(NucleusCells* nc, float cs, int T, float C, float distance_factor, int ion_channel_n);
    void add_afferent_synapse(NucleusCells* source_nucleus, float cs, int T, float C, float distance_factor, int ion_channel_n);
    bool set_afferent(float A, float D, int S, float C, int T, float cs, float distance, NucleusCells* N);
    bool add(float A, float D, float Vrev, float C, int T, float cs, float distance, NucleusCells* N);
    bool nmda(float A, float D, float Vrev, float C, int T, float cs, float distance, NucleusCells* N);
    bool add_ampa_gaba_nmda(float A, float D, float Vrev, float C, int T, float cs, float distance, NucleusCells* N, bool is_nmda);
    void set_dt(float value);
    void fast_add_g(int cell_i, int ion_channel_i, float distance_factor, float t_sent);
    void add_g(int cell_i, int ion_channel_i, float distance_factor, float t_sent);
    void update_cells();
    void set_rate(float rate, bool reset);
    void set_rate(float rate); // wrapper around previous function, with reset=true
    void prepare_next_iteration();
    void display(std::ostream& os);
    void display_summary(std::ostream& os);

    // various getter/setter
    int get_nucleus_number() {
      return this->nucleus_number;
    }

    float get_S(int cell_i) { float hz=0; for (int i=0; i<times_of_spikes[cell_i].size(); i++) {if (bg_.treal - times_of_spikes[cell_i][i] < 1.) {hz+=1.;} } return (hz); } 
    float get_meanS()
    {
      float hz=0; 
      for (int cell_i=0; cell_i<cells_n; cell_i++) { 
        for (int i=0; i<times_of_spikes[cell_i].size(); i++) {
          if (bg_.treal - times_of_spikes[cell_i][i] < 1.) {
            hz+=1.;
          }
        }
      }
      return (hz/((float)cells_n));
    }
    void debugg()
    { 
      float now = bg_.treal; 
      std::cout << "NOW=" << now << std::endl;
      std::cout << "DT=" << bg_.dt << std::endl;
      for (int cell_i=0; cell_i<cells_n; cell_i++) {
        std::cout << cell_i << " : " << last_spiked[cell_i] << std::endl;
      }
    } 
    float get_instantaneous_S3()
    {
      float hz=0; 
      for (int cell_i=0; cell_i<cells_n; cell_i++) { 
        for (int i=0; i<times_of_spikes[cell_i].size(); i++) {
          if (bg_.treal - times_of_spikes[cell_i][i] < 0.1) {
            hz+=1.;
          }
        }
      }
      return (10.*hz/((float)cells_n));
    }
    float get_instantaneous_S2()
    { 
      float now = bg_.treal; 
      int nb_discharging = 0;
      for (int cell_i=0; cell_i<cells_n; cell_i++) {
        if (now - last_spiked[cell_i] < absolute_refractory_period-dt) { // refractory period dependent of maximum rate
          nb_discharging++;
        } 
      }
      return (((float)nb_discharging) / ((float)cells_n)) * (absolute_refractory_period/dt);
    } 
    float get_instantaneous_S()
    { 
      float now = bg_.treal; 
      float isi=1./Smax; 
      float no_spiking_cells = 0.;
      for (int cell_i=0; cell_i<cells_n; cell_i++) {
        if (now - last_spiked[cell_i] >= 1./Smax + 2.*dt) { 
          isi += now - last_spiked[cell_i]; 
          no_spiking_cells = no_spiking_cells + 1.;
        } 
      }
      if (no_spiking_cells == 0.) {
        return 0.;
      }
      return (no_spiking_cells)/isi; 
    } 
    float get_meanV() { float v=0; for (int cell_i=0; cell_i<cells_n; cell_i++) { v+=V[cell_i];} return (v/((float)cells_n)); } // TODO TODO TODO TODO TODO original
    float get_meanV_no_silenced() { float v=0; for (int cell_i=0; cell_i<cells_n; cell_i++) { if (V[cell_i] != V_rest) {v+=V[cell_i];}} return (v/((float)cells_n)); } // corrected to exclude silenced
    float get_meanG() { float gg=0; for (int cell_i=0; cell_i<cells_n; cell_i++) { for (int ion_channel_i=0; ion_channel_i<ion_channels.size(); ion_channel_i++) {gg+=g[cell_i][ion_channel_i];}} return (gg/((float)cells_n*((float)ion_channels.size()))); } 
    float get_mean_inc_sum_v() { float mean_inc_sum_v=0; for (int cell_i=0; cell_i<cells_n; cell_i++) { mean_inc_sum_v+=inc_sum_v[cell_i];} return (mean_inc_sum_v/((float)cells_n)); } 
    float get_V(int cells_i) { return V[cells_i]; } 
    int get_afferents_number() { return n; } 
    char* get_name() { return id; }
    int get_cells_number() {return cells_n; }
    float get_ion_channel_time_to_add(int ion_channel_i) {return ion_channels[ion_channel_i].time_to_add; } // TODO this is dirty, all cells should share ion_channels (see struct IonChannel)

    void save();
    void restore();

  protected:
    // to access the constants of the circuit
    BCBG2& bg_;
    boost::mt19937 rng;
    boost::uniform_01<boost::mt19937> is_firing;

    // general variables of the nucleus_cells
    int cells_n;                   // number of cells
    std::vector <std::vector <float> > S; //output frequency
    float Smax;
    std::vector <float> V_threshold;
    std::vector <float> V_reset;                 // Reset potential
    std::vector <float> V;         // Mean voltage of the neurons
    std::vector <float> inc_sum_v;     // Mean input voltage of the neurons from ionic channels
    std::vector <std::vector <float> > g;
    std::vector <std::vector <float> > next_g;
    float V_rest;                  // Mean resting voltage
    float R;                       // Membrane input resistance
    float Rm;                      // Membrane resistance
    float Ri;                      // Internal resistivity
    float membrane_area;           // Surface of the neuron
    std::vector <float> distances; // used to reconstruct a simple model of dendrites
    std::vector <float> diameters; // used to reconstruct a simple model of dendrites
    float absolute_refractory_period;
    float time_to_vanish;
    float dt;
    char *id;
    int nucleus_number;

    // afferents
    float k1,k2,k3,k4,k5,k6,k7;
    std::vector <float> Vrev_in;            // Reversal Potential associated with synapse
    std::vector <float> Distance_factor_in; // Distance factor computed as a single or multi compartment with sealed end
    std::vector <NucleusCells *> N_in;
    int n; // counter used in the addition of afferents

    std::vector <float> last_spiked; // TODO
    std::vector <struct IonChannel> ion_channels; // TODO
    std::vector <std::vector <struct Target> > targets; // TODO
    std::vector <std::vector <struct Spike> > spikes; // TODO
    std::vector <std::vector <float> > times_of_spikes; // TODO
};


class BCBG2::SingleChannelNucleus
{
  public:
    friend class BCBG2;
    SingleChannelNucleus(BCBG2& bg) : bg_(bg) {}
    void initialize_single_channel_nucleus(float Smax, float Sini, float vh, float k, int damping, float* dists, float* diams, int compartment_nb, const char* id);
    void update_single_channel_nucleus_stabilize(int steps);
    void update_single_channel_nucleus_evo();
    float get_all_inputs()  // for debugging purpose
    {
      float r = 0;
      r += A_in[0] * C_in[0] * N_in[0]->get_S((bg_.tmod - T_in[0] - 1 + bg_.max_tau ) % bg_.max_tau);
      r += A_in[2] * C_in[2] * N_in[2]->get_S((bg_.tmod - T_in[2] - 1 + bg_.max_tau ) % bg_.max_tau);
      return r;      
    }
    void update_single_channel_nucleus_euler();
    void update_single_channel_nucleus_rk3();
    void update_single_channel_nucleus_evo_Hp();
    void update_single_channel_nucleus_vana();
    void update_single_channel_nucleus_damping_vana();
    void update_single_channel_nucleus();
    void set_afferent(float A, float D, int Sign, float C, int T, float distance, SingleChannelNucleus* N);
    void set_afferent(float A, float D, int Sign, float C, int T, float to_be_ignored, float distance, SingleChannelNucleus* N); // Dummy function, calls other "set_afferent"
    void set_afferent(float nu, int T, SingleChannelNucleus* N);
    float compute_distance_factor(float distance);
    void initialize_new_afferent();
    void display(std::ostream&);
    void display_Hs() {for (int i=0;i<n;i++) std::cout << "H=" << H_in_RK4[i][0] << " Hp=" << Hp_in_RK4[i][0] << std::endl; std::cout <<"-------------" << std::endl;}; // for debug
    void set_dt(float value);

    void save();
    void load();
    void save(MemorySCN& mem);
    void load(MemorySCN& mem);

    // various getter/setter
    float get_S(int tau_in) { if (is_damped) {return Phi_in_RK4[tau_in];} else {return S[tau_in];} }
    // Warning ! get_S() returns S for _last step done_ not for  _current step_
    float get_S() { if (is_damped) {return Phi_in_RK4[(bg_.tmod - 1 + bg_.max_tau) % bg_.max_tau];} else {return S[(bg_.tmod - 1 + bg_.max_tau) % bg_.max_tau];} } 
    void set_S(float value, int tau_in) { S[tau_in] = value; }
    void set_S(float value) { S[bg_.tmod] = value; }
    float get_H(int i) { return H_in_RK4[i][(bg_.tmod - 1 + bg_.max_tau) % bg_.max_tau]; } 
    float get_V() {float V=0; for (int i=0;i<n;i++) V+= H_in_RK4[i][(bg_.tmod  + bg_.max_tau) % bg_.max_tau] * Sign_in[i]; return V; } 
    float get_Hp(int i) { return Hp_in_RK4[i][(bg_.tmod - 1 + bg_.max_tau) % bg_.max_tau]; } 
    void scale_Hp(int i, float value) { Hp_in_RK4[i][(bg_.tmod - 1 + bg_.max_tau) % bg_.max_tau] *= value; } 
    int get_afferents_number() { return n; } 
    char* get_name() { return id; }


  protected:
    // to access the constants of the circuit
    BCBG2& bg_;

    // general variables of the multi_channels_nucleus
    std::vector <float> S; //output frequency
    float Smax;
    float vh;
    float k;
    float kvh;
    float dt;
    char *id;
    int is_damped;

    float Rm;                      // Membrane resistance
    float Ri;                      // Internal resistivity
    float membrane_area;           // Surface of the neuron
    std::vector <float> distances; // used to reconstruct a simple model of dendrites
    std::vector <float> diameters; // used to reconstruct a simple model of dendrites

    // afferents
    std::vector <float> H_in;
    float k1,k2,k3,k4,k5,k6,k7;
    float k1bis,k2bis,k3bis;
    std::vector <float> Phi_in_RK4;
    std::vector <float> Phip_in_RK4;
    std::vector <std::vector <float> > H_in_RK4;
    std::vector <std::vector <float> > Hp_in_RK4;
    std::vector <std::vector <float> > Hpp_in_RK4;
    std::vector <float> H_in_old;
    std::vector <float> Hp_in;
    std::vector <float> Hpp_in;
    std::vector <float> Hp_in_old;
    std::vector <float> Hpp_in_old;
    std::vector <float> A_in;
    std::vector <float> D_in;
    std::vector <float> C_in;
    std::vector <float> ADC_in;
    std::vector <float> DD_in;
    std::vector <float> D2_in;
    std::vector <float> nu_in;
    std::vector <float> dtADC_in;
    std::vector <float> dtDD_in;
    std::vector <float> dtD2_in;
    std::vector <int> T_in;
    std::vector <int> Sign_in;
    std::vector <SingleChannelNucleus *> N_in;
    int n; // counter used in the addition of afferents

    // to store/restore the state of a circuit
    std::vector <float> S_backup; //output frequency
    std::vector <std::vector <float> > H_in_RK4_backup;
    std::vector <std::vector <float> > Hp_in_RK4_backup;
};


#endif
