#include "bcbg2.hpp"

void BCBG2::MultiChannelsNucleus::initialize_multi_channels_nucleus(int ch_n, float Smax, float Sini, float vh, float k, int connectivity_type, float* dists, float* diams, int compartment_nb, const char* id)
{
  int i,j;
  // initialization of multi_channels_nucleus parameters
  this->ch_n = ch_n;
  this->Smax = Smax;
  this->vh = vh*1e-3; // convert vh from mV to V
  this->k = k*1e3; // convert vh from mV⁻¹ to V⁻¹
  this->kvh = k*vh; // kvh is dimensionless
  this->id = (char*)malloc((1+strlen(id))*sizeof(char));
  for (i=0;i<strlen(id);i++) {
    this->id[i]=id[i];
  }
  this->id[strlen(id)] = '\0';

// Computation of distance attenuation factor <
  this->Rm = 20000.*1e-4;         // HARDCODED & conversion Ohm.cm² -> Ohm.m²
  this->Ri = 200.*1e-2;           // HARDCODED & conversion Ohm.cm -> Ohm.m
  this->membrane_area = 0.;
  this->distances.resize(0); this->diameters.resize(0);
  for (i=0; i<compartment_nb; i++) {this->distances.push_back(dists[i]); this->diameters.push_back(diams[i]); this->membrane_area += 2.*3.141592654*diams[i]*dists[i];}
// Computation of distance attenuation factor >

  // initialization of counter for afferents
  n = 0;
  // allocation of arrays
  S.resize(bg_.max_tau); for(i=0;i<bg_.max_tau;i++) { S[i].assign(ch_n,Sini); }
  N_in.resize(0);
  A_in.resize(0);
  D_in.resize(0);
  C_in.resize(0);
  ADC_in.resize(0);
  DD_in.resize(0);
  D2_in.resize(0);
  dtADC_in.resize(0);
  dtDD_in.resize(0);
  dtD2_in.resize(0);
  Sign_in.resize(0);
  conn_scheme_in.resize(0);
  T_in.resize(0);
  H_in_RK4.resize(0);
  Hp_in_RK4.resize(0);
  Hpp_in_RK4.resize(0);
}

// wrapper for set_afferent with default connexion scheme "one-to-one"
void BCBG2::MultiChannelsNucleus::set_afferent(float A, float D, int Sign, float C, int T, float distance, MultiChannelsNucleus* N)
{
  set_afferent(A, D, Sign, C, T, 0., distance, N);
}

void BCBG2::MultiChannelsNucleus::set_afferent(float A, float D, int Sign, float C, int T, float connexion_scheme, float distance, MultiChannelsNucleus* N)
{
  this->initialize_new_afferent();
  A_in[n] = A*compute_distance_factor(distance)*1e-3; // convert A in V (from mV)
  D_in[n] = D*1e3; // convert D in s⁻¹ (from ms⁻¹)
  C_in[n] = C;
  ADC_in[n] = A_in[n]*D_in[n]*C;
  DD_in[n] = D_in[n]*D_in[n];
  D2_in[n] = 2.0*D_in[n];
  T_in[n] = T - bg_.max_tau;
  Sign_in[n] = Sign;
  N_in[n] = N;
  conn_scheme_in[n] = connexion_scheme;
  dtADC_in[n] = ADC_in[n]*bg_.dt;
  dtDD_in[n]  = DD_in[n]*bg_.dt;
  dtD2_in[n]  = D2_in[n]*bg_.dt;
  n++;
}

float BCBG2::MultiChannelsNucleus::compute_distance_factor(float distance)
{
  /* calculations of the distance factor */
  float end_of_pseudocompart, total_length, lambda, L, X, pos, distance_factor;
  total_length = 0.;
  for (int j=0; j<distances.size(); j++) {
    total_length += distances[j];
  }
  distance_factor = 1.;
  pos = distance*total_length; //conversion ratio -> m
  for (int i=distances.size()-1; i>=0; i--) {
    end_of_pseudocompart = 0.;
    for (int j=0; j<i; j++) {
      end_of_pseudocompart += distances[j];
    }
    if (pos > end_of_pseudocompart) {
      lambda = sqrt((diameters[i]*Rm)/(4.*Ri));
      L = (distances[i]) / lambda;
      X = (pos - end_of_pseudocompart) / lambda;
      distance_factor *= cosh(L-X) / cosh(L);
      pos = end_of_pseudocompart;
    }
  }

  return distance_factor;
}

void BCBG2::MultiChannelsNucleus::initialize_new_afferent()
{
  int i,j;
  N_in.resize(n+1);
  A_in.resize(n+1);
  D_in.resize(n+1);
  C_in.resize(n+1);
  ADC_in.resize(n+1);
  DD_in.resize(n+1);
  D2_in.resize(n+1);
  dtADC_in.resize(n+1);
  dtDD_in.resize(n+1);
  dtD2_in.resize(n+1);
  Sign_in.resize(n+1);
  conn_scheme_in.resize(n+1);
  T_in.resize(n+1);
  H_in_RK4.resize(n+1);
  H_in_RK4[n].resize(bg_.max_tau); for (i=0; i<bg_.max_tau; i++) { H_in_RK4[n][i].assign(ch_n,0); }
  Hp_in_RK4.resize(n+1);
  Hp_in_RK4[n].resize(bg_.max_tau); for (i=0; i<bg_.max_tau; i++) { Hp_in_RK4[n][i].assign(ch_n,0); }
  Hpp_in_RK4.resize(n+1);
  Hpp_in_RK4[n].resize(bg_.max_tau); for (i=0; i<bg_.max_tau; i++) { Hpp_in_RK4[n][i].assign(ch_n,0); }
}

void BCBG2::MultiChannelsNucleus::set_dt(float value)
{
	dt=value;
}

void BCBG2::MultiChannelsNucleus::update_multi_channels_nucleus_stabilize(int steps)
{
	int i,s,ch_i;
  int t = 1;
  int t_1  = 0;

  std::vector <float> virtual_S;

  float H_previous = 0.;
  float H = 0.;
  float Hp_previous = 0.;
  float Hp = 0.;

  float input_S;

  for (i=0; i<n; i++) {
    input_S = dtADC_in[i] * N_in[i]->get_S(0);
    for (s=steps; s; s--) {
      // fastest version
      Hp = Hp_previous + input_S - dtD2_in[i] * Hp_previous - dtDD_in[i] * H_previous;
      H_previous = H_previous + dt * Hp_previous;
      Hp_previous = Hp;
    }
    for (ch_i=0; ch_i<ch_n; ch_i++) {
      Hp_in_RK4[i][ch_i].assign(bg_.max_tau,Hp_previous); H_in_RK4[i][ch_i].assign(bg_.max_tau,H_previous);
    }
  }
}

void BCBG2::MultiChannelsNucleus::update_multi_channels_nucleus_evo()
{
  int t_2  = (bg_.tmod - 2 + bg_.max_tau) % bg_.max_tau;
  int t_4  = (bg_.tmod - 4 + bg_.max_tau) % bg_.max_tau;
  float dt2 = 2.*dt;
  float dt3 = 2.*dt;

  int i,ch_i,ch_j;
  int t = bg_.tmod;
  float k1,k2,k3;
  std::vector <float> sum_v;
  sum_v.resize(ch_n);
  int t_1  = (t - 1 + bg_.max_tau) % bg_.max_tau;
  float input_from_channels;
  for (ch_i=0; ch_i<ch_n; ch_i++) {
    sum_v.assign(ch_n,0);
    for (i=0; i<n; i++) {
      if (conn_scheme_in[i] == 0.) { // one to one
        input_from_channels = N_in[i]->get_S((t - T_in[i] -1 + bg_.max_tau) % bg_.max_tau, ch_i);
      } else { // one to all
        input_from_channels = 0;
        for (ch_j=0;ch_j<ch_n;ch_j++) {
          input_from_channels += N_in[i]->get_S((t - T_in[i] -1 + bg_.max_tau) % bg_.max_tau, ch_j);
        }
        input_from_channels /= ch_n;
      }
      // explicit Euler
      Hp_in_RK4[i][t][ch_i] = Hp_in_RK4[i][t_1][ch_i] +
        (
         dtADC_in[i] * input_from_channels
         - dtD2_in[i] * Hp_in_RK4[i][t_1][ch_i]
         - dtDD_in[i] * H_in_RK4[i][t_1][ch_i]
        );
      // explicit Euler
      H_in_RK4[i][t][ch_i] = H_in_RK4[i][t_1][ch_i] + dt * Hp_in_RK4[i][t_1][ch_i];

      // RK4
      //k1 = dt * (Hp_in_RK4[i][t_4][ch_i]);
      //k2 = dt * (Hp_in_RK4[i][t_2][ch_i] + dt2 * k1);
      //k3 = dt * (Hp_in_RK4[i][t_1][ch_i] + dt3 * k2);
      //H_in_RK4[i][t][ch_i] = H_in_RK4[i][t_4][ch_i] + ((8./9.) * k1 + (4./3.) * k2 + (16./9.) * k3);
    
      sum_v[ch_i] += H_in_RK4[i][t][ch_i] * Sign_in[i];
    }
  S[t][ch_i] = Smax / (1 + exp(kvh - k*sum_v[ch_i]));
  }
}

void BCBG2::MultiChannelsNucleus::display(std::ostream& os)
{
  os << "self-description of " << get_name() << " (" << ch_n << " channels nucleus)" << std::endl;
  os << "  Smax = " << Smax << std::endl;
  os << "  vh = " << vh << std::endl;
  os << "  k = " << k << std::endl;

  for (size_t i=0;i<n;i++) {
    os << "  afferent n° " << i << " is from " << N_in[i]->get_name() << std::endl;
    os << "    A = " << A_in[i] << std::endl;
    os << "    D = " << D_in[i] << std::endl;
    os << "    T = " << T_in[i] << std::endl;
    os << "    C = " << C_in[i] << std::endl;
    }
}


void BCBG2::MultiChannelsNucleus::save(MemoryMCN& mem) {
  int i,j,k;

  mem.S_backup.resize(S.size());
  for (i=0;i<S.size(); i++) {
    mem.S_backup[i].resize(S[i].size());
    for (j=0;j<S[i].size(); j++) {
      mem.S_backup[i][j] = S[i][j];
    }
  }

  mem.H_in_RK4_backup.resize(H_in_RK4.size());
  for (i=0;i<H_in_RK4.size(); i++) {
    mem.H_in_RK4_backup[i].resize(H_in_RK4[i].size());
    for (j=0;j<H_in_RK4[i].size(); j++) {
      mem.H_in_RK4_backup[i][j].resize(H_in_RK4[i][j].size());
      for (k=0;k<H_in_RK4[i][j].size(); k++) {
        mem.H_in_RK4_backup[i][j][k] = H_in_RK4[i][j][k];
      }
    }
  }

  mem.Hp_in_RK4_backup.resize(Hp_in_RK4.size());
  for (i=0;i<Hp_in_RK4.size(); i++) {
    mem.Hp_in_RK4_backup[i].resize(Hp_in_RK4[i].size());
    for (j=0;j<Hp_in_RK4[i].size(); j++) {
      mem.Hp_in_RK4_backup[i][j].resize(Hp_in_RK4[i][j].size());
      for (k=0;k<Hp_in_RK4[i][j].size(); k++) {
        mem.Hp_in_RK4_backup[i][j][k] = Hp_in_RK4[i][j][k];
      }
    }
  }
}

void BCBG2::MultiChannelsNucleus::load(MemoryMCN& mem) {
  int i,j,k;


  S.resize(mem.S_backup.size());
  for (i=0;i<mem.S_backup.size(); i++) {
    S[i].resize(mem.S_backup[i].size());
    for (j=0;j<mem.S_backup[i].size(); j++) {
      S[i][j] = mem.S_backup[i][j];
    }
  }

  H_in_RK4.resize(mem.H_in_RK4_backup.size());
  for (i=0;i<mem.H_in_RK4_backup.size(); i++) {
    H_in_RK4[i].resize(mem.H_in_RK4_backup[i].size());
    for (j=0;j<mem.H_in_RK4_backup[i].size(); j++) {
      H_in_RK4[i][j].resize(mem.H_in_RK4_backup[i][j].size());
      for (k=0;k<mem.H_in_RK4_backup[i][j].size(); k++) {
        H_in_RK4[i][j][k] = mem.H_in_RK4_backup[i][j][k];
      }
    }
  }

  Hp_in_RK4.resize(mem.Hp_in_RK4_backup.size());
  for (i=0;i<mem.Hp_in_RK4_backup.size(); i++) {
    Hp_in_RK4[i].resize(mem.Hp_in_RK4_backup[i].size());
    for (j=0;j<mem.Hp_in_RK4_backup[i].size(); j++) {
      Hp_in_RK4[i][j].resize(mem.Hp_in_RK4_backup[i][j].size());
      for (k=0;k<mem.Hp_in_RK4_backup[i][j].size(); k++) {
        Hp_in_RK4[i][j][k] = mem.Hp_in_RK4_backup[i][j][k];
      }
    }
  }
}

