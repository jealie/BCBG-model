#include "bcbg2.hpp"
#include "constants.hpp"


void BCBG2::SingleChannelNucleus::initialize_single_channel_nucleus(float Smax, float Sini, float vh, float k, int damping, float* dists, float* diams, int compartment_nb, const char* id)
{
  int i,j;
  // initialization of single_channel_nucleus parameters
  this->Smax = Smax; // in Hz
  //this->vh = vh;
  this->vh = vh*1e-3; // convert vh from mV to V
  //this->k = k;
  this->k = k*1e3; // convert vh from mV⁻¹ to V⁻¹
  this->kvh = k*vh; // kvh is dimensionless
  this->id = (char*)malloc((1+strlen(id))*sizeof(char)); 
  for (i=0;i<strlen(id);i++) {
    this->id[i]=id[i];
  }
  this->id[strlen(id)] = '\0';

  this->is_damped = damping;

// Computation of distance attenuation factor <
  this->Rm = 20000.*1e-4;         // conversion Ohm.cm² -> Ohm.m²
  this->Ri = 200.*1e-2;           // conversion Ohm.cm -> Ohm.m
  this->membrane_area = 0.;
  this->distances.resize(0); this->diameters.resize(0);
  for (i=0; i<compartment_nb; i++) {this->distances.push_back(dists[i]); this->diameters.push_back(diams[i]); this->membrane_area += 2.*3.141592654*diams[i]*dists[i];}
// Computation of distance attenuation factor >

  //nb_afferents = afferents_number;
  // initialization of counter for afferents
  n = 0;
  // allocation of arrays
  S.resize(bg_.max_tau);
  N_in.resize(0);
  A_in.resize(0);
  D_in.resize(0);
  C_in.resize(0);
  ADC_in.resize(0);
  DD_in.resize(0);
  D2_in.resize(0);
  nu_in.resize(0);
  dtADC_in.resize(0);
  dtDD_in.resize(0);
  dtD2_in.resize(0);
  Sign_in.resize(0);
  T_in.resize(0);
  H_in_RK4.resize(0);
  Hp_in_RK4.resize(0);
  // initialization to Sini for the initial discharge rate
  for(i=0;i<bg_.max_tau;i++) {
      S[i]=Sini;
  }
}

void BCBG2::SingleChannelNucleus::set_afferent(float A, float D, int Sign, float C, int T, float to_be_ignored, float distance, SingleChannelNucleus* N)
{
 SingleChannelNucleus::set_afferent(A, D, Sign, C, T, distance, N); // Dummy function, calls other "set_afferent"
}
void BCBG2::SingleChannelNucleus::set_afferent(float A, float D, int Sign, float C, int T, float distance, SingleChannelNucleus* N)
{
  this->initialize_new_afferent();
  A_in[n] = A*compute_distance_factor(distance)*1e-3; // convert A in V (from mV)
  D_in[n] = D*1e3; // convert D in s⁻¹ (from ms⁻¹)
  C_in[n] = C;
  nu_in[n] = 0.;
  ADC_in[n] = A_in[n]*D_in[n]*C;
  DD_in[n] = D_in[n]*D_in[n];
  D2_in[n] = 2.*D_in[n];
  T_in[n] = T - bg_.max_tau;
  Sign_in[n] = Sign;
  N_in[n] = N;
  dtADC_in[n] = ADC_in[n]*bg_.dt;
  dtDD_in[n]  = DD_in[n]*bg_.dt;
  dtD2_in[n]  = D2_in[n]*bg_.dt;
  n++;
}

void BCBG2::SingleChannelNucleus::set_afferent(float nu, int T, SingleChannelNucleus* N)
{
  this->initialize_new_afferent();
  A_in[n] = 0;
  D_in[n] = 0;
  C_in[n] = 0;
  nu_in[n] = nu;
  ADC_in[n] = 0;
  DD_in[n] = 0;
  D2_in[n] = 0;
  Sign_in[n] = 0;
  T_in[n] = T - bg_.max_tau;
  N_in[n] = N;
  n++;
}

float BCBG2::SingleChannelNucleus::compute_distance_factor(float distance)
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

void BCBG2::SingleChannelNucleus::initialize_new_afferent()
{
  N_in.resize(n+1);
  A_in.resize(n+1);
  D_in.resize(n+1);
  C_in.resize(n+1);
  ADC_in.resize(n+1);
  DD_in.resize(n+1);
  D2_in.resize(n+1);
  nu_in.resize(n+1);
  dtADC_in.resize(n+1);
  dtDD_in.resize(n+1);
  dtD2_in.resize(n+1);
  Sign_in.resize(n+1);
  T_in.resize(n+1);
  H_in_RK4.resize(n+1); H_in_RK4[n].assign(bg_.max_tau,0);
  Hp_in_RK4.resize(n+1); Hp_in_RK4[n].assign(bg_.max_tau,0);
}

void BCBG2::SingleChannelNucleus::set_dt(float value)
{
	dt=value;
}

void BCBG2::SingleChannelNucleus::update_single_channel_nucleus()
{
	int i;
	float sum_v = 0;
	for (i=0; i<n; i++) {
		Hp_in[i] = Hp_in_old[i] + dtADC_in[i] * N_in[i]->get_S((bg_.tmod - T_in[i]) % bg_.max_tau) - dtD2_in[i] * Hp_in_old[i] - dtDD_in[i] * H_in_old[i];
		H_in[i] = H_in_old[i] + Hp_in[i] * dt;
		H_in_old[i] = H_in[i];
		Hp_in_old[i] = Hp_in[i];
		sum_v += H_in[i] * Sign_in[i];
	}
	S[(bg_.tmod+1) % bg_.max_tau] = Smax / (1 + exp(kvh - k*sum_v));
}

void BCBG2::SingleChannelNucleus::update_single_channel_nucleus_stabilize(int steps)
{
	int i,s;
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

      ////// explicit Euler
      //Hp_in_RK4[i][t] = Hp_in_RK4[i][t_1] +
      //  (
      //   virtual_S[i]
      //   - dtD2_in[i] * Hp_in_RK4[i][t_1]
      //   - dtDD_in[i] * H_in_RK4[i][t_1]
      //  );

      ////// explicit Euler
      //H_in_RK4[i][t] = H_in_RK4[i][t_1] + dt * Hp_in_RK4[i][t_1];
      //std::cout << "just added " << dt * Hp_in_RK4[i][t_1] << " but Hp was " << Hp_in_RK4[i][t_1] << std::endl;
      //t_1++;
      //t++;
      //if (t_1 == bg_.max_tau) {
      //  t_1 = 0;
      //}
      //if (t == bg_.max_tau) {
      //  t = 0;
      //}
    }
    Hp_in_RK4[i].assign(bg_.max_tau,Hp_previous); H_in_RK4[i].assign(bg_.max_tau,H_previous);
  }
}

void BCBG2::SingleChannelNucleus::update_single_channel_nucleus_evo_Hp()
{
	int i;
	float sum_v = 0;
  int t = bg_.tmod;
  int t_1  = (t - 1 + bg_.max_tau) % bg_.max_tau;
  int t_2  = (t - 2 + bg_.max_tau) % bg_.max_tau;
  int t_3  = (t - 3 + bg_.max_tau) % bg_.max_tau;
  int t_4  = (t - 4 + bg_.max_tau) % bg_.max_tau;
	for (i=0; i<n; i++) {
    k1 = dt * (
        ADC_in[i] * N_in[i]->get_S((t-4 - T_in[i] + bg_.max_tau) % bg_.max_tau)
        - D2_in[i] * (Hp_in_RK4[i][t_4])
        - DD_in[i] * H_in_RK4[i][t_4]
        );
    k2 = dt * (
        ADC_in[i] * N_in[i]->get_S((t-2 - T_in[i] + bg_.max_tau) % bg_.max_tau)
        - D2_in[i] * (Hp_in_RK4[i][t_4]+2.*k1)
        - DD_in[i] * H_in_RK4[i][t_2]
        );
    k3 = dt * (
        ADC_in[i] * N_in[i]->get_S((t-1 - T_in[i] + bg_.max_tau) % bg_.max_tau)
        - D2_in[i] * (Hp_in_RK4[i][t_4]+3.*k2)
        - DD_in[i] * H_in_RK4[i][t_1]
        );
    Hp_in_RK4[i][t] = Hp_in_RK4[i][t_4] + (8./9.) * k1 + (4./3.) * k2 + (16./9.) * k3;

    H_in_RK4[i][t] = H_in_RK4[i][t_1];

		sum_v += H_in_RK4[i][t] * Sign_in[i];
	}
	S[t] = Smax / (1 + exp(kvh - k*sum_v));
}

void BCBG2::SingleChannelNucleus::update_single_channel_nucleus_euler()
{
	int i;
	float sum_v = 0;
  int t = bg_.tmod;
  int t_1  = (t - 1 + bg_.max_tau) % bg_.max_tau;
  
	for (i=0; i<n; i++) {
    // This gives the same discharge rates as Tsirogiannis et al 2010
    ////// explicit Euler
    Hp_in_RK4[i][t] = Hp_in_RK4[i][t_1] +
      (
       dtADC_in[i] * N_in[i]->get_S((t - T_in[i] - 1 + bg_.max_tau ) % bg_.max_tau)
       - dtD2_in[i] * Hp_in_RK4[i][t_1]
       - dtDD_in[i] * H_in_RK4[i][t_1]
      );
    ////// explicit Euler
    H_in_RK4[i][t] = H_in_RK4[i][t_1] + dt * Hp_in_RK4[i][t_1];

		sum_v += H_in_RK4[i][t] * Sign_in[i];
	}
	S[t] = Smax / (1 + exp(kvh - k*sum_v));
}

void BCBG2::SingleChannelNucleus::update_single_channel_nucleus_rk3()
{
	int i;
	float sum_v = 0;
  int t = bg_.tmod;
  int t_1  = (t - 1 + bg_.max_tau) % bg_.max_tau;
  int t_2  = (t - 2 + bg_.max_tau) % bg_.max_tau;
  int t_4  = (t - 4 + bg_.max_tau) % bg_.max_tau;

	for (i=0; i<n; i++) {
    k1 = dt * (
        ADC_in[i] * N_in[i]->get_S((t-4 - T_in[i] + bg_.max_tau) % bg_.max_tau)
        - D2_in[i] * (Hp_in_RK4[i][t_4])
        - DD_in[i] * H_in_RK4[i][t_4]
        );
    k1bis = dt * (Hp_in_RK4[i][t_4]);
    k2 = dt * (
        ADC_in[i] * N_in[i]->get_S((t-2 - T_in[i] + bg_.max_tau) % bg_.max_tau)
        - D2_in[i] * (Hp_in_RK4[i][t_4]+2.*k1)
        - DD_in[i] * (H_in_RK4[i][t_4]+2.*k1bis)
        );
    k2bis = dt * (Hp_in_RK4[i][t_4]+2.*k1bis);
    k3 = dt * (
        ADC_in[i] * N_in[i]->get_S((t-1 - T_in[i] + bg_.max_tau) % bg_.max_tau)
        - D2_in[i] * (Hp_in_RK4[i][t_4]+3.*k2)
        - DD_in[i] * (H_in_RK4[i][t_4]+3.*k2bis)
        );
    k3bis = dt * (Hp_in_RK4[i][t_4]+3.*k2bis);
    Hp_in_RK4[i][t] = Hp_in_RK4[i][t_4] + (8./9.) * k1 + (4./3.) * k2 + (16./9.) * k3;
    H_in_RK4[i][t] = H_in_RK4[i][t_4] + (8./9.) * k1bis + (4./3.) * k2bis + (16./9.) * k3bis;

		sum_v += H_in_RK4[i][t] * Sign_in[i];
	}
	S[t] = Smax / (1 + exp(kvh - k*sum_v));
}


void BCBG2::SingleChannelNucleus::update_single_channel_nucleus_evo()
{
	int i;
	float sum_v = 0;
  int t = bg_.tmod;
  int t_1  = (t - 1 + bg_.max_tau) % bg_.max_tau;
  int t_2  = (t - 2 + bg_.max_tau) % bg_.max_tau;
  int t_4  = (t - 4 + bg_.max_tau) % bg_.max_tau;

	for (i=0; i<n; i++) {

#ifdef INTEGRATIONEULER
    // This gives the same discharge rates as Tsirogiannis et al 2010
    ////// explicit Euler
    Hp_in_RK4[i][t] = Hp_in_RK4[i][t_1] +
      (
       dtADC_in[i] * N_in[i]->get_S((t - T_in[i] - 1 + bg_.max_tau ) % bg_.max_tau)
       - dtD2_in[i] * Hp_in_RK4[i][t_1]
       - dtDD_in[i] * H_in_RK4[i][t_1]
      );
    ////// explicit Euler
    H_in_RK4[i][t] = H_in_RK4[i][t_1] + dt * Hp_in_RK4[i][t_1];
#endif

    //// This gives the same discharge rates as Tsirogiannis et al 2010
    //////// explicit Euler
    //H_in_RK4[i][t] = H_in_RK4[i][t_1] + dt * Hp_in_RK4[i][t_1];
    //////// explicit (?) Euler
    //Hp_in_RK4[i][t] = Hp_in_RK4[i][t_1] +
    //  (
    //   dtADC_in[i] * N_in[i]->get_S(t_1)
    //   - dtD2_in[i] * Hp_in_RK4[i][t_1]
    //   - dtDD_in[i] * H_in_RK4[i][t_1]
    //  );
    
    //// midpoint rule
    //k1 = dt * (
    //    ADC_in[i] * N_in[i]->get_S((t-2 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_2])
    //    - DD_in[i] * H_in_RK4[i][t_2]
    //    );
    //k2 = dt * (
    //    ADC_in[i] * N_in[i]->get_S((t-1 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_2]+k1)
    //    - DD_in[i] * H_in_RK4[i][t_1]
    //    );
    //Hp_in_RK4[i][t] = Hp_in_RK4[i][t_2] + 2. * k2;
    //k1 = dt * (Hp_in_RK4[i][t_2]);
    //k2 = dt * (Hp_in_RK4[i][t_1]);
    //H_in_RK4[i][t] = H_in_RK4[i][t_2] + 2. * k2;

    //// Bogacki-Shampine
    //// work on a simple example, not on the whole circuit
    //k1 = dt * (
    //    ADC_in[i] * N_in[i]->get_S((t-4 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_4])
    //    - DD_in[i] * H_in_RK4[i][t_4]
    //    );
    //k2 = dt * (
    //    ADC_in[i] * N_in[i]->get_S((t-2 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_4]+2.*k1)
    //    - DD_in[i] * H_in_RK4[i][t_2]
    //    );
    //k3 = dt * (
    //    ADC_in[i] * N_in[i]->get_S((t-1 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_4]+3.*k2)
    //    - DD_in[i] * H_in_RK4[i][t_1]
    //    );
    //Hp_in_RK4[i][t] = Hp_in_RK4[i][t_4] + (8./9.) * k1 + (4./3.) * k2 + (16./9.) * k3;
    //k1 = dt * (Hp_in_RK4[i][t_4]);
    //k2 = dt * (Hp_in_RK4[i][t_2]);
    //k3 = dt * (Hp_in_RK4[i][t_1]);
    //H_in_RK4[i][t] = H_in_RK4[i][t_4] + (8./9.) * k1 + (4./3.) * k2 + (16./9.) * k3;

#ifdef INTEGRATIONRK3
    // Bogacki-Shampine
    k1 = dt * (
        ADC_in[i] * N_in[i]->get_S((t-4 - T_in[i] + bg_.max_tau) % bg_.max_tau)
        - D2_in[i] * (Hp_in_RK4[i][t_4])
        - DD_in[i] * H_in_RK4[i][t_4]
        );
    k1bis = dt * (Hp_in_RK4[i][t_4]);
    k2 = dt * (
        ADC_in[i] * N_in[i]->get_S((t-2 - T_in[i] + bg_.max_tau) % bg_.max_tau)
        - D2_in[i] * (Hp_in_RK4[i][t_4]+2.*k1)
        - DD_in[i] * (H_in_RK4[i][t_4]+2.*k1bis)
        );
    k2bis = dt * (Hp_in_RK4[i][t_4]+2.*k1bis);
    k3 = dt * (
        ADC_in[i] * N_in[i]->get_S((t-1 - T_in[i] + bg_.max_tau) % bg_.max_tau)
        - D2_in[i] * (Hp_in_RK4[i][t_4]+3.*k2)
        - DD_in[i] * (H_in_RK4[i][t_4]+3.*k2bis)
        );
    k3bis = dt * (Hp_in_RK4[i][t_4]+3.*k2bis);
    Hp_in_RK4[i][t] = Hp_in_RK4[i][t_4] + (8./9.) * k1 + (4./3.) * k2 + (16./9.) * k3;
    H_in_RK4[i][t] = H_in_RK4[i][t_4] + (8./9.) * k1bis + (4./3.) * k2bis + (16./9.) * k3bis;
#endif

    //// Newmark (tested only on a subproblem)
    //float gamma=0.5;
    //float beta=0.25;
    //Hpp_in_RK4[i][t] = 
    //  (
    //   dtADC_in[i] * N_in[i]->get_S((t - T_in[i] - 1 + bg_.max_tau ) % bg_.max_tau)
    //   - dtD2_in[i] * Hp_in_RK4[i][t_1]
    //   - dtDD_in[i] * H_in_RK4[i][t_1]
    //  );
    //H_in_RK4[i][t] = H_in_RK4[i][t_1] + dt * Hp_in_RK4[i][t_1] + ((dt * dt)/2) * ((1-2*beta)*Hpp_in_RK4[i][t_1] + 2*beta*Hpp_in_RK4[i][t]);
    //Hp_in_RK4[i][t] = Hp_in_RK4[i][t_1] + dt * ((1-gamma)*Hpp_in_RK4[i][t_1] + gamma*Hpp_in_RK4[i][t]);

    // Below are other integration schemes of high-order:

    //k1 = dt * (
    //    ADC_in[i] * N_in[i]->get_S((t-3 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_3])
    //    - DD_in[i] * H_in_RK4[i][t_3]
    //    );
    //k2 = dt * (
    //    ADC_in[i] * N_in[i]->get_S((t-2 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_3]+k1)
    //    - DD_in[i] * H_in_RK4[i][t_2]
    //    );
    //k3 = dt * (
    //    ADC_in[i] * N_in[i]->get_S((t-2 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_3]+k2)
    //    - DD_in[i] * H_in_RK4[i][t_2]
    //    );
    //k4 = dt * (
    //    ADC_in[i] * N_in[i]->get_S((t-1 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_3]+2.*k3)
    //    - DD_in[i] * H_in_RK4[i][t_1]
    //    );
    //Hp_in_RK4[i][t] = Hp_in_RK4[i][t_3] + (1./3.) * (k1 + 2. * k2 + 2. * k3 + k4);
    //k1 = dt * (Hp_in_RK4[i][t_3]);
    //k2 = dt * (Hp_in_RK4[i][t_2]);
    //k3 = dt * (Hp_in_RK4[i][t_2]);
    //k4 = dt * (Hp_in_RK4[i][t_1]);
    //H_in_RK4[i][t] = H_in_RK4[i][t_1] + (1./6.) * (k1 + 2. * k2 + 2. * k3 + k4);

    //k1 = dt * (
    //    ADC_in[i] * N_in[i]->get_S((t-4 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_4])
    //    - DD_in[i] * H_in_RK4[i][t_4]
    //    );
    //k2 = dt * (
    //    ADC_in[i] * N_in[i]->get_S((t-2 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_2])
    //    - DD_in[i] * H_in_RK4[i][t_2]
    //    );
    //k3 = dt * (
    //    ADC_in[i] * N_in[i]->get_S((t-1 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_1])
    //    - DD_in[i] * H_in_RK4[i][t_1]
    //    );
    //Hp_in_RK4[i][t] = Hp_in_RK4[i][t_4] + ((8./9.) * k1 + (4./3.) * k2 + (16./9.) * k3);

    //k1 = dt90 * (
    //    ADC_in[i] * N_in[i]->get_S((t-90 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_90])
    //    - DD_in[i] * H_in_RK4[i][t_90]
    //    );
    //k2 = dt90 * (
    //    ADC_in[i] * N_in[i]->get_S((t-72 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_90] + 1./5. * k1)
    //    - DD_in[i] * H_in_RK4[i][t_72]
    //    );
    //k3 = dt90 * (
    //    ADC_in[i] * N_in[i]->get_S((t-63 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_90] + 3./40. * k1 + 9./40. * k2)
    //    - DD_in[i] * H_in_RK4[i][t_63]
    //    );
    //k4 = dt90 * (
    //    ADC_in[i] * N_in[i]->get_S((t-18 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_90] + 44./45. * k1 - 56./15. * k2 + 32./9. * k3)
    //    - DD_in[i] * H_in_RK4[i][t_18]
    //    );
    //k5 = dt90 * (
    //    ADC_in[i] * N_in[i]->get_S((t-10 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_90] + 19372./6561. * k1 - 25360./2187. * k2 + 64448./6561. * k3 - 212./729. * k4)
    //    - DD_in[i] * H_in_RK4[i][t_10]
    //    );
    //k6 = dt90 * (
    //    ADC_in[i] * N_in[i]->get_S((t - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - D2_in[i] * (Hp_in_RK4[i][t_90] + 9017./3168. * k1 - 355./33. * k2 + 46732./5247. * k3 + 49./176. * k4 - 5103./18656. * k5)
    //    - DD_in[i] * (H_in_RK4[i][t])
    //    );
    //Hp_in_RK4[i][t] = 
    //  Hp_in_RK4[i][t_90]
    //  + 35./384. * k1 + 500./1113. * k3 + 125./192. * k4 - 2187./6784. * k5 + 11./84. * k6;

		sum_v += H_in_RK4[i][t] * Sign_in[i];
	}
	S[t] = Smax / (1 + exp(kvh - k*sum_v));
}

void BCBG2::SingleChannelNucleus::update_single_channel_nucleus_damping_vana()
{
  int t = bg_.tmod;
  int t_1  = (t - 1 + bg_.max_tau) % bg_.max_tau;
  int t_2  = (t - 2 + bg_.max_tau) % bg_.max_tau;
  int t_4  = (t - 4 + bg_.max_tau) % bg_.max_tau;
  if (is_damped) {
    //// Bogacki-Shampine (ode23 in matlab) for damping
    //k1 = dt * (Phip_in_RK4[t_4]);
    //k2 = dt * (Phip_in_RK4[t_2] + 2. * dt * k1);
    //k3 = dt * (Phip_in_RK4[t_1] + 3. * dt * k2);
    //Phi_in_RK4[t] = Phi_in_RK4[t_4] + ((8./9.) * k1 + (4./3.) * k2 + (16./9.) * k3);
    
    // explicit Euler for damping Phi
    Phi_in_RK4[t] = Phi_in_RK4[t_1] + dt * Phip_in_RK4[t_1];
    // explicit Euler for damping Phip
    Phip_in_RK4[t] = Phip_in_RK4[t_1] + dt *
      (
       15625 * this->get_S((t-1 +bg_.max_tau) % bg_.max_tau)
       - 250 * Phip_in_RK4[t_1]
       - 15625 * Phi_in_RK4[t_1]
      );
  }
}

void BCBG2::SingleChannelNucleus::update_single_channel_nucleus_vana()
{
  // should be optimised by setting the dt factor inside parameters
	int i;
	float sum_v = 0;
  float gamma=0.5;
  float beta=0.25;

  int t = bg_.tmod;
  int t_1  = (t - 1 + bg_.max_tau) % bg_.max_tau;
  int t_2  = (t - 2 + bg_.max_tau) % bg_.max_tau;
  int t_3  = (t - 3 + bg_.max_tau) % bg_.max_tau;
  int t_4  = (t - 4 + bg_.max_tau) % bg_.max_tau;
  int t_5  = (t - 5 + bg_.max_tau) % bg_.max_tau;
  int t_6  = (t - 6 + bg_.max_tau) % bg_.max_tau;
  int t_7  = (t - 7 + bg_.max_tau) % bg_.max_tau;

  int c=-1;
  int t_90c = (t - 90 + c + bg_.max_tau) % bg_.max_tau;
  int t_72c = (t - 72 + c + bg_.max_tau) % bg_.max_tau;
  int t_63c = (t - 63 + c + bg_.max_tau) % bg_.max_tau;
  int t_18c  = (t - 18 + c + bg_.max_tau) % bg_.max_tau;
  int t_10c = (t - 10 + c + bg_.max_tau) % bg_.max_tau;
  int tc = (t - 1 + c + bg_.max_tau) % bg_.max_tau;

  int t_90 = (t - 90 + bg_.max_tau) % bg_.max_tau;
  int t_72 = (t - 72 + bg_.max_tau) % bg_.max_tau;
  int t_63 = (t - 63 + bg_.max_tau) % bg_.max_tau;
  int t_18  = (t - 18 + bg_.max_tau) % bg_.max_tau;
  int t_10 = (t - 10 + bg_.max_tau) % bg_.max_tau;

  float dt90 = 90*dt;

	for (i=0; i<n; i++) {
    //// explicit Euler
		//H_in_RK4[i][t] = H_in_RK4[i][t_1] + dt * Hp_in_RK4[i][t_1];
    //// explicit trapezoidal (Heun)
    //H_in_RK4[i][t] = H_in_RK4[i][t_1] + 0.5 * dt *
    //  (
    //   Hp_in_RK4[i][t_1]
    //   +
    //   (H_in_RK4[i][t_1] + dt * Hp_in_RK4[i][t_1])
    //  );
    // Bogacki-Shampine (ode23)
    k1 = dt * (Hp_in_RK4[i][t_4]);
    k2 = dt * (Hp_in_RK4[i][t_2] + 2. * dt * k1);
    k3 = dt * (Hp_in_RK4[i][t_1] + 3. * dt * k2);
    H_in_RK4[i][t] = H_in_RK4[i][t_4] + ((8./9.) * k1 + (4./3.) * k2 + (16./9.) * k3);
    //// 2-steps Adams-Bashforth
		//H_in_RK4[i][t] = H_in_RK4[i][t_1] + dt * ((3./2.)*Hp_in_RK4[i][t_1] - (1./2.)*Hp_in_RK4[i][t_2]);
    //// 3-steps Adams-Bashforth ALT
		//H_in_RK4[i][t] = H_in_RK4[i][t_1] + dt * ((23./12.)*Hp_in_RK4[i][t_1] - (4./3.)*Hp_in_RK4[i][t_2] + (5./12.)*Hp_in_RK4[i][t_3]);
    //// 4-steps Adams-Bashforth
		//H_in_RK4[i][t] = H_in_RK4[i][t_1] + dt * ((55./24.)*Hp_in_RK4[i][t_1] - (59./24.)*Hp_in_RK4[i][t_2] + (37./24.)*Hp_in_RK4[i][t_3] - (3./8.)*Hp_in_RK4[i][t_4]);
    //// 5-steps Adams-Bashforth
		//H_in_RK4[i][t] = H_in_RK4[i][t_1] + dt * ((1901./720.)*Hp_in_RK4[i][t_1] - (1387./360.)*Hp_in_RK4[i][t_2] + (109./30.)*Hp_in_RK4[i][t_3] - (637./360.)*Hp_in_RK4[i][t_4] + (251./720.)*Hp_in_RK4[i][t_5]);
    //// Dormand Prince (approximated)
    //k1 = dt90 * (Hp_in_RK4[i][t_90]);
    //k2 = dt90 * (Hp_in_RK4[i][t_72] + 1./5. * k1);
    //k3 = dt90 * (Hp_in_RK4[i][t_63] + 3./40. * k1 + 9./40. * k2);
    //k4 = dt90 * (Hp_in_RK4[i][t_18] + 44./45. * k1 - 56./15. * k2 + 32./9. * k3);
    //k5 = dt90 * (Hp_in_RK4[i][t_10] + 19372./6561. * k1 - 25360./2187. * k2 + 64448./6561. * k3 - 212./729. * k4);
    //k6 = dt90 * (Hp_in_RK4[i][t_1] + 9017./3168. * k1 - 355./33. * k2 + 46732./5247. * k3 + 49./176. * k4 - 5103./18656. * k5);
    //H_in_RK4[i][t] = H_in_RK4[i][t_90] + 35./384. * k1 + 500./1113. * k3 + 125./192. * k4 - 2187./6784. * k5 + 11./84. * k6;
    //// Dormand Prince
    //k1 = dt90 * (Hp_in_RK4[i][t_90c]);
    //k2 = dt90 * (Hp_in_RK4[i][t_72c] + 1./5. * k1);
    //k3 = dt90 * (Hp_in_RK4[i][t_63c] + 3./40. * k1 + 9./40. * k2);
    //k4 = dt90 * (Hp_in_RK4[i][t_18c] + 44./45. * k1 - 56./15. * k2 + 32./9. * k3);
    //k5 = dt90 * (Hp_in_RK4[i][t_10c] + 19372./6561. * k1 - 25360./2187. * k2 + 64448./6561. * k3 - 212./729. * k4);
    //k6 = dt90 * (Hp_in_RK4[i][tc] + 9017./3168. * k1 - 355./33. * k2 + 46732./5247. * k3 + 49./176. * k4 - 5103./18656. * k5);
    //H_in_RK4[i][t] = H_in_RK4[i][t_90c] + 35./384. * k1 + 500./1113. * k3 + 125./192. * k4 - 2187./6784. * k5 + 11./84. * k6;

    //// explicit Euler
    //Hp_in_RK4[i][bg_.tmod] = Hp_in_RK4[i][t_1] + dt *
    //  (
    //   102400 * nu_in[i] * N_in[i]->get_S((t-1 - T_in[i]) % bg_.max_tau)
    //   - 800 * Hp_in_RK4[i][t_1]
    //   - 102400 * H_in_RK4[i][t_1]
    //  );
    //// explicit Euler with better estimation of H (needs prior computation of H)
    //Hp_in_RK4[i][bg_.tmod] = Hp_in_RK4[i][t_1] + dt *
    //  (
    //   102400 * nu_in[i] * N_in[i]->get_S((t-1 - T_in[i]) % bg_.max_tau)
    //   - 800 * Hp_in_RK4[i][t_1]
    //   - 102400 * H_in_RK4[i][t]
    //  );
    //// "improved" explicit Euler (mix current input with current H and former Hp)
    //Hp_in_RK4[i][bg_.tmod] = Hp_in_RK4[i][t_1] + dt *
    //  (
    //   102400 * nu_in[i] * N_in[i]->get_S((t - T_in[i]) % bg_.max_tau)
    //   - 800 * Hp_in_RK4[i][t_1]
    //   - 102400 * H_in_RK4[i][t]
    //  );
    //// Heun (explicit trapezoidal)
    //Hp_in_RK4[i][t] = Hp_in_RK4[i][t_1] + 0.5 * dt *
    //  (
    //   (102400 * nu_in[i] * N_in[i]->get_S((t-1 - T_in[i]) % bg_.max_tau) - 800 *
    //     Hp_in_RK4[i][t_1]
    //     - 102400 *
    //     H_in_RK4[i][t_1])
    //   +
    //   (102400 * nu_in[i] * N_in[i]->get_S((t - T_in[i]) % bg_.max_tau) - 800 *
    //    (
    //     Hp_in_RK4[i][t_1] + dt *
    //     (102400 * nu_in[i] * N_in[i]->get_S((t-1 - T_in[i]) % bg_.max_tau) - 800 * Hp_in_RK4[i][t_1] - 102400 * H_in_RK4[i][t_1])
    //     )
    //    - 102400 *
    //    H_in_RK4[i][t_1])
    //  );
    ////bsf Heun (explicit trapezoidal) with better estimation of H (needs prior computation of H)
    //Hp_in_RK4[i][t] = Hp_in_RK4[i][t_1] + 0.5 * dt *
    //  (
    //   (102400 * nu_in[i] * N_in[i]->get_S((t-1 - T_in[i]) % bg_.max_tau) - 800 *
    //     Hp_in_RK4[i][t_1]
    //     - 102400 *
    //     H_in_RK4[i][t_1])
    //   +
    //   (102400 * nu_in[i] * N_in[i]->get_S((t - T_in[i]) % bg_.max_tau) - 800 *
    //    (
    //     Hp_in_RK4[i][t_1] + dt *
    //     (102400 * nu_in[i] * N_in[i]->get_S((t-1 - T_in[i]) % bg_.max_tau) - 800 * Hp_in_RK4[i][t_1] - 102400 * H_in_RK4[i][t_1])
    //     )
    //    - 102400 *
    //    H_in_RK4[i][t])
    //  );
    //// RK4 with better estimation of Hp (needs prior computation of Hp)
    //k1 = 102400 * nu_in[i] * N_in[i]->get_S((t-2 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //  - 800 * (Hp_in_RK4[i][t_2])
    //  - 102400 * H_in_RK4[i][t_2];
    //k2 = 102400 * nu_in[i] * N_in[i]->get_S((t-1 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //  - 800 * (Hp_in_RK4[i][t_2] + (dt * k1))
    //  - 102400 * H_in_RK4[i][t_1];
    //k3 = 102400 * nu_in[i] * N_in[i]->get_S((t-1 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //  - 800 * (Hp_in_RK4[i][t_2] + (dt * k2))
    //  - 102400 * H_in_RK4[i][t_1];
    //k4 = 102400 * nu_in[i] * N_in[i]->get_S((t - T_in[i]) % bg_.max_tau)
    //  - 800 * (Hp_in_RK4[i][t_2] + 2 * dt * k3)
    //  - 102400 * H_in_RK4[i][t];
    //Hp_in_RK4[i][bg_.tmod] = Hp_in_RK4[i][t_2] + dt * (k1 + 2*k2 + 2*k3 + k4)/3;
    //// 5-steps Adams-Bashforth
    //Hp_in_RK4[i][t] = Hp_in_RK4[i][t_1] + dt * (
    //    (1901./720.)*
    //    (102400 * nu_in[i] * N_in[i]->get_S((t-1 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //     - 800 * (Hp_in_RK4[i][t_1])
    //     - 102400 * H_in_RK4[i][t_1])
    //    - (1387./360.)*
    //    (102400 * nu_in[i] * N_in[i]->get_S((t-2 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //     - 800 * (Hp_in_RK4[i][t_2])
    //     - 102400 * H_in_RK4[i][t_2])
    //    + (109./30.)*
    //    (102400 * nu_in[i] * N_in[i]->get_S((t-3 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //     - 800 * (Hp_in_RK4[i][t_3])
    //     - 102400 * H_in_RK4[i][t_3])
    //    - (637./360.)*
    //    (102400 * nu_in[i] * N_in[i]->get_S((t-4 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //     - 800 * (Hp_in_RK4[i][t_4])
    //     - 102400 * H_in_RK4[i][t_4])
    //    + (251./720.)*
    //    (102400 * nu_in[i] * N_in[i]->get_S((t-5 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //     - 800 * (Hp_in_RK4[i][t_5])
    //     - 102400 * H_in_RK4[i][t_5])
    //    );
    //// Dormand Prince (approximated)
    //k1 = dt90 * (
    //    102400 * nu_in[i] * N_in[i]->get_S((t-90 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - 800 * (Hp_in_RK4[i][t_90])
    //    - 102400 * H_in_RK4[i][t_90]
    //    );
    //k2 = dt90 * (
    //    102400 * nu_in[i] * N_in[i]->get_S((t-72 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - 800 * (Hp_in_RK4[i][t_90] + 1./5. * k1)
    //    - 102400 * H_in_RK4[i][t_72]
    //    );
    //k3 = dt90 * (
    //    102400 * nu_in[i] * N_in[i]->get_S((t-63 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - 800 * (Hp_in_RK4[i][t_90] + 3./40. * k1 + 9./40. * k2)
    //    - 102400 * H_in_RK4[i][t_63]
    //    );
    //k4 = dt90 * (
    //    102400 * nu_in[i] * N_in[i]->get_S((t-18 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - 800 * (Hp_in_RK4[i][t_90] + 44./45. * k1 - 56./15. * k2 + 32./9. * k3)
    //    - 102400 * H_in_RK4[i][t_18]
    //    );
    //k5 = dt90 * (
    //    102400 * nu_in[i] * N_in[i]->get_S((t-10 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - 800 * (Hp_in_RK4[i][t_90] + 19372./6561. * k1 - 25360./2187. * k2 + 64448./6561. * k3 - 212./729. * k4)
    //    - 102400 * H_in_RK4[i][t_10]
    //    );
    //k6 = dt90 * (
    //    102400 * nu_in[i] * N_in[i]->get_S((t-1 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - 800 * (Hp_in_RK4[i][t_90] + 9017./3168. * k1 - 355./33. * k2 + 46732./5247. * k3 + 49./176. * k4 - 5103./18656. * k5)
    //    - 102400 * (H_in_RK4[i][t-1])
    //    );
    //Hp_in_RK4[i][t] = Hp_in_RK4[i][t_90] + 35./384. * k1 + 500./1113. * k3 + 125./192. * k4 - 2187./6784. * k5 + 11./84. * k6;
    // Dormand Prince
    k1 = dt90 * (
        102400 * nu_in[i] * N_in[i]->get_S((t-90 - T_in[i] + bg_.max_tau) % bg_.max_tau)
        - 800 * (Hp_in_RK4[i][t_90])
        - 102400 * H_in_RK4[i][t_90]
        );
    k2 = dt90 * (
        102400 * nu_in[i] * N_in[i]->get_S((t-72 - T_in[i] + bg_.max_tau) % bg_.max_tau)
        - 800 * (Hp_in_RK4[i][t_90] + 1./5. * k1)
        - 102400 * H_in_RK4[i][t_72]
        );
    k3 = dt90 * (
        102400 * nu_in[i] * N_in[i]->get_S((t-63 - T_in[i] + bg_.max_tau) % bg_.max_tau)
        - 800 * (Hp_in_RK4[i][t_90] + 3./40. * k1 + 9./40. * k2)
        - 102400 * H_in_RK4[i][t_63]
        );
    k4 = dt90 * (
        102400 * nu_in[i] * N_in[i]->get_S((t-18 - T_in[i] + bg_.max_tau) % bg_.max_tau)
        - 800 * (Hp_in_RK4[i][t_90] + 44./45. * k1 - 56./15. * k2 + 32./9. * k3)
        - 102400 * H_in_RK4[i][t_18]
        );
    k5 = dt90 * (
        102400 * nu_in[i] * N_in[i]->get_S((t-10 - T_in[i] + bg_.max_tau) % bg_.max_tau)
        - 800 * (Hp_in_RK4[i][t_90] + 19372./6561. * k1 - 25360./2187. * k2 + 64448./6561. * k3 - 212./729. * k4)
        - 102400 * H_in_RK4[i][t_10]
        );
    k6 = dt90 * (
        102400 * nu_in[i] * N_in[i]->get_S((t - T_in[i] + bg_.max_tau) % bg_.max_tau)
        - 800 * (Hp_in_RK4[i][t_90] + 9017./3168. * k1 - 355./33. * k2 + 46732./5247. * k3 + 49./176. * k4 - 5103./18656. * k5)
        - 102400 * (H_in_RK4[i][t])
        );
    Hp_in_RK4[i][t] = 
      Hp_in_RK4[i][t_90]
      + 35./384. * k1 + 500./1113. * k3 + 125./192. * k4 - 2187./6784. * k5 + 11./84. * k6;



    //// Runge-Kutta integration scheme (as defined in http://mymathlib.webtrellis.net/c_source/diffeq/second_order/runge_kutta_2nd_order.c)
    ////k1 = h * f(x[i], y[i], yp[i])
    //k1 = dt * (
    //    102400 * nu_in[i] * N_in[i]->get_S((t-2 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - 800 * (Hp_in_RK4[i][t_2])
    //    - 102400 * H_in_RK4[i][t_2]
    //    );
    ////k2 = h * f(x[i]+h/2, y[i]+(h/2)*yp[i], yp[i]+k1/2)
    //k2 = dt * (
    //    102400 * nu_in[i] * N_in[i]->get_S((t-1 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - 800 * (Hp_in_RK4[i][t_2] + k1/2.)
    //    - 102400 * (H_in_RK4[i][t_2] + dt * Hp_in_RK4[i][t_2])
    //    );
    ////k3 = h * f(x[i]+h/2, y[i]+(h/2)*yp[i]+(h/4)*k1, yp[i]+k2/2)
    //k3 = dt * (
    //    102400 * nu_in[i] * N_in[i]->get_S((t-1 - T_in[i] + bg_.max_tau) % bg_.max_tau)
    //    - 800 * (Hp_in_RK4[i][t_2] + k2/2.)
    //    - 102400 * (H_in_RK4[i][t_2] + dt * Hp_in_RK4[i][t_2] + (dt/2.)*k1)
    //    );
    ////y[i+1] = y[i] + h * ( yp[i] + 1/6 (k1 + k2 + k3 )
    //H_in_RK4[i][t] = H_in_RK4[i][t_2] + dt * (2. * Hp_in_RK4[i][t_2] + k1 + k2 + k3)/3.;
    ////yp[i+1] = yp[i] + 1/6 * ( k1 + 2 k2 + 2 k3 + k4 )
    //Hp_in_RK4[i][t] = Hp_in_RK4[i][t_2] + (k1 + 2. * k2 + 2. * k3 + k4)/6.;

    //// Newmark integration scheme
    //Hpp_in_RK4[i][t] = (102400 * nu_in[i] * N_in[i]->get_S(t) - 800 * Hp_in_RK4[i][t_1] - 102400 * H_in_RK4[i][t_1]);
    //H_in_RK4[i][t] = H_in_RK4[i][t_1] + dt * Hp_in_RK4[i][t_1] + ((dt * dt)/2) * ((1-2*beta)*Hpp_in_RK4[i][t_1] + 2*beta*Hpp_in_RK4[i][t]);
    //Hp_in_RK4[i][t] = Hp_in_RK4[i][t_1] + dt * ((1-gamma)*Hpp_in_RK4[i][t_1] + gamma*Hpp_in_RK4[i][t]);

		sum_v += H_in_RK4[i][t];
	}
	S[(bg_.tmod) % bg_.max_tau] = Smax / (1 + exp(kvh - k*sum_v));
}

void BCBG2::SingleChannelNucleus::display(std::ostream& os)
{
  os << "self-description of " << get_name() << std::endl;
  os << "  Smax = " << Smax << std::endl;
  os << "  vh = " << vh << std::endl;
  os << "  k = " << k << std::endl;
  for (size_t i=0;i<n;i++) {
    os << "  afferent n° " << i << " is from " << N_in[i]->get_name() << std::endl;
    os << "    A = " << A_in[i] << std::endl;
    os << "    D = " << D_in[i] << std::endl;
    os << "    T = " << T_in[i] << std::endl;
    os << "    C = " << C_in[i] << std::endl;
    os << "    nu = " << nu_in[i] << std::endl;
    os << "    dtADC = " << dtADC_in[i] << std::endl;
    }
}

void BCBG2::SingleChannelNucleus::save() {
  int i,j;

  S_backup.resize(S.size());
  for (i=0;i<S.size(); i++) {
    S_backup[i] = S[i];
  }


  H_in_RK4_backup.resize(H_in_RK4.size());
  for (i=0;i<H_in_RK4.size(); i++) {
    H_in_RK4_backup[i].resize(H_in_RK4[i].size());
    for (j=0;j<H_in_RK4[i].size(); j++) {
      H_in_RK4_backup[i][j] = H_in_RK4[i][j];
    }
  }

  Hp_in_RK4_backup.resize(Hp_in_RK4.size());
  for (i=0;i<Hp_in_RK4.size(); i++) {
    Hp_in_RK4_backup[i].resize(Hp_in_RK4[i].size());
    for (j=0;j<Hp_in_RK4[i].size(); j++) {
      Hp_in_RK4_backup[i][j] = Hp_in_RK4[i][j];
    }
  }
}

void BCBG2::SingleChannelNucleus::save(MemorySCN& mem) {
  int i,j;

  mem.S_backup.resize(S.size());
  for (i=0;i<S.size(); i++) {
    mem.S_backup[i] = S[i];
  }

  mem.H_in_RK4_backup.resize(H_in_RK4.size());
  for (i=0;i<H_in_RK4.size(); i++) {
    mem.H_in_RK4_backup[i].resize(H_in_RK4[i].size());
    for (j=0;j<H_in_RK4[i].size(); j++) {
      mem.H_in_RK4_backup[i][j] = H_in_RK4[i][j];
    }
  }

  mem.Hp_in_RK4_backup.resize(Hp_in_RK4.size());
  for (i=0;i<Hp_in_RK4.size(); i++) {
    mem.Hp_in_RK4_backup[i].resize(Hp_in_RK4[i].size());
    for (j=0;j<Hp_in_RK4[i].size(); j++) {
      mem.Hp_in_RK4_backup[i][j] = Hp_in_RK4[i][j];
    }
  }
}

void BCBG2::SingleChannelNucleus::load() {
  int i,j;

  S.resize(S_backup.size());
  for (i=0;i<S.size(); i++) {
    S[i] = S_backup[i];
  }

  H_in_RK4.resize(H_in_RK4_backup.size());
  for (i=0;i<H_in_RK4_backup.size(); i++) {
    H_in_RK4[i].resize(H_in_RK4_backup[i].size());
    for (j=0;j<H_in_RK4_backup[i].size(); j++) {
      H_in_RK4[i][j] = H_in_RK4_backup[i][j];
    }
  }

  Hp_in_RK4.resize(Hp_in_RK4_backup.size());
  for (i=0;i<Hp_in_RK4_backup.size(); i++) {
    Hp_in_RK4[i].resize(Hp_in_RK4_backup[i].size());
    for (j=0;j<Hp_in_RK4_backup[i].size(); j++) {
      Hp_in_RK4[i][j] = Hp_in_RK4_backup[i][j];
    }
  }
}

void BCBG2::SingleChannelNucleus::load(MemorySCN& mem) {
  int i,j;

  S.resize(mem.S_backup.size());
  for (i=0;i<S.size(); i++) {
    S[i] = mem.S_backup[i];
  }

  H_in_RK4.resize(mem.H_in_RK4_backup.size());
  for (i=0;i<mem.H_in_RK4_backup.size(); i++) {
    H_in_RK4[i].resize(mem.H_in_RK4_backup[i].size());
    for (j=0;j<mem.H_in_RK4_backup[i].size(); j++) {
      H_in_RK4[i][j] = mem.H_in_RK4_backup[i][j];
    }
  }

  Hp_in_RK4.resize(mem.Hp_in_RK4_backup.size());
  for (i=0;i<mem.Hp_in_RK4_backup.size(); i++) {
    Hp_in_RK4[i].resize(mem.Hp_in_RK4_backup[i].size());
    for (j=0;j<mem.Hp_in_RK4_backup[i].size(); j++) {
      Hp_in_RK4[i][j] = mem.Hp_in_RK4_backup[i][j];
    }
  }
}

