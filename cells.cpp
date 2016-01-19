#include "bcbg2.hpp"

// the implementation of nucleus cell is not finished, and code for it is not released. The skeleton of functions is left here as to not break the compilation...

void BCBG2::NucleusCells::initialize_cells(int cells_n, float Smax, float Sini, float V_threshold, float k, float V_rest, float R, float Rm, float Ri, float* dists, float* diams, int compartment_nb, const char* id, int nucleus_number)
{
}

int BCBG2::NucleusCells::affect_ion_channel(float A, float D, float Vrev, bool is_nmda)
{
  return(0);
}

void BCBG2::NucleusCells::add_afferent_axon(NucleusCells* nc, float cs, int T, float C, float distance_factor, int ion_channel_n)
{
}

void BCBG2::NucleusCells::add_efferent_axon(NucleusCells* nc, float cs, int T, float C, float distance_factor, int ion_channel_n)
{
}

void BCBG2::NucleusCells::add_afferent_synapse(NucleusCells* source_nucleus, float cs, int T, float C, float distance_factor, int ion_channel_n)
{
}

bool BCBG2::NucleusCells::set_afferent(float A, float D, int S, float C, int T, float cs, float distance, NucleusCells* N)
{
  return(0);
}

bool BCBG2::NucleusCells::add(float A, float D, float Vrev, float C, int T, float cs, float distance, NucleusCells* N)
{
  return(0);
}

bool BCBG2::NucleusCells::nmda(float A, float D, float Vrev, float C, int T, float cs, float distance, NucleusCells* N)
{
  return(0);
}

bool BCBG2::NucleusCells::add_ampa_gaba_nmda(float A, float D, float Vrev, float C, int T, float cs, float distance, NucleusCells* N, bool is_nmda)
{
  return(0);
}

void BCBG2::NucleusCells::set_dt(float value)
{
}


void BCBG2::NucleusCells::fast_add_g(int cell_i, int ion_channel_i, float distance_factor, float t_sent)
{
}


void BCBG2::NucleusCells::add_g(int cell_i, int ion_channel_i, float distance_factor, float t_sent)
{
}


void BCBG2::NucleusCells::update_cells()
{
}

void BCBG2::NucleusCells::set_rate(float rate)
{
}

void BCBG2::NucleusCells::set_rate(float rate, bool reset)
{
}

void BCBG2::NucleusCells::prepare_next_iteration()
{
}

void BCBG2::NucleusCells::display_summary(std::ostream& os)
{
}


void BCBG2::NucleusCells::display(std::ostream& os)
{
}

void BCBG2::NucleusCells::save()
{
}

void BCBG2::NucleusCells::restore()
{
}
