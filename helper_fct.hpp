#ifndef HELPER_FCT_HPP
#define HELPER_FCT_HPP

#include <vector>

float param2hz(float value);
float param2boutons(float value);
float param2boutons2(float value);

float _do_trial(std::vector<float>& means, std::vector<float>& ref, std::vector<float>& params, std::vector<int>& delays, std::vector<float>& activations, int nucleus, float proportional_change, float proportional_radius, float sim_time, MemoryBCBG2& mem);
float calc_score_desactivation_other(std::vector <float>& means, std::vector <float>& params, std::vector <int>& delays, float desactivation_level, float sim_time, MemoryBCBG2& mem, bool verbose);
float calc_score_desactivation(std::vector <float>& means, std::vector <float>& params, std::vector <int>& delays, float desactivation_level, float sim_time, MemoryBCBG2& mem, bool verbose);
float calc_score_selective_boutons(std::vector <float>& params, bool verbose, int studied_nucleus);
float calc_score_selective_axons(std::vector <float>& params, bool verbose, int studied_nucleus);

int is_in(std::vector <int>& c, int i);
float _has_changed_near_tronqgaussian(float reference, float mean, float proportional_change, float proportional_radius);
float _has_changed_near_gaussian(float reference, float mean, float proportional_change, float proportional_radius);
float _is_near_gaussian(float reference, float mean, float absolute_radius);
float _is_near_tronqgaussian(float reference, float mean, float absolute_radius);
float _is_near_step(float reference, float mean, float absolute_radius);
#endif
