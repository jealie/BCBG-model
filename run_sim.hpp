#ifndef RUN_SIM_HPP
#define RUN_SIM_HPP

#include <list>

int _run_sim( // encapsulation for one channel nucleus
    float max_duration,
    float bunch_duration,
    float dt,
    std::vector<float> &activ,
    std::vector<float> &params,
    std::vector<int> &delay,
    std::vector<float> &means,
    int verbose);

int _run_sim(
    float max_duration,
    float bunch_duration,
    float dt,
    std::vector<float> &activ,
    std::vector<float> &cs,
    std::vector<float> &params,
    std::vector<int> &delay,
    std::vector<float> &means,
    int verbose,
    int ch_n,
    int msn_separation,
    int do_checks,
    MemoryBCBG2& mem,
    int integration_method);

int _run_sim_tsirogiannis_2010(
    float max_duration,
    float bunch_duration,
    float dt,
    std::vector<float> &activations,
    std::vector<float> &params,
    std::vector<int> &delay,
    std::vector<float> &means,
    int verbose);

#endif
