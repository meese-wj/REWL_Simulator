#ifndef REWL_PARAMETERS
/* Header file for the REWL parameters. */

// Create an alias for the observable
// data type. 
using ENERGY_TYPE = float;
using LOGDOS_TYPE = double;
using OBS_TYPE = double;

namespace REWL_Parameters
{
    int num_walkers = 1;
    size_t replicas_per_window = 1;
    float window_overlap = static_cast<float>( single_bin_overlap ); 

    constexpr size_t sweeps_per_check = 1000;
    constexpr LOGDOS_TYPE final_increment = 5e-7;
    constexpr float flatness_criterion = 0.3;
}

#endif
