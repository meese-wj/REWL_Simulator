#ifndef REWL_PARAMETERS
#define REWL_PARAMETERS
/* Header file for the REWL parameters. */

#include <cmath>

// Create an alias for the observable
// data type. 
using ENERGY_TYPE = float;
using LOGDOS_TYPE = double;
using OBS_TYPE = double;

#if MPI_ON
// TODO: Are these macros really necessary for type aliasing?
#define MPI_ENERGY_TYPE MPI_FLOAT
#define MPI_LOGDOS_TYPE MPI_DOUBLE
#define MPI_OBS_TYPE MPI_DOUBLE
#endif

namespace REWL_Parameters
{
    int num_walkers = 1;
    size_t replicas_per_window = 1;
    //float window_overlap = static_cast<float>( single_bin_overlap ); 
    float window_overlap = 0.5;

    constexpr size_t sweeps_per_check = 1000;
    constexpr size_t sweeps_per_exchange = 100;
    constexpr LOGDOS_TYPE final_increment = 1e-7;
    constexpr float flatness_criterion = 0.3;
}

#endif
