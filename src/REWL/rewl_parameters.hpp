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
    float window_overlap = 0.75;

    constexpr size_t sweeps_per_check = 1000;
    constexpr size_t sweeps_per_exchange = 1;
    constexpr LOGDOS_TYPE final_increment = 1e-6;
    //constexpr LOGDOS_TYPE final_increment = 0.25;
    
    // This does not matter if the 1/t algorithm is
    // used to converge the logDoS
    constexpr float flatness_criterion = 0.3;

    // How many iterations to run through before
    // outputting any data to stdout
    constexpr size_t iterations_per_stdout = 20;
}

#endif
