#ifndef REWL_SIMULATION
/* This file contains the templated simulation
 * struct that will initialize the simulation
 * and contain all the relevant types. */

#include <histogram_index.hpp>
#include <glazier.hpp>
#include <rewl_walker.hpp>

// Create an alias for the observable
// data type. 
using ENERGY_TYPE = float;
using LOGDOS_TYPE = double;
using OBS_TYPE = float;

namespace REWL_Parameters
{
    int num_walkers = 1;
    size_t replicas_per_window = 1;
    float window_overlap = static_cast<float>( single_bin_overlap ); 

    constexpr size_t sweeps_per_check = 1000;
    constexpr LOGDOS_TYPE final_increment = 1e-6;
    constexpr float flatness_criterion = 0.3;
}

struct REWL_simulation
{
#if MPI_ON
    int my_world_rank = REWL_MASTER_PROC;
#else
    int my_world_rank = 0;
#endif
    bool i_am_the_master = true;

    // TODO: Make a glazier and start this stuff up.
    glazier<ENERGY_TYPE, histogram_index<ENERGY_TYPE> > * window_maker = nullptr;
    
    REWL_Walker<ENERGY_TYPE, LOGDOS_TYPE, OBS_TYPE, histogram_index<ENERGY_TYPE> > * my_walker = nullptr;

    REWL_simulation();

    ~REWL_simulation()
    { 
        if (window_maker != nullptr) delete window_maker;
        if (my_walker != nullptr) delete my_walker;
    }

    void simulate() const;
};

#endif
