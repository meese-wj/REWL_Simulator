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
    int my_world_rank = REWL_MASTER_PROC;
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

    simulate();
};

REWL_simulation::REWL_simulation()
{
#if MPI_ON
    MPI_Comm_rank( &my_world_rank, MPI_COMM_WORLD );
    MPI_Comm_size( &num_walkers, MPI_COMM_WORLD );
#endif
    
    i_am_the_master = ( my_world_rank == REWL_MASTER_PROC );

    // Construct the glazier
    window_maker = new glazier<ENERGY_TYPE, histogram_index<ENERGY_TYPE> >
                        (System_Parameters::energy_min, System_Parameters::energy_max,
                         System_Parameters::energy_bin_size, 
                         static_cast<size_t>(REWL_Parameters::num_walkers),
                         REWL_Parameters::replicas_per_window, 
                         static_cast<ENERGY_TYPE>(REWL_Parameters::window_overlap));
    
    window_maker -> construct_windows();

    // Construct the walker
    ENERGY_TYPE walker_min = window_maker -> all_windows[my_world_rank].minimum;
    ENERGY_TYPE walker_max = window_maker -> all_windows[my_world_rank].maximum;
    ENERGY_TYPE walker_bin_size = window_maker -> all_windows[my_world_rank].bin_size;
    size_t walker_num_bins = window_maker -> all_windows[my_world_rank].num_bins;

    // TODO: Set up the timer as the seed.
    std::uint32_t walker_seed = static_cast<std::uint32_t> (1);

    my_walker = new REWL_Walker<ENERGY_TYPE, LOGDOS_TYPE, OBS_TYPE, histogram_index<ENERGY_TYPE> >
                (walker_min, walker_max, walker_bin_size, walker_seed);

}

#endif
