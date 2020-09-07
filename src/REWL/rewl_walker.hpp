#ifndef REWL_WALKER
/* Create a struct which is really
 * just a wrapper around all of the 
 * relevant data types for a single
 * MPI process in REWL. */

#include <model_hamiltonians.hpp>
#include "wang_landau.hpp"

#if MPI_ON
#include <mpi_rewl_helpers.hpp>
#endif

template<typename energy_t, typename logdos_t, typename obs_t, class histogram_index_functor>
struct REWL_Walker
{
    int walker_world_rank = 0;
    int i_am_the_master = 0;
   
    logdos_t incrementer = 1.;
    histogram_index_functor hist_idx;
    rng<float> random;
    Wang_Landau<energy_t, logdos_t, Hamiltonian_t<obs_t>, 
                Observables_t<obs_t>, State_t<obs_t>, histogram_index_functor> wl_walker;
    Hamiltonian_t<obs_t> system;
    Observables_t<obs_t> system_obs;


    REWL_Walker(const energy_t _min, const energy_t _max, const energy_t _bsize, const size_t _nbins, const std::uint32_t _seed);
    ~REWL_Walker() {}

    void wang_landau_walk(const size_t num_sweeps);

};

template<typename energy_t,
         typename logdos_t, 
         typename obs_t, 
         class histogram_index_functor>
REWL_Walker<energy_t, 
            logdos_t, 
            obs_t, 
            histogram_index_functor>::
            REWL_Walker(const energy_t _min, 
                        const energy_t _max, 
                        const energy_t _bsize, 
                        const size_t _nbins, 
                        const std::uint32_t _seed)
                      : 
                        hist_idx(_min, _max, _bsize),
                        random(_seed),
                        wl_walker(_min, _max, _bsize, _nbins),
                        system(),
                        system_obs(_nbins)
{
    bool in_range = hist_idx.energy_in_range(system.current_state.energy);
    while ( !in_range )
    {
       size_t site = static_cast<size_t>( random() * System_Parameters::N );
       State_t<obs_t> temporary_state;
       system.change_state(site, temporary_state);
       
       if ( hist_idx.energy_too_low(system.current_state.energy) && temporary_state.energy > system.current_state.energy )
           system.set_state(site, temporary_state);
       else if ( hist_idx.energy_too_high(system.current_state.energy) && temporary_state.energy <= system.current_state.energy )
           system.set_state(site, temporary_state);

       in_range = hist_idx.energy_in_range(system.current_state.energy);
    }

#if MPI_ON
    MPI_Comm_rank( &walker_world_rank, MPI_COMM_WORLD );
    i_am_the_master = ( walker_world_rank == REWL_MASTER_PROC ? 1 : 0 );
#endif 
}

template<typename energy_t, typename logdos_t, typename obs_t, class histogram_index_functor>
void REWL_Walker<energy_t, logdos_t, obs_t, histogram_index_functor>::wang_landau_walk(const size_t num_sweeps)
{
    size_t system_size = System_Parameters::N;
    for ( size_t sweep = 0; sweep != num_sweeps; ++sweep )
    {
        wl_walker.wang_landau_sweep(system_size, incrementer, &system, &system_obs, random, hist_idx);    
    }
}

#endif
