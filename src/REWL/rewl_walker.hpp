#ifndef REWL_WALKER
/* Create a struct which is really
 * just a wrapper around all of the 
 * relevant data types for a single
 * MPI process in REWL. */

#include <glazier.hpp>
#include <model_hamiltonians.hpp>
#include "wang_landau.hpp"

#if MPI_ON
#include <mpi_rewl_helpers.hpp>
#endif

template<typename logdos_t, typename obs_t>
struct REWL_Walker
{
    int walker_world_rank = 0;
    int i_am_the_master = 0;
   
    logdos_t incrementer = 1.;
    // TODO: change the histogram indexer type to float.
    histogram_index<logdos_t> hist_idx;
    rng<float> random;
    Wang_Landau<logdos_t, Hamiltonian_t<obs_t>, Observables_t<obs_t>, State_t<obs_t> > wl_walker;
    Hamiltonian_t<obs_t> system;
    Observables_t<obs_t> system_obs;


    REWL_Walker(const logdos_t _min, const logdos_t _max, const logdos_t _bsize, const std::uint32_t _seed);
    ~REWL_Walker() {}

    void wang_landau_walk(const size_t num_sweeps) const;

};

template<typename logdos_t, typename obs_t>
REWL_Walker<logdos_t, obs_t>::REWL_Walker(const logdos_t _min, const logdos_t _max, const logdos_t _bsize, const std::uint32_t _seed)
                                      : hist_idx(_min, _max, _bsize),
                                        random(_seed),   
                                        wl_walker(_min, _max, _bsize),
                                        system_obs(wl_walker.wl_histograms.num_bins)
{
#if MPI_ON
    MPI_Comm_rank( &walker_world_rank, MPI_COMM_WORLD );
    i_am_the_master = ( walker_world_rank == REWL_MASTER_PROC ? 1 : 0 );
#endif 
}

template<typename logdos_t, typename obs_t>
void REWL_Walker<logdos_t, obs_t>::wang_landau_walk(const size_t num_sweeps) const
{
    size_t system_size = System_Parameters::N;
    for ( size_t sweep = 0; sweep != num_sweeps; ++sweep )
    {
        wl_walker.wang_landau_sweep(system_size, incrementer, &system, &system_obs, random, hist_idx);    
    }
}

#endif
