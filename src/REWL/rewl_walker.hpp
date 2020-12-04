#ifndef REWL_WALKER
#define REWL_WALKER
/* Create a struct which is really
 * just a wrapper around all of the 
 * relevant data types for a single
 * MPI process in REWL. */

#include <model_hamiltonians.hpp>
#include "wang_landau.hpp"

#if MPI_ON
#include <mpi_rewl_definitions.hpp>
#endif

template<typename energy_t, typename logdos_t, typename obs_t, class histogram_index_functor>
struct REWL_Walker
{
    int walker_world_rank = 0;
    int i_am_the_master = 0;
   
    logdos_t incrementer = 1.;
    histogram_index_functor hist_idx;
    rng<float> random;
    Hamiltonian_t<obs_t> system;
    Observables_t<obs_t> system_obs;
    Wang_Landau<energy_t, logdos_t, Hamiltonian_t<obs_t>, 
                Observables_t<obs_t>, State_t<obs_t>, histogram_index_functor> wl_walker;

    REWL_Walker(const energy_t _min, const energy_t _max, const energy_t _bsize, const size_t _nbins, const std::uint32_t _seed);
    ~REWL_Walker(){}

    void adjust_state_to_range();
    energy_t get_min_energy() const;

#if SAMPLE_AFTER
    void wang_landau_walk(const size_t num_sweeps, const bool sample_observables);
#else
    void wang_landau_walk(const size_t num_sweeps);
#endif

    /* ************************************** */
    /* Include necessary REWL functions here  */

    bool energy_in_range( const energy_t energy ) const;
    float get_rand(); 
    logdos_t get_logdos( const energy_t energy ) const;
    energy_t current_energy() const;
    State_t<obs_t> * current_state();
    obs_t * DoFs() const;

    void update_histograms();
#if SAMPLE_AFTER
    void update_observables( const bool sample_observables );
#else
    void update_observables();
#endif

    /* ************************************** */

    void export_energy_bins( energy_t *& data_arr )
    {
        size_t nbins = wl_walker.wl_histograms.num_bins;
        data_arr = new energy_t [ nbins ];
        for ( size_t bin = 0; bin != nbins; ++bin )
            data_arr[ bin ] = hist_idx.get_bin_min( bin );
    }

};


// Move the system into its given 
// energy range
template<typename energy_t,
         typename logdos_t, 
         typename obs_t, 
         class histogram_index_functor>
REWL_Walker<energy_t, 
            logdos_t, 
            obs_t, 
            histogram_index_functor>::
            adjust_state_to_range()
{
#if MPI_ON
    MPI_Comm_rank( MPI_COMM_WORLD, &walker_world_rank );
    i_am_the_master = ( walker_world_rank == REWL_MASTER_PROC ? 1 : 0 );
#endif 

    bool in_range = hist_idx.energy_in_range(system.current_state.energy);
    size_t site = static_cast<size_t>( random() * System_Parameters::N );
    while ( !in_range )
    {
        ++site;

        if ( site == System_Parameters::N )
        {
            site = 0;
            system.switch_flavor_to_update();
        }

        State_t<obs_t> temporary_state;
        system.change_state(site, temporary_state);

	//printf("\nID %d: current energy = %e, temporary energy = %e", walker_world_rank, system.current_state.energy, temporary_state.energy);

        if ( hist_idx.energy_too_low(system.current_state.energy) && temporary_state.energy > system.current_state.energy )
        {
            system.set_state(site, temporary_state);
        }
        else if ( hist_idx.energy_too_high(system.current_state.energy) && temporary_state.energy <= system.current_state.energy )
        {
            system.set_state(site, temporary_state);
        }

        in_range = hist_idx.energy_in_range(system.current_state.energy);
        //printf("\nID %d: current energy = %e, in_range = %s\n", walker_world_rank, system.current_state.energy, ( in_range ? "true" : "false" ));
    }
#if MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

// Return the minimum energy
// from the master processor
template<typename energy_t,
         typename logdos_t, 
         typename obs_t, 
         class histogram_index_functor>
REWL_Walker<energy_t, 
            logdos_t, 
            obs_t, 
            histogram_index_functor>::
            get_min_energy()
{
    energy_t min_energy = system.current_state.energy;
#if MPI_ON
    // Broadcast the master processor min energy to everyone
    MPI_Bcast( &min_energy, 1, MPI_ENERGY_TYPE, REWL_MASTER_PROC, MPI_COMM_WORLD );
#endif
    return min_energy;
}

// REWL_Walker constructor
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
    // Don't do anything here...
}

// Reinitialize the walker's binning
// and set the state.
template<typename energy_t,
         typename logdos_t, 
         typename obs_t, 
         class histogram_index_functor>
REWL_Walker<energy_t, 
            logdos_t, 
            obs_t, 
            histogram_index_functor>::
            reinitialize_energies( const energy_t min, const energy_t max,
                                   const energy_t binsize, const size_t num_bins )
{
   hist_idx = histogram_index_functor(min, max, binsize);
   system_obs = Observables_t<obs_t>(num_bins);
   wl_walker = Wang_Landau<energy_t, logdos_t, Hamiltonian_t<obs_t>, Observables_t<obs_t>, State_t<obs_t>, histogram_index_functor>(min, max, binsize, num_bins);

   // Finally, readjust the state to the proper range.
   adjust_state_to_range();
}

template<typename energy_t, typename logdos_t, typename obs_t, class histogram_index_functor>
void REWL_Walker<energy_t, logdos_t, obs_t, histogram_index_functor>::wang_landau_walk(const size_t num_sweeps
#if SAMPLE_AFTER
        , const bool sample_observables 
#endif
    )
{
    size_t system_size = System_Parameters::N;
    size_t num_flavors = System_Parameters::num_DoF / system_size;
    for ( size_t sweep = 0; sweep != num_sweeps; ++sweep )
    {
        wl_walker.wang_landau_sweep(system_size, num_flavors,
#if SAMPLE_AFTER
                sample_observables,
#endif 
                incrementer, &system, &system_obs, random, hist_idx);    
    }
}

/* ******************************************************** */
/* Required REWL functions go here                          */
/* ******************************************************** */

// Wrapper around the histogram indexing functor
template<typename energy_t, typename logdos_t, typename obs_t, class histogram_index_functor>
bool REWL_Walker<energy_t, logdos_t, obs_t, histogram_index_functor>::energy_in_range( const energy_t energy ) const
{
    return hist_idx.energy_in_range( energy );
}

// Wrapper around the rng
template<typename energy_t, typename logdos_t, typename obs_t, class histogram_index_functor>
float REWL_Walker<energy_t, logdos_t, obs_t, histogram_index_functor>::get_rand()
{
    return random();
}

// Wrapper around the rewl histograms. Get the value of the
// logDoS at a particular energy
template<typename energy_t, typename logdos_t, typename obs_t, class histogram_index_functor>
logdos_t REWL_Walker<energy_t, logdos_t, obs_t, histogram_index_functor>::get_logdos( const energy_t energy ) const
{
    return wl_walker.wl_histograms.get_logdos( hist_idx(energy) );
}

// Wrapper for the current energy
template<typename energy_t, typename logdos_t, typename obs_t, class histogram_index_functor>
energy_t REWL_Walker<energy_t, logdos_t, obs_t, histogram_index_functor>::current_energy() const
{
    return system.current_state.energy;
}

// Get the pointer to the current state
template<typename energy_t, typename logdos_t, typename obs_t, class histogram_index_functor>
State_t<obs_t> * REWL_Walker<energy_t, logdos_t, obs_t, histogram_index_functor>::current_state() 
{
    return &( system.current_state );
}

// Get the pointer to the front of the degrees of freedom
template<typename energy_t, typename logdos_t, typename obs_t, class histogram_index_functor>
obs_t * REWL_Walker<energy_t, logdos_t, obs_t, histogram_index_functor>::DoFs() const
{
    return system.get_front_DoFs();
}

// Wrapper around the update histograms functionality
template<typename energy_t, typename logdos_t, typename obs_t, class histogram_index_functor>
void REWL_Walker<energy_t, logdos_t, obs_t, histogram_index_functor>::update_histograms()
{
    size_t bin = hist_idx(system.current_state.energy);
    wl_walker.wl_histograms.increment_count( bin );
    wl_walker.wl_histograms.increment_logdos( bin, incrementer );
}

// Wrapper around the update observables
template<typename energy_t, typename logdos_t, typename obs_t, class histogram_index_functor>
#if SAMPLE_AFTER
void REWL_Walker<energy_t, logdos_t, obs_t, histogram_index_functor>::update_observables( const bool sample_observables )
{
    if ( sample_observables )
    {
        size_t bin = hist_idx(system.current_state.energy);
        system.update_observables( bin, &system_obs );
    }
}
#else
void REWL_Walker<energy_t, logdos_t, obs_t, histogram_index_functor>::update_observables()
{
    size_t bin = hist_idx(system.current_state.energy);
    system.update_observables( bin, &system_obs );
}
#endif


#endif
