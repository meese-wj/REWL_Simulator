#ifndef WANG_LANDAU
#define WANG_LANDAU
/* Use this file to define the Wang Landau
 * operations for a single random walker.
 * These functions will be used by every
 * MPI process running in the Replica-
 * Exchange Wang Landau algorithm. */

#include <cmath>
#include <random_number_generators.hpp>
#include "rewl_histograms.hpp"

template<typename data_t>
using rng = random_number_generator<data_t>;

template<typename energy_t, typename logdos_t, class Hamiltonian_t, 
         class Observables_t, class State_t, class histogram_index_functor>
struct Wang_Landau
{
    rewl_histograms<logdos_t> wl_histograms;
    
    Wang_Landau(const energy_t _min, const energy_t _max, const energy_t _bsize, const size_t _nbins) : wl_histograms(_min, _max, _bsize, _nbins) 
    {}

    ~Wang_Landau() {}

#if SAMPLE_AFTER
    void wang_landau_update(const size_t idx, const bool sample_observables, const logdos_t incrementer, Hamiltonian_t * const ham,
                            Observables_t * const ham_obs, rng<double> & random, const histogram_index_functor & hist_idx);

    void wang_landau_sweep(const size_t system_size, const size_t num_flavors, const bool sample_observables, const logdos_t incrementer, Hamiltonian_t * const ham, 
                           Observables_t * const ham_obs, rng<double> & random, const histogram_index_functor & hist_idx );
#else
    void wang_landau_update(const size_t idx, const logdos_t incrementer, Hamiltonian_t * const ham,
                            Observables_t * const ham_obs, rng<double> & random, const histogram_index_functor & hist_idx);

    void wang_landau_sweep(const size_t system_size, const size_t num_flavors, const logdos_t incrementer, Hamiltonian_t * const ham, 
                           Observables_t * const ham_obs, rng<double> & random, const histogram_index_functor & hist_idx );
#endif


    bool is_flat(const float tolerance) const;
    void reset_histogram() const;
};

// Run a single update step using the
// normal Wang Landau method
template<typename energy_t, typename logdos_t, class Hamiltonian_t, 
         class Observables_t, class State_t, class histogram_index_functor>
void Wang_Landau<energy_t, logdos_t, Hamiltonian_t, Observables_t, State_t, histogram_index_functor>::wang_landau_update(const size_t idx, 
#if SAMPLE_AFTER
        const bool sample_observables,
#endif
        const logdos_t incrementer,
                                                                                    Hamiltonian_t * const ham, Observables_t * const ham_obs, 
                                                                                    rng<double> & random, const histogram_index_functor & hist_idx)
{
    //printf("\n\nSite %ld of %ld: spin = %e", idx, System_Parameters::N-1, ham -> spin_array[idx]);
    //printf("\n\nDoF %ld = %c before", idx, ham -> spin_array[idx] == 1. ? '+' : '-');
    //ham -> print_lattice();
    State_t temporary_state;
    //printf("\n\ncurrent state");
    //print(ham -> current_state);
    ham -> change_state(idx, temporary_state);
    //printf("\n\ntemporary state");
    //print(temporary_state);
    
    size_t current_bin = hist_idx(ham -> current_state.energy);

    /*
    if ( current_bin >= System_Parameters::num_bins -5 )
    {
        printf("\n\nSite %ld of %ld: spin = %e", idx, System_Parameters::N-1, ham -> spin_array[idx]);
        printf("\n\nDoF %ld = %c before", idx, ham -> spin_array[idx] == 1. ? '+' : '-');
        ham -> print_lattice();
        printf("\n\ncurrent state");
        print(ham -> current_state);
        printf("\n\ntemporary state");
        print(temporary_state);
    }
    */

    if ( hist_idx.energy_in_range( temporary_state.energy ) )
    {
        const size_t new_bin = hist_idx(temporary_state.energy);

        float entropy_change =   wl_histograms.get_logdos( new_bin )
                               - wl_histograms.get_logdos( current_bin );

        //printf("\ncurrent energy, new energy = %e, %e", ham -> current_state.energy, temporary_state.energy);
        //printf("\ncurrent bin, new bin, entropy change = %ld, %ld, %e", current_bin, new_bin, entropy_change);

        // Check if the move is allowed
        if ( entropy_change <= 0. || random() < static_cast<double>(exp(-entropy_change)) )
        {
            // If allowed, change the state to the temporary one
            // and change the current bin to the new bin. Note that
            // the state MUST contain the new degree of freedom in it.
            //printf("\nChanging State...\n");
            ham -> set_state(idx, temporary_state);
            current_bin = new_bin;
        }
    }
   
    // Increment the REWL histograms at the final bin
    // and update the observables.
    wl_histograms.increment_count( current_bin );
    wl_histograms.increment_logdos( current_bin, incrementer );
#if SAMPLE_AFTER
    if ( sample_observables ) ham -> update_observables( current_bin, ham_obs );
#else
    ham -> update_observables( current_bin, ham_obs );
#endif
  
    /*
    if ( current_bin >= System_Parameters::num_bins -6)
    {
        printf("\n\nDoF %ld = %c after", idx, ham -> spin_array[idx] == 1. ? '+' : '-');
        ham -> print_lattice();

        printf("\n\n");
        printf("\n*************************************************************************\n");
        printf("\n*************************************************************************\n");
        printf("\n*************************************************************************\n");
    }
    */
}


// Run a sweep using the standard
// Wang Landau method on the system
template<typename energy_t, typename logdos_t, class Hamiltonian_t, 
         class Observables_t, class State_t, class histogram_index_functor>
void Wang_Landau<energy_t, logdos_t, Hamiltonian_t, Observables_t, State_t, histogram_index_functor>::wang_landau_sweep(const size_t system_size, const size_t num_flavors, 
#if SAMPLE_AFTER
        const bool sample_observables,
#endif
        const logdos_t incrementer, Hamiltonian_t * const ham, Observables_t * const ham_obs, 
                                                                                   rng<double> & random, const histogram_index_functor & hist_idx)
{
    // Update each flavor of the degree of freedom independently
    // while all others are quenched for a single sweep.
    for ( size_t flavor = 0; flavor != num_flavors; ++flavor )
    {
        for ( size_t dof_idx = 0; dof_idx != system_size; ++dof_idx )
        {
            wang_landau_update(dof_idx, 
#if SAMPLE_AFTER
                    sample_observables,
#endif
                    incrementer, ham, ham_obs, random, hist_idx);
        }
    }
}


template<typename energy_t, typename logdos_t, class Hamiltonian_t, 
         class Observables_t, class State_t, class histogram_index_functor>
bool Wang_Landau<energy_t, logdos_t, Hamiltonian_t, Observables_t, State_t, histogram_index_functor>::is_flat(const float tolerance) const
{
    //printf("\n\nTolerance = %e, histogram flatness = %e", tolerance, wl_histograms.count_flatness());
    return ( wl_histograms.count_flatness() <= tolerance );
}


template<typename energy_t, typename logdos_t, class Hamiltonian_t, 
         class Observables_t, class State_t, class histogram_index_functor>
void Wang_Landau<energy_t, logdos_t, Hamiltonian_t, Observables_t, State_t, histogram_index_functor>::reset_histogram() const
{
    wl_histograms.reset_counts();
}

#endif
