#ifndef WANG_LANDAU
/* Use this file to define the Wang Landau
 * operations for a single random walker.
 * These functions will be used by every
 * MPI process running in the Replica-
 * Exchange Wang Landau algorithm. */

#include <cmath>
#include <random_number_generator.hpp>
#include "../rewl_histograms.hpp"

template<typename data_t>
using rng = random_number_generator<data_t>;

template<typename data_t, class Hamiltonian_t, class Observables_t, class State_t>
struct Wang_Landau
{
    rewl_histograms<data_t> wl_histograms;
    
    Wang_Landau(const data_t _min, const data_t _max, const data_t _bsize) : wl_histograms(_min, _max, _bsize) 
    {}

    ~Wang_Landau() {}

    void wang_landau_update(const size_t idx, const data_t incrementer, Hamiltonian_t * ham,
                            Observables_t * ham_obs, rng<float> & random, const histogram_index & hist_idx) const;

    void wang_landau_sweep(const size_t system_size, const data_t incrementer, Hamiltonian_t * ham, 
                           Observables_t * ham_obs, rng<float> & random, const histogram_index & hist_idx ) const;

    bool is_flat(const float tolerance) const;
};

// Run a single update step using the
// normal Wang Landau method
template<typename data_t, class Hamiltonian_t, class Observables_t, class State_t>
void Wang_Landau<data_t, Hamiltonian_t, Observables_t, State_t>::wang_landau_update(const size_t idx, const data_t incrementer,
                                                                                    Hamiltonian_t * ham, Observables_t * ham_obs, 
                                                                                    rng<float> & random, const histogram_index & hist_idx) const
{
    State_t temporary_state;
    change_state(idx, temporary_state);
    
    const size_t new_bin = hist_idx(temporary_state.energy);
    size_t current_bin = hist_idx(ham -> current_state.energy);

    float entropy_change =   wl_histograms.get_logdos( new_bin )
                           - wl_histograms.get_logdos( current_bin );

    // Check if the move is allowed
    if ( entropy_change < 0. || random() < exp(-entropy_change) )
    {
        // If allowed, change the state to the temporary one
        // and change the current bin to the new bin. Note that
        // the state MUST contain the new degree of freedom in it.
        ham -> set_state(idx, temporary_state);
        current_bin = new_bin;
    }
   
    // Increment the REWL histograms at the final bin
    // and update the observables.
    wl_histograms -> increment_count( current_bin );
    wl_histograms -> increment_logdos( current_bin, incrementer );
    ham -> update_observables( current_bin, ham_obs );
}

// Run a sweep using the standard
// Wang Landau method on the system
template<typename data_t, class Hamiltonian_t, class Observables_t, class State_t>
void Wang_Landau<data_t, Hamiltonian_t, Observables_t, State_t>::wang_landau_sweep(const size_t system_size, const data_t incrementer,
                                                                                   Hamiltonian_t * ham, Observables_t * ham_obs, 
                                                                                   rng<float> & random, const histogram_index & hist_idx) const
{
    for ( size_t dof_idx = 0; dof_idx != system_size; ++dof_idx )
    {
        wang_landau_update(dof_idx, incrementer, ham, ham_obs, random, hist_idx);
    }
}

#endif
