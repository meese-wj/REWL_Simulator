#ifndef ASKIN_TELLER2D_HAMILTONIAN
#define ASKIN_TELLER2D_HAMILTONIAN

// Include the parameters 
#include "ashkin_teller2d_parameters.cxx"
#include "ashkin_teller2d_parameter_string.hpp"

// Include the observables enum class
#include "ashkin_teller2d_observables.hpp"

// Include the grid setup from the 
// cmake include directories.
#include <grid_setup.hpp>

#if CORRELATION_LENGTHS
// Include the correlation functionality.
#include "../Correlations/fourier_correlator.hpp"
#endif

#if RFAT_BAXTER
// Include the random fields
#include "../Disorder/random_fields.hpp"
#if MPI_ON 
#include <mpi.h>
#endif
#endif

#if AT_DENSITIES
#include "Density_Plots/ashkin_teller_densities.hpp"
#endif

#if SIMULATED_ANNEALING
#include <random_number_generators.hpp>
#endif

constexpr float divisor = 1./static_cast<float>( Ashkin_Teller2d_Parameters::N );

enum spin_type
{
    sigma, tau, NUM_SPIN_TYPES
};

// TODO: Upgrade the energies to allow for 
// double calculations. 
template<typename data_t>
struct State
{
    short which_to_update = spin_type::sigma;      
    float energy = 0;
    data_t sigma_magnetization = 0;
    data_t tau_magnetization = 0;
    data_t baxter = 0.;
    data_t DoF [spin_type::NUM_SPIN_TYPES] = { 1., 1. }; // Store the local degree of freedom
};

template<typename data_t>
void print(const State<data_t> & stat)
{
    printf("\nState:");
    printf("\n\twhich to update     = %s", stat.which_to_update == spin_type::sigma ? "sigma" : "tau");
    printf("\n\tenergy              = %e", stat.energy);
    printf("\n\tsigma magnetization = %e", stat.sigma_magnetization);
    printf("\n\ttau magnetization   = %e", stat.tau_magnetization);
    printf("\n\tbaxter          = %e", stat.baxter);
    printf("\n\tDoF                 = %e\n", stat.DoF);
}

template<typename data_t>
struct Ashkin_Teller2d
{
    State<data_t> current_state;

    data_t * spin_array = nullptr;
#if RFAT_BAXTER
    float * field_array = nullptr;
#endif
    size_t * neighbor_array = nullptr;

    // Add some Hamiltonian dependent functions
    float local_energy(const size_t idx, const data_t * const spins) const;
    void recalculate_state();
    void change_state(const size_t idx, State<data_t> & temp_state);
    void set_state(const size_t idx, const State<data_t> & _state);
    void update_observables(const size_t bin, Ashkin_Teller2d_Obs<data_t> * obs_ptr) const;

    // Finally add the constructor and destructor.
    Ashkin_Teller2d()
    {
        spin_array = new data_t [ static_cast<size_t>(spin_type::NUM_SPIN_TYPES) * Ashkin_Teller2d_Parameters::N ];
        // TODO: change this to be randomized?
        for ( size_t idx = 0; idx != static_cast<size_t>(spin_type::NUM_SPIN_TYPES) * Ashkin_Teller2d_Parameters::N; ++idx )
            spin_array[idx] = 1.;

        // TODO: Generalize this to different types 
        // of grids.
        // Allocate the neighbor array
        define_2d_square_periodic_neighbors(Ashkin_Teller2d_Parameters::L,
                                            Ashkin_Teller2d_Parameters::L,
                                            Ashkin_Teller2d_Parameters::num_neighbors_i,
                                            neighbor_array);
#if RFAT_BAXTER
        generate_random_field<float>( Ashkin_Teller2d_Parameters::N, 
                                      Ashkin_Teller2d_Parameters::h, 
                                      field_array, disorder_distribution::uniform );
#endif
 
        recalculate_state();
    }

    // Get sigma or tau at a site
    data_t * spin_at_site(const size_t site, const spin_type type) const {   return &( spin_array[ static_cast<size_t>(spin_type::NUM_SPIN_TYPES) * site + static_cast<size_t>(type)] ); }
    data_t * spins_address(const size_t site) const {   return &( spin_array[ static_cast<size_t>(spin_type::NUM_SPIN_TYPES) * site ] ); }

    ~Ashkin_Teller2d()
    {
        delete [] spin_array;
        delete [] neighbor_array;
#if RFAT_BAXTER
        delete [] field_array;
#endif
    }

    void switch_flavor_to_update() 
    {
        if ( current_state.which_to_update == spin_type::sigma )
            current_state.which_to_update = spin_type::tau;
        else
            current_state.which_to_update = spin_type::sigma;
    }
    void print_lattice() const;

    data_t * get_front_DoFs() const { return spin_array; }

    void import_DoFs( const data_t * const array )
    {
        for ( size_t idx = 0; idx != Ashkin_Teller2d_Parameters::num_DoF; ++idx )
            spin_array[idx] = array[idx];

        recalculate_state();
    }

#if SIMULATED_ANNEALING
    void randomize_dofs() 
    {
        std::uint64_t seed = static_cast<std::uint64_t>( std::chrono::high_resolution_clock::now().time_since_epoch().count() );
        random_number_generator<float> rng (seed);

        for ( size_t idx = 0; idx != Ashkin_Teller2d_Parameters::num_DoF; ++idx )
            spin_array[ idx ] = ( rng() < 0.5 ? 1. : -1. );

        recalculate_state();
    }
#endif

#if RFAT_BAXTER
    void import_disorder( const float * const disorder )
    {
        for ( size_t idx = 0; idx != Ashkin_Teller2d_Parameters::N; ++idx )
            field_array[idx] = disorder[idx];

        recalculate_state();
    }
#endif
};

template<typename data_t>
float Ashkin_Teller2d<data_t>::local_energy(const size_t idx, const data_t * const spins ) const
{
    // Calculate positive energies everywhere because
    // it will be negated upon return!
    float en = 0.;
    data_t sigma_idx = *(spins + spin_type::sigma);
    data_t tau_idx   = *(spins + spin_type::tau);
    size_t neighbor = 0;
    for ( size_t nidx = 0; nidx != Ashkin_Teller2d_Parameters::num_neighbors_i; ++nidx )
    {
        neighbor = neighbor_array[ idx * Ashkin_Teller2d_Parameters::num_neighbors_i + nidx ];
        en += Ashkin_Teller2d_Parameters::J * static_cast<float>( sigma_idx * (*spin_at_site(neighbor, spin_type::sigma)) );
        en += Ashkin_Teller2d_Parameters::J * static_cast<float>( tau_idx   * (*spin_at_site(neighbor, spin_type::tau)  ) );
        en += Ashkin_Teller2d_Parameters::K * static_cast<float>( sigma_idx * tau_idx
                                                                * (*spin_at_site(neighbor, spin_type::sigma))
                                                                * (*spin_at_site(neighbor, spin_type::tau)  )   );

    }

#if RFAT_BAXTER
    en += field_array[idx] * static_cast<float>( sigma_idx * tau_idx );
#endif
    return -en;
}

template<typename data_t>
void Ashkin_Teller2d<data_t>::recalculate_state()
{
    float temp_energy      = 0.;
    data_t temp_sigma_mag  = 0.;
    data_t temp_tau_mag    = 0.;
    data_t temp_baxter = 0.;

    for ( size_t idx = 0; idx != Ashkin_Teller2d_Parameters::N; ++idx )
    {
       temp_sigma_mag += *spin_at_site(idx, spin_type::sigma);
       temp_tau_mag   += *spin_at_site(idx, spin_type::tau);
       temp_baxter += (*spin_at_site(idx, spin_type::sigma)) * (*spin_at_site(idx, spin_type::tau));

#if RFAT_BAXTER
       temp_energy += -0.5 * field_array[idx] * static_cast<float>( (*spin_at_site(idx, spin_type::sigma)) * (*spin_at_site(idx, spin_type::tau)) );
#endif
       temp_energy += 0.5 * local_energy(idx, spins_address(idx));
    }

    current_state.energy              = temp_energy;
    current_state.sigma_magnetization = temp_sigma_mag;
    current_state.tau_magnetization   = temp_tau_mag;
    current_state.baxter              = temp_baxter;
}

// Change the state by changing the value of the spin
// at a particular site idx. Return what the new state
// would be if the change is accepted (the state is
// passed by reference to avoid extra copies).
template<typename data_t>
void Ashkin_Teller2d<data_t>::change_state(const size_t idx, State<data_t> & temp_state)
{
    temp_state.DoF[spin_type::sigma] = *spin_at_site(idx, spin_type::sigma);
    temp_state.DoF[spin_type::tau]   = *spin_at_site(idx, spin_type::tau);

    // Change sigma and tau only if they're supposed 
    // to according to the current state. 
    temp_state.DoF[spin_type::sigma] += -2. * temp_state.DoF[spin_type::sigma] * ( current_state.which_to_update == spin_type::sigma );
    temp_state.DoF[spin_type::tau]   += -2. * temp_state.DoF[spin_type::tau] * ( current_state.which_to_update == spin_type::tau );

    // Store which spin was updated in case they 
    // need to switch flavors.
    temp_state.which_to_update =   spin_type::sigma * ( current_state.which_to_update == spin_type::tau ) 
                                 + spin_type::tau * ( current_state.which_to_update == spin_type::sigma );
   

    // Update the state now
    temp_state.sigma_magnetization = current_state.sigma_magnetization + temp_state.DoF[spin_type::sigma] - (*spin_at_site(idx, spin_type::sigma));
    temp_state.tau_magnetization   = current_state.tau_magnetization   + temp_state.DoF[spin_type::tau]   - (*spin_at_site(idx, spin_type::tau));
    temp_state.baxter              = current_state.baxter + (temp_state.DoF[spin_type::sigma] * temp_state.DoF[spin_type::tau]) - ( *spin_at_site(idx, spin_type::sigma) * (*spin_at_site(idx, spin_type::tau)) );

    temp_state.energy              = current_state.energy + local_energy(idx, &(temp_state.DoF[0])) - local_energy(idx, spins_address(idx) );

    // After a single sweep, switch the spin type being updated
    current_state.which_to_update =   current_state.which_to_update * ( idx != Ashkin_Teller2d_Parameters::N - 1 )
                                    + temp_state.which_to_update    * ( idx == Ashkin_Teller2d_Parameters::N - 1 ); 
}

// Set the current state
template<typename data_t>
void Ashkin_Teller2d<data_t>::set_state(const size_t idx, const State<data_t> & _state)
{
    current_state.energy                 = _state.energy;

    current_state.sigma_magnetization    = _state.sigma_magnetization;
    current_state.tau_magnetization      = _state.tau_magnetization;
    current_state.baxter                 = _state.baxter;

    *spin_at_site(idx, spin_type::sigma) = _state.DoF[spin_type::sigma];
    *spin_at_site(idx, spin_type::tau)   = _state.DoF[spin_type::tau];
}

// Update the non-energetic observables
template<typename data_t>
void Ashkin_Teller2d<data_t>::update_observables(const size_t bin, Ashkin_Teller2d_Obs<data_t> * obs_ptr) const
{
    const data_t sigma_val  = Ashkin_Teller2d_Parameters::divide_N * abs( current_state.sigma_magnetization );
    const data_t tau_val    = Ashkin_Teller2d_Parameters::divide_N * abs( current_state.tau_magnetization );
    const data_t baxter_val = Ashkin_Teller2d_Parameters::divide_N * abs( current_state.baxter );
    const data_t phi_val2   = sigma_val * sigma_val + tau_val * tau_val;
    
    obs_ptr -> update_observable_average(sigma_val, Obs::enum_names::sigma_mag, bin);
    obs_ptr -> update_observable_average(sigma_val * sigma_val, Obs::enum_names::sigma_mag2, bin);
    obs_ptr -> update_observable_average(sigma_val * sigma_val * sigma_val * sigma_val, Obs::enum_names::sigma_mag4, bin);

    obs_ptr -> update_observable_average(tau_val, Obs::enum_names::tau_mag, bin);
    obs_ptr -> update_observable_average(tau_val * tau_val, Obs::enum_names::tau_mag2, bin);
    obs_ptr -> update_observable_average(tau_val * tau_val * tau_val * tau_val, Obs::enum_names::tau_mag4, bin);

    obs_ptr -> update_observable_average(sqrt(phi_val2), Obs::enum_names::phi, bin);
    obs_ptr -> update_observable_average(phi_val2,  Obs::enum_names::phi2, bin);
    obs_ptr -> update_observable_average(phi_val2 * phi_val2, Obs::enum_names::phi4, bin);

    obs_ptr -> update_observable_average(baxter_val, Obs::enum_names::baxter_mag, bin);
    obs_ptr -> update_observable_average(baxter_val * baxter_val, Obs::enum_names::baxter_mag2, bin);
    obs_ptr -> update_observable_average(baxter_val * baxter_val * baxter_val * baxter_val, Obs::enum_names::baxter_mag4, bin);

#if AT_DENSITIES
    update_densities( obs_ptr -> density_histograms, bin, Ashkin_Teller2d_Parameters::divide_N * current_state.sigma_magnetization, Ashkin_Teller2d_Parameters::divide_N * current_state.tau_magnetization );
#endif

#if CORRELATION_LENGTHS
    if ( static_cast<size_t> (obs_ptr -> get_observable(Obs::enum_names::counts_per_bin, bin)) % counts_per_transform == 0 )
    {
        data_t sigma_Gq_value = Ashkin_Teller2d_Parameters::divide_N * Ashkin_Teller2d_Parameters::divide_N * ( obs_ptr -> correlator.compute_correlator( spin_array, spin_type::sigma, spin_type::tau, spin_type::NUM_SPIN_TYPES, false ) );
        data_t tau_Gq_value =   Ashkin_Teller2d_Parameters::divide_N * Ashkin_Teller2d_Parameters::divide_N * ( obs_ptr -> correlator.compute_correlator( spin_array, spin_type::sigma, spin_type::tau, spin_type::NUM_SPIN_TYPES, false ) );
        data_t baxter_Gq_value = Ashkin_Teller2d_Parameters::divide_N * Ashkin_Teller2d_Parameters::divide_N * ( obs_ptr -> correlator.compute_correlator( spin_array, spin_type::sigma, spin_type::tau, spin_type::NUM_SPIN_TYPES, true ) );

        // Sigma correlator
        obs_ptr -> update_qmin_correlator( sigma_Gq_value, Obs::enum_names::sigma_corr_qmin, bin,
                                           obs_ptr -> get_observable(Obs::enum_names::counts_per_bin, bin) / counts_per_transform );
        // Tau correlator
        obs_ptr -> update_qmin_correlator( tau_Gq_value, Obs::enum_names::tau_corr_qmin, bin,
                                           obs_ptr -> get_observable(Obs::enum_names::counts_per_bin, bin) / counts_per_transform );
        // Phi correlator
        obs_ptr -> update_qmin_correlator( sigma_Gq_value + tau_Gq_value, Obs::enum_names::phi_corr_qmin, bin,
                                           obs_ptr -> get_observable(Obs::enum_names::counts_per_bin, bin) / counts_per_transform );
        // Baxter correlator
        obs_ptr -> update_qmin_correlator( baxter_Gq_value, Obs::enum_names::baxter_corr_qmin, bin,
                                           obs_ptr -> get_observable(Obs::enum_names::counts_per_bin, bin) / counts_per_transform );
    }
#endif
    
    obs_ptr -> increment_counts_per_bin(bin);
}

template<typename data_t>
void Ashkin_Teller2d<data_t>::print_lattice() const
{
    printf("\n\nSigmas +/-");
    for ( size_t idx = 0; idx != Ashkin_Teller2d_Parameters::L; ++idx )
    {
        printf("\n%ld\t    ", idx * Ashkin_Teller2d_Parameters::L);
        for ( size_t jdx = 0; jdx != Ashkin_Teller2d_Parameters::L; ++jdx )
        {
            printf("%c ", *spin_at_site( idx * Ashkin_Teller2d_Parameters::L + jdx, spin_type::sigma ) == 1. ? '+' : '-');   
        }
    }
    printf("\nTaus */o");
    for ( size_t idx = 0; idx != Ashkin_Teller2d_Parameters::L; ++idx )
    {
        printf("\n%ld\t    ", idx * Ashkin_Teller2d_Parameters::L);
        for ( size_t jdx = 0; jdx != Ashkin_Teller2d_Parameters::L; ++jdx )
        {
            printf("%c ", *spin_at_site( idx * Ashkin_Teller2d_Parameters::L + jdx, spin_type::tau ) == 1. ? '*' : 'o');   
        }
    }
}

#if MPI_ON
#ifndef INDEPENDENT_WALKERS
#include <rewl_parameters.hpp>   // Include this for my MPI datatype macros

enum State_SendRecv
{
    energy = 444, sigma_mag, tau_mag, bax, state_dof, all_dof
};

// TODO: This should be specialized because it will break probably
template<class state_t>
void mpi_exchange_state( state_t * const state, const int partner_index, const int comm_id, const MPI_Comm * const local_communicators, MPI_Status * const status  )
{
    MPI_Sendrecv_replace( &( state -> energy  ), 1, MPI_ENERGY_TYPE, partner_index, State_SendRecv::energy, partner_index, State_SendRecv::energy, local_communicators[ comm_id  ], status  );

    MPI_Sendrecv_replace( &( state -> sigma_magnetization  ), 1, MPI_OBS_TYPE, partner_index, State_SendRecv::sigma_mag, partner_index, State_SendRecv::sigma_mag, local_communicators[ comm_id  ], status  );

    MPI_Sendrecv_replace( &( state -> tau_magnetization  ), 1, MPI_OBS_TYPE, partner_index, State_SendRecv::tau_mag, partner_index, State_SendRecv::tau_mag, local_communicators[ comm_id  ], status  );
    
    MPI_Sendrecv_replace( &( state -> baxter ), 1, MPI_OBS_TYPE, partner_index, State_SendRecv::bax, partner_index, State_SendRecv::bax, local_communicators[ comm_id  ], status  );
    
    MPI_Sendrecv_replace( &( state -> DoF ), static_cast<int>(spin_type::NUM_SPIN_TYPES), MPI_OBS_TYPE, partner_index, State_SendRecv::state_dof, partner_index, State_SendRecv::state_dof, local_communicators[ comm_id  ], status  );

}

template<typename obs_t>
void mpi_exchange_DoFs( obs_t * const front_dofs, const size_t num_dof, const int partner_index, const int comm_id, const MPI_Comm * const local_communicators, MPI_Status * const status  );

// Specialize this templated functions
template<>
void mpi_exchange_DoFs<OBS_TYPE>( OBS_TYPE * const front_dofs, const size_t num_dof, const int partner_index, const int comm_id, const MPI_Comm * const local_communicators, MPI_Status * const status  )
{
    MPI_Sendrecv_replace( front_dofs, static_cast<int>(num_dof), MPI_OBS_TYPE, partner_index, State_SendRecv::all_dof, partner_index, State_SendRecv::all_dof, local_communicators[ comm_id  ], status  );
}

#endif
#endif

#endif


