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
    data_t nematicity = 0.;
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
    printf("\n\tnematicity          = %e", stat.nematicity);
    printf("\n\tDoF                 = %e\n", stat.DoF);
}

template<typename data_t>
struct Ashkin_Teller2d
{
    State<data_t> current_state;

    data_t * spin_array = nullptr;
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

        recalculate_state();
    }

    // Get sigma or tau at a site
    data_t * spin_at_site(const size_t site, const spin_type type) const {   return &( spin_array[ static_cast<size_t>(spin_type::NUM_SPIN_TYPES) * site + static_cast<size_t>(type)] ); }
    data_t * spins_address(const size_t site) const {   return &( spin_array[ static_cast<size_t>(spin_type::NUM_SPIN_TYPES) * site ] ); }

    ~Ashkin_Teller2d()
    {
        delete [] spin_array;
        delete [] neighbor_array;
    }

    void print_lattice() const;
};

template<typename data_t>
float Ashkin_Teller2d<data_t>::local_energy(const size_t idx, const data_t * const spins ) const
{
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
    return -en;
}

template<typename data_t>
void Ashkin_Teller2d<data_t>::recalculate_state()
{
    float temp_energy      = 0.;
    data_t temp_sigma_mag  = 0.;
    data_t temp_tau_mag    = 0.;
    data_t temp_nematicity = 0.;

    for ( size_t idx = 0; idx != Ashkin_Teller2d_Parameters::N; ++idx )
    {
       temp_sigma_mag += *spin_at_site(idx, spin_type::sigma);
       temp_tau_mag   += *spin_at_site(idx, spin_type::tau);
       temp_nematicity += (*spin_at_site(idx, spin_type::sigma)) * (*spin_at_site(idx, spin_type::tau));
       temp_energy += 0.5 * local_energy(idx, spins_address(idx));
    }

    current_state.energy              = temp_energy;
    current_state.sigma_magnetization = temp_sigma_mag;
    current_state.tau_magnetization   = temp_tau_mag;
    current_state.nematicity          = temp_nematicity;
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

    switch (current_state.which_to_update)
    {
        case spin_type::sigma:
        {
            temp_state.DoF[spin_type::sigma] *= -1.;
            temp_state.which_to_update = spin_type::tau;    // Which spin to switch to IFF a switch is made    
            break;
        }
        case spin_type::tau:
        {
            temp_state.DoF[spin_type::tau]   *= -1.;
            temp_state.which_to_update = spin_type::sigma;  // Which spin to switch to IFF a switch is made
            break;
        }
    }

    temp_state.sigma_magnetization = current_state.sigma_magnetization + temp_state.DoF[spin_type::sigma] - (*spin_at_site(idx, spin_type::sigma));
    temp_state.tau_magnetization   = current_state.tau_magnetization   + temp_state.DoF[spin_type::tau]   - (*spin_at_site(idx, spin_type::tau));
    temp_state.nematicity          = current_state.nematicity + (temp_state.DoF[spin_type::sigma] * temp_state.DoF[spin_type::tau]) - ( *spin_at_site(idx, spin_type::sigma) * (*spin_at_site(idx, spin_type::tau)) );

    temp_state.energy              = current_state.energy + local_energy(idx, &(temp_state.DoF[0])) - local_energy(idx, spins_address(idx) );

    // After a single sweep, switch the spin type being updated
    if ( idx == Ashkin_Teller2d_Parameters::N - 1 )
        current_state.which_to_update = temp_state.which_to_update;
        
}

// Set the current state
template<typename data_t>
void Ashkin_Teller2d<data_t>::set_state(const size_t idx, const State<data_t> & _state)
{
    current_state.energy                 = _state.energy;

    current_state.sigma_magnetization    = _state.sigma_magnetization;
    current_state.tau_magnetization      = _state.tau_magnetization;
    current_state.nematicity             = _state.nematicity;

    *spin_at_site(idx, spin_type::sigma) = _state.DoF[spin_type::sigma];
    *spin_at_site(idx, spin_type::tau)   = _state.DoF[spin_type::tau];
}

// Update the non-energetic observables
template<typename data_t>
void Ashkin_Teller2d<data_t>::update_observables(const size_t bin, Ashkin_Teller2d_Obs<data_t> * obs_ptr) const
{
    const data_t sigma_val  = abs( current_state.sigma_magnetization );
    const data_t tau_val    = abs( current_state.tau_magnetization );
    const data_t order_val2 = sigma_val * sigma_val + tau_val * tau_val;
    const data_t nem_val    = abs( current_state.nematicity );
    
    obs_ptr -> update_observable_average(sigma_val, Obs::enum_names::sigma_mag, bin);
    obs_ptr -> update_observable_average(sigma_val * sigma_val, Obs::enum_names::sigma_mag2, bin);
    obs_ptr -> update_observable_average(sigma_val * sigma_val * sigma_val * sigma_val, Obs::enum_names::sigma_mag4, bin);

    obs_ptr -> update_observable_average(tau_val, Obs::enum_names::tau_mag, bin);
    obs_ptr -> update_observable_average(tau_val * tau_val, Obs::enum_names::tau_mag2, bin);
    obs_ptr -> update_observable_average(tau_val * tau_val * tau_val * tau_val, Obs::enum_names::tau_mag4, bin);

    obs_ptr -> update_observable_average(sqrt(order_val2), Obs::enum_names::order_param, bin);
    obs_ptr -> update_observable_average(order_val2,  Obs::enum_names::order_param2, bin);
    obs_ptr -> update_observable_average(order_val2 * order_val2, Obs::enum_names::order_param4, bin);

    obs_ptr -> update_observable_average(nem_val, Obs::enum_names::nem_mag, bin);
    obs_ptr -> update_observable_average(nem_val * nem_val, Obs::enum_names::nem_mag2, bin);
    obs_ptr -> update_observable_average(nem_val * nem_val * nem_val * nem_val, Obs::enum_names::nem_mag4, bin);
    
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
    energy = 444, sigma_mag, tau_mag, nem, state_dof, all_dof
};

// TODO: This should be specialized because it will break probably
template<class state_t>
void mpi_exchange_state( state_t * const state, const int partner_index, const int comm_id, const MPI_Comm * const local_communicators, MPI_Status * const status  )
{
    MPI_Sendrecv_replace( &( state -> energy  ), 1, MPI_ENERGY_TYPE, partner_index, State_SendRecv::energy, partner_index, State_SendRecv::energy, local_communicators[ comm_id  ], status  );

    MPI_Sendrecv_replace( &( state -> sigma_magnetization  ), 1, MPI_OBS_TYPE, partner_index, State_SendRecv::sigma_mag, partner_index, State_SendRecv::sigma_mag, local_communicators[ comm_id  ], status  );

    MPI_Sendrecv_replace( &( state -> tau_magnetization  ), 1, MPI_OBS_TYPE, partner_index, State_SendRecv::tau_mag, partner_index, State_SendRecv::tau_mag, local_communicators[ comm_id  ], status  );
    
    MPI_Sendrecv_replace( &( state -> nematicity  ), 1, MPI_OBS_TYPE, partner_index, State_SendRecv::nem, partner_index, State_SendRecv::nem, local_communicators[ comm_id  ], status  );
    
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


