#ifndef ISING2D_HAMILTONIAN
#define ISING2D_HAMILTONIAN

// Include the parameters 
#include "ising2d_parameters.cxx"
#include "ising2d_parameter_string.hpp"

// Include the observables enum class
#include "ising2d_observables.hpp"

// Include the address book
// #include "../Address_Books/Square_Periodic_Lattices/square_2d_nearest_neighbors.hpp" 
#include "../Address_Books/Square_Periodic_Lattices/square_2d_nearest_neighbor_functor.cpp" // need the cpp to avoid linker errors

#if CORRELATION_LENGTHS
// Include the correlation functionality.
#include "../Correlations/fourier_correlator.hpp"
#endif

#if MPI_ON 
#include <mpi.h>
#endif

#if RFIM
// Include the random fields
#include "../Disorder/random_fields.hpp"
   
    // Define the disorder type
    #if UNIFORM_DISORDER
        constexpr disorder_distribution disorder_type = disorder_distribution::uniform;
    #elif GAUSSIAN_DISORDER
        constexpr disorder_distribution disorder_type = disorder_distribution::gaussian;
    #endif

#endif

#if SIMULATED_ANNEALING
#include <random_number_generators.hpp>
#endif

#if PHONON_MEDIATED_NEMATIC_INTERACTIONS
#include "PM_Dipolar_Interactions/pmd_interactions.hpp"
#include "PM_Dipolar_Interactions/pmd_ising2d_observables.hpp"  // TODO: Not sure if this is necessary...
#endif // PHONON_MEDIATED_NEMATIC_INTERACTIONS

// TODO: Upgrade the energies to allow for 
// double calculations. 
template<typename data_t>
struct State
{
    float energy = 0;
    data_t magnetization = 0;
    data_t DoF = 0.;                 // Store the local degree of freedom
};

template<typename data_t>
void print(const State<data_t> & stat)
{
    printf("\nState:");
    printf("\n\tenergy        = %e", stat.energy);
    printf("\n\tmagnetization = %e", stat.magnetization);
    printf("\n\tDoF           = %e\n", stat.DoF);
}

template<typename data_t>
struct Ising2d
{
    State<data_t> current_state;

    data_t * spin_array = nullptr;
#if RFIM
    // TODO: change the energies to doubles.
    float * field_array = nullptr;
#endif
    
    Square_2D_Nearest_Neighbor_Functor<Ising2d_Parameters::L> address_book;

#if PHONON_MEDIATED_NEMATIC_INTERACTIONS
    PMDN_Interactions<float, data_t> * pmd_interaction;
#endif // PHONON_MEDIATED_NEMATIC_INTERACTIONS

    // Add some Hamiltonian dependent functions
    float local_field(const size_t idx) const;
    float local_energy(const size_t idx, const data_t spin_value) const;
    void recalculate_state();
    void change_state(const size_t idx, State<data_t> & temp_state) const;
    void set_state(const size_t idx, const State<data_t> & _state);
    void update_observables(const size_t bin, Ising2d_Obs<data_t> * obs_ptr) const;

    // Finally add the constructor and destructor.
    // Ising2d() : address_book(Ising2d_Parameters::L, Ising2d_Parameters::L)
    Ising2d() : address_book()
    {
        spin_array = new data_t [ Ising2d_Parameters::N ];
        // TODO: change this to be randomized?
        for ( size_t idx = 0; idx != Ising2d_Parameters::N; ++idx )
            spin_array[idx] = 1.;

        address_book.initialize();
        address_book.print();

#if RFIM
        generate_random_field<float>( Ising2d_Parameters::N, Ising2d_Parameters::h, 
                                      field_array, disorder_type );
#endif

#if PHONON_MEDIATED_NEMATIC_INTERACTIONS
        pmd_interaction = new PMDN_Interactions<float, data_t>( Ising2d_Parameters::PMNI_Coupling, Ising2d_Parameters::L, Ising2d_Parameters::L );
#endif // PHONON_MEDIATED_NEMATIC_INTERACTIONS
        recalculate_state();
    }

    ~Ising2d()
    {
        delete [] spin_array;
#if RFIM
        delete [] field_array;
#endif

#if PHONON_MEDIATED_NEMATIC_INTERACTIONS
        delete pmd_interaction;
#endif // PHONON_MEDIATED_NEMATIC_INTERACTIONS
    }

    void switch_flavor_to_update() { return; };

    void print_lattice() const;
    
    data_t * get_front_DoFs() const { return spin_array; }

    void import_DoFs( const data_t * const array )
    {
        for ( size_t idx = 0; idx != Ising2d_Parameters::num_DoF; ++idx )
            spin_array[idx] = array[idx];

        recalculate_state();
    }

#if SIMULATED_ANNEALING
    void randomize_dofs()
    {
        std::uint64_t seed = static_cast<std::uint64_t>( std::chrono::high_resolution_clock::now().time_since_epoch().count() );
        random_number_generator<float> rng (seed);

        for ( size_t idx = 0; idx != Ising2d_Parameters::num_DoF; ++idx )
            spin_array[ idx ] = ( rng() < 0.5 ? 1. : -1. );

        recalculate_state();

        print(current_state);
        print_lattice();
    }
#endif

#if RFIM
    void import_disorder( const float * const disorder )
    {
        for ( size_t idx = 0; idx != Ising2d_Parameters::N; ++idx )
            field_array[idx] = disorder[idx];

        recalculate_state();
    }
#endif
};

template<typename data_t>
float Ising2d<data_t>::local_field(const size_t idx) const
{
#if RFIM
    float field = 0.;
#else
    float field = Ising2d_Parameters::h;
#endif
    for ( auto nn_itr = address_book.neighbor_begin(idx); nn_itr != address_book.neighbor_end(idx); ++nn_itr )
    {
        field += Ising2d_Parameters::J * static_cast<float>( spin_array[ *nn_itr ] );
    }

#if RFIM
    field += field_array[ idx ];
#endif 

    return field;
}

template<typename data_t>
float Ising2d<data_t>::local_energy(const size_t idx, const data_t spin_value) const
{
    float local_en = -1. * static_cast<float>( spin_value ) * local_field(idx);
#if PHONON_MEDIATED_NEMATIC_INTERACTIONS
    local_en += pmd_interaction -> calculate_energy_per_spin( idx, spin_value, spin_array );
#endif // PHONON_MEDIATED_NEMATIC_INTERACTIONS
    return local_en;
}

template<typename data_t>
void Ising2d<data_t>::recalculate_state()
{
    float temp_energy = 0.;
    data_t temp_magnetization = 0.;

    for ( size_t idx = 0; idx != Ising2d_Parameters::N; ++idx )
    {
       temp_magnetization += spin_array[idx];
       // First term is necessary to account for h field loss
       // when 0.5 multiplies the energy to avoid double counting.
#if RFIM
       temp_energy += -0.5 * field_array[idx] * static_cast<float>(spin_array[idx]) + 0.5 * local_energy(idx, spin_array[idx]);
#else
       temp_energy += -0.5 * Ising2d_Parameters::h * static_cast<float>(spin_array[idx]) + 0.5 * local_energy(idx, spin_array[idx]);
#endif
    }

    current_state.energy = temp_energy;
    current_state.magnetization = temp_magnetization;
}

// Change the state by changing the value of the spin
// at a particular site idx. Return what the new state
// would be if the change is accepted (the state is
// passed by reference to avoid extra copies).
template<typename data_t>
void Ising2d<data_t>::change_state(const size_t idx, State<data_t> & temp_state) const
{
    temp_state.DoF = -spin_array[idx];
    temp_state.magnetization = current_state.magnetization + temp_state.DoF - spin_array[idx];
    temp_state.energy = current_state.energy + local_energy(idx, temp_state.DoF) - local_energy(idx, spin_array[idx]);
}

// Set the current state
template<typename data_t>
void Ising2d<data_t>::set_state(const size_t idx, const State<data_t> & _state)
{
    current_state.energy = _state.energy;
    current_state.magnetization = _state.magnetization;
    spin_array[idx] = _state.DoF;
}

// Update the non-energetic observables
template<typename data_t>
void Ising2d<data_t>::update_observables(const size_t bin, Ising2d_Obs<data_t> * obs_ptr) const
{
    const data_t mag_val = Ising2d_Parameters::divide_N * abs( current_state.magnetization );
/*
    std::cout << "\nmag_val = " << std::scientific << mag_val << ", total magnetization = " << abs(current_state.magnetization) << "\n"; 
    std::cout << std::scientific << "magnetization / mag_val = " << abs(current_state.magnetization) / mag_val << "\n";
*/
    obs_ptr -> update_observable_average(mag_val, Obs::enum_names::mag, bin);
    obs_ptr -> update_observable_average(mag_val * mag_val, Obs::enum_names::mag2, bin);
    obs_ptr -> update_observable_average(mag_val * mag_val * mag_val * mag_val, Obs::enum_names::mag4, bin);
    
#if CORRELATION_LENGTHS
    if ( static_cast<size_t> (obs_ptr -> get_observable(Obs::enum_names::counts_per_bin, bin)) % counts_per_transform == 0 )
    {
        data_t Gq_value = Ising2d_Parameters::divide_N * Ising2d_Parameters::divide_N * ( obs_ptr -> correlator.compute_correlator( spin_array, 0, 0, 1 ) );
        obs_ptr -> update_qmin_correlator( Gq_value, Obs::enum_names::corr_qmin, bin,
                                           obs_ptr -> get_observable(Obs::enum_names::counts_per_bin, bin) / counts_per_transform );
    }
#endif

    obs_ptr -> increment_counts_per_bin(bin);
}

template<typename data_t>
void Ising2d<data_t>::print_lattice() const
{
    printf("\n");
    for ( size_t idx = 0; idx != Ising2d_Parameters::L; ++idx )
    {
        printf("\n%ld\t    ", idx * Ising2d_Parameters::L);
        for ( size_t jdx = 0; jdx != Ising2d_Parameters::L; ++jdx )
        {
           if ( spin_array[idx * Ising2d_Parameters::L + jdx] == 1. )
               printf("+ ");
           else
               printf("- ");
        }
    }
}

#if MPI_ON
#ifndef INDEPENDENT_WALKERS
#include <rewl_parameters.hpp>   // Include this for my MPI datatype macros

enum State_SendRecv
{
    energy = 444, mag, state_dof, all_dof
};

// TODO: This should be specialized because it will break probably
template<class state_t>
void mpi_exchange_state( state_t * const state, const int partner_index, const int comm_id, const MPI_Comm * const local_communicators, MPI_Status * const status )
{
    MPI_Sendrecv_replace( &( state -> energy ), 1, MPI_ENERGY_TYPE, partner_index, State_SendRecv::energy, partner_index, State_SendRecv::energy, local_communicators[ comm_id ], status );
    
    MPI_Sendrecv_replace( &( state -> magnetization ), 1, MPI_OBS_TYPE, partner_index, State_SendRecv::mag, partner_index, State_SendRecv::mag, local_communicators[ comm_id ], status );

    MPI_Sendrecv_replace( &( state -> DoF ), 1, MPI_OBS_TYPE, partner_index, State_SendRecv::state_dof, partner_index, State_SendRecv::state_dof, local_communicators[ comm_id ], status );
}

template<typename obs_t>
void mpi_exchange_DoFs( obs_t * const front_dofs, const size_t num_dof, const int partner_index, const int comm_id, const MPI_Comm * const local_communicators, MPI_Status * const status );

// Specialize this templated functions
template<>
void mpi_exchange_DoFs<OBS_TYPE>( OBS_TYPE * const front_dofs, const size_t num_dof, const int partner_index, const int comm_id, const MPI_Comm * const local_communicators, MPI_Status * const status )
{
    MPI_Sendrecv_replace( front_dofs, static_cast<int>(num_dof), MPI_OBS_TYPE, partner_index, State_SendRecv::all_dof, partner_index, State_SendRecv::all_dof, local_communicators[ comm_id ], status );
}

#endif
#endif

#endif


