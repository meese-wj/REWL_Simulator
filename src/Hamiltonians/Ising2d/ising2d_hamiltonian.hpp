#ifndef ISING2D_HAMILTONIAN

// Include the parameters through the 
// string header to avoid redefinitions of
// the constexpr values.
#include "ising2d_parameter_string.hpp"

// Include the observables enum class
#include "ising2d_observables.hpp"

// Include the grid setup from the 
// cmake include directories.
#include <grid_setup.hpp>

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
    size_t * neighbor_array = nullptr;

    // Add some Hamiltonian dependent functions
    float local_field(const size_t idx) const;
    float local_energy(const size_t idx, const data_t spin_value) const;
    void recalculate_state();
    void change_state(const size_t idx, State<data_t> & temp_state) const;
    void set_state(const size_t idx, const State<data_t> & _state);
    void update_observables(const size_t bin, Ising2d_Obs<data_t> * obs_ptr) const;

    // Finally add the constructor and destructor.
    Ising2d()
    {
        spin_array = new data_t [ Ising2d_Parameters::N ];
        neighbor_array = new size_t [ Ising2d_Parameters::N * Ising2d_Parameters::num_neighbors_i ];
        // Initialize all spins to 1.
        // TODO: change this to be randomized?
        for ( size_t idx = 0; idx != Ising2d_Parameters::N; ++idx )
            spin_array[idx] = 1.;

        // TODO: Generalize this to different types 
        // of grids.
        define_2d_square_periodic_neighbors(Ising2d_Parameters::L,
                                            Ising2d_Parameters::L,
                                            Ising2d_Parameters::num_neighbors_i,
                                            neighbor_array);

        recalculate_state();
    }

    ~Ising2d()
    {
       if ( spin_array != nullptr ) delete [] spin_array;
       if ( neighbor_array != nullptr ) delete [] neighbor_array;
    } 

    void print_lattice() const;
};

template<typename data_t>
float Ising2d<data_t>::local_field(const size_t idx) const
{
    float field = Ising2d_Parameters::h;
    for ( size_t nidx = 0; nidx != Ising2d_Parameters::num_neighbors_i; ++nidx )
    {
        field += Ising2d_Parameters::J * static_cast<float>( spin_array[ neighbor_array[ idx * Ising2d_Parameters::num_neighbors_i + nidx ] ] );
    }

    return field;
}

template<typename data_t>
float Ising2d<data_t>::local_energy(const size_t idx, const data_t spin_value) const
{
    return -1. * static_cast<float>( spin_value ) * local_field(idx);
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
       temp_energy += -0.5 * Ising2d_Parameters::h * static_cast<float>(spin_array[idx]) + 0.5 * local_energy(idx, spin_array[idx]);
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
    const data_t mag_val = abs( current_state.magnetization );
    obs_ptr -> update_observable_average(mag_val, Obs::enum_names::mag, bin);
    obs_ptr -> update_observable_average(mag_val * mag_val, Obs::enum_names::mag2, bin);
    obs_ptr -> update_observable_average(mag_val * mag_val * mag_val * mag_val, Obs::enum_names::mag4, bin);
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

#endif


