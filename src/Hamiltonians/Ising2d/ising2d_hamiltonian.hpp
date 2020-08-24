#ifndef ISING2D

// Include the parameters. These will be
// edited by the user before each run.
#include "ising2d_parameters.cxx"

// Include the observables enum class
#include "ising2d_observables.hpp"

// Include the grid setup
#include "grid_setup.hpp"

// Create a namespace to house the 
// parameters.
namespace Ising2d_Parameters
{
    // Hamiltonian parameters
    constexpr size_t L = SYSTEM_SIZE_L;
    constexpr size_t N = L * L;
    constexpr float J = SYSTEM_EXCHANGE_J;
    constexpr float h = SYSTEM_FIELD_H;
    constexpr size_t num_neighbors_i = SYSTEM_NUM_NEAREST_NEIGHBORS_int;
    constexpr float num_neighbors_f = SYSTEM_NUM_NEAREST_NEIGHBORS_float; 

    // Wang Landau energy parameters
    constexpr float energy_min = SYSTEM_ENERGY_MIN;
    constexpr float energy_max = SYSTEM_ENERGY_MAX;
    constexpr float energy_bin_size = SYSTEM_ENERGY_BIN_WIDTH;
    constexpr size_t num_bins = SYSTEM_NUM_BINS;
}

template<typename data_t>
struct State
{
    float energy = 0;
    data_t magnetization = 0;
};

template<typename data_t>
struct Ising2d
{
    State<data_t> current_state;

    float * spin_array = nullptr;
    size_t * neighbor_array = nullptr;

    Ising2d()
    {
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
    }

    ~Ising2d()
    {
       if ( spin_array != nullptr ) delete [] spin_array;
       if ( neighbor_array != nullptr ) delete [] neighbor_array;
    }
};


#endif


