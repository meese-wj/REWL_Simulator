// Text file enumerating the Ising-2d 
// parameters. For now, this Hamiltonian
// is on a periodic grid at zero field.
// TODO: Add ghost cells for fixed 
// boundary conditions.
// TODO: Add nonzero field.

#ifndef ISING2D_PARAMETER_DEFINITIONS
#define ISING2D_PARAMETER_DEFINITIONS

#include <cmath>

// Linear size of the lattice
constexpr size_t SYSTEM_SIZE_L = 64;

// Exchange constant. J > 0 is ferromagnetic
constexpr float SYSTEM_EXCHANGE_J = 1.0;

// External field strength
constexpr float SYSTEM_FIELD_H = 3.0;

// Number of nearest neighbors
constexpr size_t SYSTEM_NUM_NEAREST_NEIGHBORS_int = 4;
constexpr float SYSTEM_NUM_NEAREST_NEIGHBORS_float = 4;

// Ising Hamiltonian Ground State Degeneracy
// For h == 0, the ground state is doubly-degenerate.
constexpr double GROUND_STATE_DEGENERACY = 2.;

// Wang Landau Energy Minimum
constexpr float SYSTEM_ENERGY_MIN = -0.5 * SYSTEM_EXCHANGE_J * SYSTEM_NUM_NEAREST_NEIGHBORS_float * SYSTEM_SIZE_L * SYSTEM_SIZE_L;

// Wang Landau Energy Maximum
constexpr float SYSTEM_ENERGY_MAX = 0.;

// Wang Landau Bin Width
constexpr float SYSTEM_ENERGY_BIN_WIDTH = SYSTEM_EXCHANGE_J * SYSTEM_NUM_NEAREST_NEIGHBORS_float * 2;

// Wang Landau Bin Number
constexpr size_t SYSTEM_NUM_BINS = static_cast<size_t>( (SYSTEM_ENERGY_MAX - SYSTEM_ENERGY_MIN) / SYSTEM_ENERGY_BIN_WIDTH );

// Create a namespace to house the 
// parameters.
namespace Ising2d_Parameters
{
    // Hamiltonian parameters
    constexpr size_t L = SYSTEM_SIZE_L;
    constexpr size_t N = L * L;
    constexpr size_t num_DoF = N;        // Number of independent degrees of freedom
    constexpr float J = SYSTEM_EXCHANGE_J;
    constexpr float h = SYSTEM_FIELD_H;
    constexpr size_t num_neighbors_i = SYSTEM_NUM_NEAREST_NEIGHBORS_int;
    constexpr float num_neighbors_f = SYSTEM_NUM_NEAREST_NEIGHBORS_float; 
    constexpr double ground_state_degeneracy = GROUND_STATE_DEGENERACY;

    // Wang Landau energy parameters
    constexpr float energy_min = SYSTEM_ENERGY_MIN;
    constexpr float energy_max = SYSTEM_ENERGY_MAX;
    constexpr float energy_bin_size = SYSTEM_ENERGY_BIN_WIDTH;
    constexpr size_t num_bins = SYSTEM_NUM_BINS;
}

#endif
