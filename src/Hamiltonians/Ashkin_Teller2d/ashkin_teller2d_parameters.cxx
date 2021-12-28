// Text file enumerating the AshkinTeller2d 
// parameters. For now, this Hamiltonian
// is on a periodic grid at zero field.
// TODO: Add ghost cells for fixed 
// boundary conditions.
// TODO: Add nonzero field.

#ifndef ASKIN_TELLER2D_PARAMETER_DEFINITIONS
#define ASKIN_TELLER2D_PARAMETER_DEFINITIONS

#include <cmath>

// Linear size of the lattice
constexpr size_t SYSTEM_SIZE_L = 24;

// Exchange constant. J > 0 is ferromagnetic
constexpr float SYSTEM_EXCHANGE_J = 1.0;

// Mixed spin coupling K > 0 is ferromagnetic
constexpr float SYSTEM_EXCHANGE_K = 0.5;

#if RFAT_BAXTER
// Random field strength
constexpr float BAXTER_FIELD_H = 0.8;
#endif

// Number of nearest neighbors
constexpr size_t SYSTEM_NUM_NEAREST_NEIGHBORS_int = 4;
constexpr float SYSTEM_NUM_NEAREST_NEIGHBORS_float = 4;

// Askin-Teller Hamiltonian Ground State Degeneracy
constexpr double GROUND_STATE_DEGENERACY = 4.;

// Wang Landau Energy Minimum
constexpr float SYSTEM_ENERGY_MIN = -0.5 * (2 * SYSTEM_EXCHANGE_J + SYSTEM_EXCHANGE_K) * SYSTEM_NUM_NEAREST_NEIGHBORS_float * SYSTEM_SIZE_L * SYSTEM_SIZE_L; 

// Wang Landau Energy Maximum
constexpr float SYSTEM_ENERGY_MAX = 0.;

// Wang Landau Bin Width
constexpr float SYSTEM_ENERGY_BIN_WIDTH = (SYSTEM_EXCHANGE_J + SYSTEM_EXCHANGE_K) * SYSTEM_NUM_NEAREST_NEIGHBORS_float * 2;

// Wang Landau Bin Number
constexpr size_t SYSTEM_NUM_BINS = static_cast<size_t>( (SYSTEM_ENERGY_MAX - SYSTEM_ENERGY_MIN) / SYSTEM_ENERGY_BIN_WIDTH );

// Create a namespace to house the 
// parameters.
namespace Ashkin_Teller2d_Parameters
{
    // Hamiltonian parameters
    constexpr size_t L = SYSTEM_SIZE_L;
    constexpr size_t N = L * L;
    constexpr double divide_N = 1. / static_cast<double>(N);
    constexpr size_t num_DoF = 2 * N;     // Number of independent degrees of freedom.
    constexpr float J = SYSTEM_EXCHANGE_J;
    constexpr float K = SYSTEM_EXCHANGE_K;
#if RFAT_BAXTER
    // Field strength for the Baxter phase
    constexpr float h = BAXTER_FIELD_H;
#endif
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
