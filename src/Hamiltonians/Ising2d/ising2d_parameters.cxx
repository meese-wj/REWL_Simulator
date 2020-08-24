// Text file enumerating the Ising-2d 
// parameters. For now, this Hamiltonian
// is on a periodic grid at zero field.
// TODO: Add ghost cells for fixed 
// boundary conditions.
// TODO: Add nonzero field.

#if ISING2D

#include <cmath>

// Linear size of the lattice
constexpr size_t SYSTEM_SIZE_L = 8;

// Exchange constant. J > 0 is ferromagnetic
constexpr float SYSTEM_EXCHANGE_J = 1.0;

// External field strength
constexpr float SYSTEM_FIELD_H = 0.;

// Number of nearest neighbors
constexpr size_t SYSTEM_NUM_NEAREST_NEIGHBORS_int = 4;
constexpr float SYSTEM_NUM_NEAREST_NEIGHBORS_float = 4;

// Wang Landau Energy Minimum
constexpr float SYSTEM_ENERGY_MIN = -0.5 * SYSTEM_EXCHANGE_J * SYSTEM_NUM_NEAREST_NEIGHBORS_float * SYSTEM_SIZE_L * SYSTEM_SIZE_L;

// Wang Landau Energy Maximum
constexpr float SYSTEM_ENERGY_MAX = 0.;

// Wang Landau Bin Width
constexpr float SYSTEM_ENERGY_BIN_WIDTH = SYSTEM_EXCHANGE_J * SYSTEM_NUM_NEAREST_NEIGHBORS_float; 

// Wang Landau Bin Number
constexpr size_t SYSTEM_NUM_BINS = static_cast<size_t>( (SYSTEM_ENERGY_MAX - SYSTEM_ENERGY_MIN) / SYSTEM_ENERGY_BIN_WIDTH );

#endif
