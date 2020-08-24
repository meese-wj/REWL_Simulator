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

#endif
