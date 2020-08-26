// Use this file to read in various
// Hamiltonians at compile time.

#if ISING2D
    #include "Ising2d/ising2d_hamiltonian.hpp"
    template<typename data_t>
    using Hamiltonian = Ising2d<data_t>;
#endif
