// Use this file to read in various
// Hamiltonians at compile time.

#if ISING2D
    #include "Ising2d/ising2d_hamiltonian.hpp"
    template<typename data_t>
    using Observables_t = Ising2d_Obs<data_t>;
    
    template<typename data_t>
    using Hamiltonian_t = Ising2d<data_t>;
#endif
