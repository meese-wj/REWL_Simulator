#ifndef MODEL_HAMILTONIANS
#define MODEL_HAMILTONIANS
// Use this file to read in various
// Hamiltonians at compile time.

#if ISING2D
    #include "Ising2d/ising2d_hamiltonian.hpp"
    namespace System_Parameters = Ising2d_Parameters;

    using System_Strings = Ising2d_Parameter_String;

    template<typename data_t>
    using State_t = State<data_t>;
   
    namespace System_Obs = Obs;
    using System_Obs_enum_t = System_Obs::enum_names;

    template<typename data_t>
    using Observables_t = Ising2d_Obs<data_t>;
    
    template<typename data_t>
    using Hamiltonian_t = Ising2d<data_t>;
#endif
#endif
