add_library(ising2d
            ising2d_parameters.cxx
            ising2d_parameter_string.hpp
            ising2d_observables.hpp
            ising2d_hamiltonian.hpp
            ising_site.hpp)

add_subdirectory(source)
add_subdirectory(PM_Dipolar_Interactions)

target_link_libraries(ising2d ising2d_source pmd_interactions)
