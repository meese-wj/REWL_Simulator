add_library(ashkin_teller2d
            ashkin_teller2d_parameters.cxx
            ashkin_teller2d_parameter_string.hpp
            ashkin_teller2d_observables.hpp
            ashkin_teller2d_hamiltonian.hpp)

add_subdirectory(source)
add_subdirectory(Density_Plots)

target_link_libraries(ashkin_teller2d ashkin_teller2d_source AT_densities)
