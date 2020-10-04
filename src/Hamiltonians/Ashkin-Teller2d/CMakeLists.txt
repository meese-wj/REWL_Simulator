add_library(ashkin-teller2d
            ashkin-teller2d_parameters.cxx
            ashkin-teller2d_parameter_string.hpp
            ashkin-teller2d_observables.hpp
            ashkin-teller2d_hamiltonian.hpp)

add_subdirectory(source)

target_link_libraries(ashkin-teller2d ashkin-teller2d_source)
