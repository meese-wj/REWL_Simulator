add_library(utilities
            array_shift.hpp
            numerical_derivatives.hpp
            random_number_generators.hpp
            write_microcanonical_observables.hpp
            mpi_rewl_helpers.hpp
            mpi_rewl_definitions.hpp)

set_target_properties(utilities PROPERTIES LINKER_LANGUAGE CXX)
