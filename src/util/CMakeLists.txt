add_library(utilities
            array_shift.hpp
            numerical_derivatives.hpp
            vector_matching.hpp
            random_number_generators.hpp
            write_microcanonical_observables.hpp
            mpi_rewl_helpers.hpp
            mpi_rewl_definitions.hpp
            concatenate_rewl_tables.hpp
            command_line_arguments.cpp
            command_line_arguments.hpp)

set_target_properties(utilities PROPERTIES LINKER_LANGUAGE CXX)
