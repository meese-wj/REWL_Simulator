add_library(utilities
            array_shift.hpp
            random_number_generators.hpp
            mpi_rewl_helpers.hpp)

set_target_properties(utilities PROPERTIES LINKER_LANGUAGE CXX)
