add_library(utilities
            random_number_generators.hpp
            mpi_rewl_helpers.hpp)

set_target_properties(utilities PROPERTIES LINKER_LANGUAGE CXX)
