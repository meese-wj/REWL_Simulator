add_library(simulated_annealer
            simulated_annealer.hpp
            simulated_annealing_utilities.hpp
            simulated_annealing_parameter_strings.hpp
            simulated_annealing_parameters.cxx)

set_target_properties(simulated_annealer PROPERTIES LINKER_LANGUAGE CXX)
