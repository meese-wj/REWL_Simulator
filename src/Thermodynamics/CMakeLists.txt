add_library(thermodynamics
            order_parameter_cumulants.hpp
            write_self_averaged_observables.hpp
            thermodynamics.hpp)

set_target_properties(thermodynamics PROPERTIES LINKER_LANGUAGE CXX)
