add_library(thermodynamics
            order_parameter_cumulants.hpp
            self_averaged_observables_writer.hpp
            thermodynamics.hpp)

set_target_properties(thermodynamics PROPERTIES LINKER_LANGUAGE CXX)
