add_library(thermodynamics
            order_parameter_cumulants.hpp
            thermodynamics.hpp)

set_target_properties(thermodynamics PROPERTIES LINKER_LANGUAGE CXX)
