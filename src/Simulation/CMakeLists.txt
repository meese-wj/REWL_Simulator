add_library(simulation
            rewl_simulation.hpp)

set_target_properties(simulation PROPERTIES LINKER_LANGUAGE CXX)

add_executable(REWL_Simulator main.cpp)
