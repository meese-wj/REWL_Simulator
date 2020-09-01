add_library(simulation
            rewl_simulation.hpp
            rewl_simulation.cpp)

set_target_properties(simulation PROPERTIES LINKER_LANGUAGE CXX)

add_executable(REWL_Simulator main.cpp)
