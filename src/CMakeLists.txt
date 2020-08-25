
add_subdirectory(Hamiltonians)
add_subdirectory(REWL)

add_executable(REWL_Simulator main.cpp)

target_link_libraries(REWL_Simulator PUBLIC
                      grid_functions
                      ising2d
                      rewl_functions)
