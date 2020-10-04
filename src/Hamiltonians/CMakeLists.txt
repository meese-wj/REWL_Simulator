# The grid directory must be the first
# subdirectory added so as to add the
# subdirectory as an include directory
# properly for all of the Hamiltonians.
add_subdirectory(grid)
add_subdirectory(Ising2d)
add_subdirectory(Ashkin_Teller2d)

add_library(hamiltonians
            model_hamiltonians.hpp)

target_link_libraries(hamiltonians PUBLIC
                      grid_functions
                      ising2d
                      ashkin_teller2d)
