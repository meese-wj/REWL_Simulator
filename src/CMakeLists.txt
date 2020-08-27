
set(INCLUDE_DIRS
    "${CMAKE_CURRENT_SOURCE_DIR}/Hamiltonians/grid"
    "${CMAKE_CURRENT_SOURCE_DIR}/Hamiltonians"
    "${CMAKE_CURRENT_SOURCE_DIR}/util")

include_directories("${INCLUDE_DIRS}")

add_subdirectory(util)
add_subdirectory(Hamiltonians)
add_subdirectory(REWL)

add_executable(REWL_Simulator main.cpp)

target_link_libraries(REWL_Simulator PUBLIC
                      utilities
                      hamiltonians
                      rewl_functions)
