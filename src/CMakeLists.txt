
set(INCLUDE_DIRS "${INCLUDE_DIRS}"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_management"
    "${CMAKE_CURRENT_SOURCE_DIR}/Hamiltonians/grid"
    "${CMAKE_CURRENT_SOURCE_DIR}/Hamiltonians"
    "${CMAKE_CURRENT_SOURCE_DIR}/util"
    "${CMAKE_CURRENT_SOURCE_DIR}/Thermodynamics"
    "${CMAKE_CURRENT_SOURCE_DIR}/REWL")

message(STATUS "Include directories: ${INCLUDE_DIRS}")

include_directories("${INCLUDE_DIRS}")

add_subdirectory(file_management)
add_subdirectory(util)
add_subdirectory(Thermodynamics)
add_subdirectory(Hamiltonians)
add_subdirectory(REWL)
add_subdirectory(Simulation)

target_link_libraries(REWL_Simulator PUBLIC
                      files
                      utilities
                      thermodynamics
                      hamiltonians
                      rewl_functions
                      simulation)

if (MPI_ON AND MPI_CXX_FOUND)
    target_include_directories(REWL_Simulator PRIVATE ${MPI_CXX_INCLUDE_PATH})
    target_compile_options(REWL_Simulator PRIVATE ${MPI_CXX_COMPILE_FLAGS})
    target_link_libraries(REWL_Simulator PUBLIC ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
endif()

