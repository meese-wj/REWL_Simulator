
set(input_config "build.config")
set(output_config "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/build_config.txt")
configure_file(${input_config} ${output_config})


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

