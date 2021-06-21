
# Add the build config file
set(input_config "build.config")
set(output_config "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/build_config.txt")
configure_file(${input_config} ${output_config})
message(STATUS "Build config file written to ${output_config}")

# Now copy the relevant parameter files
# into the the same directory
set(rewl_parameters "REWL/rewl_parameters.hpp")
set(rewl_param_output "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/rewl_parameters.txt")
configure_file(${rewl_parameters} ${rewl_param_output} COPYONLY)
message(STATUS "REWL Parameter file copied to ${rewl_param_output}")

set(model_parameters "")
if (ISING2D)
    set(model_parameters "Hamiltonians/Ising2d/ising2d_parameters.cxx")
elseif (ASHKIN_TELLER2D)
    set(model_parameters "Hamiltonians/Ashkin_Teller2d/ashkin_teller2d_parameters.cxx")
endif()
set(model_param_output "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/model_parameters.txt")
configure_file(${model_parameters} ${model_param_output} COPYONLY)
message(STATUS "Model Parameter file copied to ${model_param_output}")

if(SIMULATED_ANNEALING)
    set(sa_parameters "Hamiltonians/Simulated_Annealer/simulated_annealing_parameters.cxx")
    set(sa_param_output "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/simulated_annealing_parameters.txt")
    configure_file(${sa_parameters} ${sa_param_output} COPYONLY)
    message(STATUS "Simulated Annealing Parameter file copied to ${model_param_output}")
endif()


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

