cmake_minimum_required(VERSION 3.16)
project(REWL VERSION 0.1.0)
set(CMAKE_BUILD_TYPE Release)

string(TIMESTAMP cmake_time "%Y-%m-%dT%H:%M:%S")
message(STATUS "CMAKE Configure Time: ${cmake_time}")

# Get Git information
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "Git Branch: ${GIT_BRANCH}")
message(STATUS "Git Commit: ${GIT_COMMIT_HASH}")


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
message(STATUS "CMAKE_BINARY_DIR = ${CMAKE_BINARY_DIR}")

add_compile_options(-std=c++17 -Wall -Werror -pedantic)

# *********************************************************************
# Add compile-time options below for the Global Simulation

option(MPI_ON "Turn on MPI. This engages parallelized a Wang Landau simulation." ON)
option(COLLECT_TIMINGS "Print timings collected by STL Chrono." ON)
option(PRINT_HISTOGRAM "Print the energy histogram to a file periodically." OFF)
option(REDUCE_LOGDOS "Truncate the logDoS after each iteration to heighten numerical sensitivity." ON)
option(SAMPLE_AFTER "Sample observables after one round of WL updates" OFF)
option(TRAPEZOIDAL_RULE "Use trapezoidal rule to evaluate the thermodynamical integrals" ON)
option(INDEPENDENT_WALKERS "Do not allow for exchange or for intra-window walkers to average" OFF)
option(DIFFERENT_SEEDS "Use timings to set the seeds for the walkers." ON)
option(EQUAL_WINDOWS "Construct equally-sized windows in REWL." ON)
option(JOB_ARRAYS "Allow for job arrays to run without data over-write." OFF)
option(CORRELATION_LENGTHS "Compute the correlation length periodically in each energy bin." ON)
option(SIMULATED_ANNEALING "Perform a simulated annealing run to find the ground state configuration at the start." OFF)
option(ONE_OVER_T_ALGORITHM "Use the 1/t algorithm as the basis for WL updates instead of exponentially decaying updates." ON)
option(NORMALIZE_BY_STATES "Normalize the final logDoS by the total number of states rather than by fixing the ground state degeneracy." ON)

# Add compile-time options to specify the Hamiltonian
option(ISING2D "Run a 2d Ising model on a periodic square lattice." OFF)
option(ASHKIN_TELLER2D "Run a 2d Askin-Teller model on a periodic square lattice." OFF)

# Additions to Hamiltonians
option(RFIM "Include random fields for the Ising model." OFF)
option(RFAT_BAXTER "Include random fields for the Baxter spin in the Ashkin-Teller model." OFF)
option(PMND_INTERACTIONS "Include phonon-mediated dipolar interactions for Ising-nematic degrees of freedom." OFF)
option(GAUSSIAN_DISORDER "Use a Gaussian random field distribution centered about zero. The default is a uniform distribution." ON)
set(UNIFORM_DISORDER OFF) # This will be turned on if necessary later on.
option(AT_DENSITIES "Keep track of sigma and tau densities per energy bin." OFF)

# *********************************************************************
# Add consequences of the options below.

if (MPI_ON)
    add_definitions(-DMPI_ON)
    find_package(MPI REQUIRED)

    message(STATUS "MPI Flags:")
    list(APPEND CMAKE_MESSAGE_INDENT "  ")
    message(STATUS "MPIEXEC = ${MPIEXEC}")
    message(STATUS "MPIEXEC_NUMPROC_FLAG = ${MPIEXEC_NUMPROC_FLAG}")
    message(STATUS "MPIEXEC_MAX_NUMPROCS = ${MPIEXEC_MAX_NUMPROCS}")
    message(STATUS "MPIEXEC_PREFLAGS = ${MPIEXEC_PREFLAGS}")
    message(STATUS "MPIEXEC_POSTFLAGS = ${MPIEXEC_POSTFLAGS}")
    list(POP_BACK CMAKE_MESSAGE_INDENT)
endif()

# *********************************************************************
# Error out if more than one Hamiltonian is activated.
# TODO: This is probably unsustainable for many models...
if (ISING2D AND ASHKIN_TELLER2D)
    message(FATAL_ERROR "More than one Hamiltonian is activated.")
endif()

# Error out if RFIM and Ashkin-Teller model
if (RFIM AND ASHKIN_TELLER2D)
    message(FATAL_ERROR "RFIM was activated for the Ashkin-Teller model.")
endif()

# Error out if RFAT_BAXTER and Ising model
if (RFAT_BAXTER AND ISING2D)
    message(FATAL_ERROR "RFAT_BAXTER was activated for the Ising model.")
endif()

# Error out if AT_DENSITIES and Ising model
if (AT_DENSITIES AND ISING2D)
    message(FATAL_ERROR "AT_DENSITIES was activated for the Ising model.")
endif()

# Error out if AT_DENSITIES and PMND_INTERACTIONS
if (AT_DENSITIES AND PMND_INTERACTIONS)
    message(FATAL_ERROR "PMND_INTERACTIONS are not implemented for the Ashkin-Teller model.")
endif()
# *********************************************************************

if (ISING2D)
    add_definitions(-DISING2D)
endif()

if (RFIM)
    message(STATUS "RFIM engaged.")
    add_definitions(-DRFIM)
    add_definitions(-DRANDOM_DISORDER)
endif()

if (ASHKIN_TELLER2D)
    add_definitions(-DASHKIN_TELLER2D)
endif()

if (AT_DENSITIES)
    add_definitions(-DAT_DENSITIES)
endif()

if (RFAT_BAXTER)
    add_definitions(-DRFAT_BAXTER)
    add_definitions(-DRANDOM_DISORDER)
endif()

if (PMND_INTERACTIONS)
    add_definitions(-DPHONON_MEDIATED_NEMATIC_INTERACTIONS)
endif()

# Set up the disorder compile-time options
if (RFIM OR RFAT_BAXTER)
    if (GAUSSIAN_DISORDER)
        add_definitions(-DGAUSSIAN_DISORDER)
        set(UNIFORM_DISORDER OFF)
    else()
        add_definitions(-DUNIFORM_DISORDER)
        set(UNIFORM_DISORDER ON)
        set(GAUSSIAN_DISORDER OFF)
    endif()
else()
    set(UNIFORM_DISORDER OFF)
    set(GAUSSIAN_DISORDER OFF)
endif()

if (COLLECT_TIMINGS)
    add_definitions(-DCOLLECT_TIMINGS)
endif()

if (PRINT_HISTOGRAM)
    add_definitions(-DPRINT_HISTOGRAM)
endif()

if (REDUCE_LOGDOS)
    add_definitions(-DREDUCE_LOGDOS)
endif()

if (SAMPLE_AFTER)
    add_definitions(-DSAMPLE_AFTER)
endif()

if (TRAPEZOIDAL_RULE)
    add_definitions(-DTRAPEZOIDAL_RULE)
endif()

if (DIFFERENT_SEEDS)
    add_definitions(-DDIFFERENT_SEEDS)
endif()

if (EQUAL_WINDOWS)
    add_definitions(-DEQUAL_WINDOWS)
endif()

if (JOB_ARRAYS)
    add_definitions(-DJOB_ARRAYS)
endif()

if (CORRELATION_LENGTHS)
    add_definitions(-DCORRELATION_LENGTHS)
endif()

if (SIMULATED_ANNEALING)
    add_definitions(-DSIMULATED_ANNEALING)
endif()

if (ONE_OVER_T_ALGORITHM)
    add_definitions(-DONE_OVER_T_ALGORITHM)
endif()

if (NORMALIZE_BY_STATES)
    add_definitions(-DNORMALIZE_BY_STATES)
endif()

# Now continue with the build
add_subdirectory(src)

if (MPI_ON AND MPI_CXX_FOUND)
    target_include_directories(REWL_Simulator PRIVATE ${MPI_CXX_INCLUDE_PATH})
    target_compile_options(REWL_Simulator PRIVATE ${MPI_CXX_COMPILE_FLAGS})
    target_link_libraries(REWL_Simulator PUBLIC ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
    message(STATUS "")
    message(STATUS "Run code with")
    list(APPEND CMAKE_MESSAGE_INDENT "  ")
    message(STATUS "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ./bin/REWL_Simulator ${MPIEXEC_POSTFLAGS}")
    message(STATUS "")
    list(POP_BACK CMAKE_MESSAGE_INDENT)
endif()

