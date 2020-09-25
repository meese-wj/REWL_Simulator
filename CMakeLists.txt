cmake_minimum_required(VERSION 3.16)
project(REWL VERSION 0.1.0)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_compile_options(-std=c++17 -Wall -Werror -pedantic)

# *********************************************************************
# Add compile-time options below.

option(MPI_ON "Turn on MPI. This engages parallelized a Wang Landau simulation." ON)
option(ISING2D "Run a 2d Ising model on a periodic square lattice." ON)
option(COLLECT_TIMINGS "Print timings collected by STL Chrono." ON)
option(PRINT_HISTOGRAM "Print the energy histogram to a file periodically." ON)

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

if (ISING2D)
    add_definitions(-DISING2D)
endif()

if (COLLECT_TIMINGS)
    add_definitions(-DCOLLECT_TIMINGS)
endif()

if (PRINT_HISTOGRAM)
    add_definitions(-DPRINT_HISTOGRAM)
endif()

add_subdirectory(src)
