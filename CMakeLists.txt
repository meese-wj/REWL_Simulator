cmake_minimum_required(VERSION 3.16)
project(REWL VERSION 0.1.0)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_compile_options(-Wall -Wextra -Werror -pedantic)

option(ISING2D "Run a 2d Ising model on a periodic square lattice." ON)

if (ISING2D)
    add_definitions(-DISING2D)
endif()

add_subdirectory(src)
