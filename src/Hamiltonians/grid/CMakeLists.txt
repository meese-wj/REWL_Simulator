add_library(grid_functions
            grid_setup.hpp
            grid_setup.cpp)

target_include_directories(grid_functions PUBLIC
                           "${CMAKE_CURRENT_SOURCE_DIR}")
