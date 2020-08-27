add_subdirectory(Wang_Landau)

add_library(rewl_functions
            glazier.hpp
            histogram_index.hpp
            rewl_histograms.hpp
            REWL.cpp)

target_link_libraries(rewl_functions wang_landau_functions)
