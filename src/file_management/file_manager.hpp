#ifndef FILE_MANAGER
/* Set up a file to take care of the file management
 * through C++. The idea is to streamline data output
 * and minimize data-overwrite. */
#include <string>
#include <sstream>
#include <chrono>
#include <iostream>
#include <filesystem>   // Requires C++17

/* ********************************************************************************************************** */
/* Define some constant strings                                                                               */
const std::string data_path = "Simulation_Data";          // Location in build parent to house all data
const std::string histogram_subfolder = "Histograms";     // Location within dated parent subfolder 
const std::string DoF_subfolder = "DoF-";                 // Location for size-dependent files
/*                                                                                                            */
/* ********************************************************************************************************** */

// Get the date string from chrono
std::string get_todays_date();

// Create a directory along the relative path 
// unless that directory already exists.
void create_directory( const std::string & path_to_dir );

// Build the output file path and return it
std::filesystem::path create_output_path( const std::string & model_name );

#endif
