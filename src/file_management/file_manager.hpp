#ifndef FILE_MANAGER
#define FILE_MANAGER
/* Set up a file to take care of the file management
 * through C++. The idea is to streamline data output
 * and minimize data-overwrite. */
#include <string>
#include <sstream>
#include <chrono>
#include <iostream>
#include <filesystem>   // Requires C++17

// Get the date string from chrono
std::string get_todays_date();

// Create a directory along the relative path 
// unless that directory already exists.
void create_directory( const std::string & path_to_dir );

// Build the output file path and return it
#if JOB_ARRAYS
std::filesystem::path create_output_path( const std::string & data_parent_directory, 
                                          const std::string & model_name, 
                                          const std::string & todays_date, const std::string & size_string, 
                                          const std::string & job_id_string );
#else
std::filesystem::path create_output_path( const std::string & data_parent_directory, 
                                          const std::string & model_name, 
                                          const std::string & todays_date, const std::string & size_string );
#endif

#endif
