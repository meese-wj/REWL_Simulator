#ifndef FILE_MANAGER
/* Set up a file to take care of the file management
 * through C++. The idea is to streamline data output
 * and minimize data-overwrite. */
#include <string>
#include <iostream>
#include <filesystem>   // Requires C++17

// Create a directory along the relative path 
// unless that directory already exists.
void create_directory( const std::string & dirname );

#endif
