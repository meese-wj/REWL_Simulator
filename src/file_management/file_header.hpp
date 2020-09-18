#ifndef FILE_HEADER
/* Construct the file header string. */
#include <string>
#include <vector>

std::string create_file_header( const std::string & system_header,
                                const std::string & rewl_header );

std::vector<std::string> concatenate_vector_string ( const std::vector<std::string> & vec1,
                                                     const std::vector<std::string> & vec2 );

#endif
