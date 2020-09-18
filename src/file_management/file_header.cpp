#include "file_header.hpp"

std::string create_file_header( const std::string & system_header,
                                const std::string & rewl_header )
{
    std::string output_header = system_header + "\n#\n" + rewl_header + "\n#\n";
    return output_header;
}


std::vector<std::string> concatenate_vector_string ( const std::vector<std::string> & vec1,
                                                     const std::vector<std::string> & vec2 )
{
    std::vector<std::string> result (vec1);
    result.insert( result.end(), vec2.begin(), vec2.end() );
    return result;
}
