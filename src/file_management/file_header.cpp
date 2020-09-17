#include "file_header.hpp"

std::string create_file_header( const std::string & system_header,
                                const std::string & rewl_header )
{
    std::string output_header = system_header + "\n#\n" + rewl_header + "\n#\n";
    return output_header;
}
