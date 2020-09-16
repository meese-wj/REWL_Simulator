#include "file_manager.hpp"

namespace FS = std::filesystem;

void create_directory( const std::string & path_to_dir )
{
    FS::path filepath = path_to_dir; 
    if ( FS::exists( filepath ) )
    {
        std::cout << "\nInactice filesystem:\n\tDirectory " << filepath << " already exists.\n\n";
        return;
    }

    FS::create_directory( filepath );
}
