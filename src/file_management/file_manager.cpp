#include "file_manager.hpp"

namespace FS = std::filesystem;

// Get the date string from chrono
std::string get_todays_date()
{
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    std::time_t time_now = std::chrono::system_clock::to_time_t( now );
    std::stringstream ss_date;
    ss_date << std::put_time( std::localtime( &time_now ), "%B_%d_%Y" );    // Full Month _ Day _ Year
    std::string date = ss_date.str();
    return date;
}

// Create a directory along the relative path 
// unless that directory already exists.
void create_directory( const std::string & path_to_dir )
{
    FS::path filepath = path_to_dir; 
    if ( FS::exists( filepath ) )
    {
        std::cout << "\nInactive filesystem:\n\tDirectory " << filepath << " already exists.\n\n";
        return;
    }
    
    std::cout << "\nCreating directory...\n\tDirectory " << filepath << "\n";
    FS::create_directory( filepath );
}

// Build the output file path and return it
FS::path create_output_path( const std::string & model_name, const std::string & size_string )
{
    FS::path output_path = FS::current_path();
    
    // Put data in outsid of the build directory
    // Path = build parent / data_path
    output_path = output_path.parent_path() / data_path;

    // Determine if the path exists and create it
    // if it does not
    create_directory( output_path.string() );

    // Add which simulation type it is
    // Path = build parent / data_path / model_name
    output_path /= model_name;
    create_directory( output_path.string() );

    // Create a new directory with today's date
    // Path = build parent / data_path / model_name / today's date 
    output_path /= get_todays_date(); 
    create_directory( output_path.string() );

    // Finally create a folder for the histograms
    // for a given system size.
    // Path = build parent / data_path / model_name / today's date 
    // Subpath = build parent / data_path / model_name / today's date / histogram_subfolder
    FS::path sub_path = output_path / histogram_subfolder;
    create_directory( sub_path.string() );

    // Path = build parent / data_path / model_name / today's date 
    // Subpath = build parent / data_path / model_name / today's date / histogram_subfolder / size_string_subfolder
    sub_path /= size_string;
    create_directory( sub_path.string() );
 
    return output_path;
}
