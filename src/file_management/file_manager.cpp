#include "file_manager.hpp"

namespace FS = std::filesystem;

/* ************************************************************************************************************ */
/* Define some constant strings                                                                                 */
const std::string OUTPUT_DATA_PATH = "Simulation_Data";          // Location in build parent to house all data
const std::string HISTOGRAM_SUBFOLDER = "Histograms";            // Location within dated parent subfolder 
/* ************************************************************************************************************ */

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
#if JOB_ARRAYS
FS::path create_output_path( const std::string & model_name, const std::string & todays_date, const std::string & size_string, const std::string & job_id_string )
#else
FS::path create_output_path( const std::string & model_name, const std::string & todays_date, const std::string & size_string )
#endif
{
    FS::path output_path = FS::current_path();
    
    // Put data in outsid of the build directory
    // Path = build parent / OUTPUT_DATA_PATH
    output_path = output_path.parent_path() / OUTPUT_DATA_PATH;

    // Determine if the path exists and create it
    // if it does not
    create_directory( output_path.string() );

    // Add which simulation type it is
    // Path = build parent / OUTPUT_DATA_PATH / model_name
    output_path /= model_name;
    create_directory( output_path.string() );

    // Create a new directory with today's date
    // Path = build parent / OUTPUT_DATA_PATH / model_name / today's date 
    output_path /= todays_date; 
    create_directory( output_path.string() );

#if JOB_ARRAYS
    // Add a subdirectory for the Job Arrays
    output_path /= std::string("Job_Arrays");
    create_directory( output_path.string() );
#endif

    // Finally create a folder for the histograms
    // for a given system size.
    // Path = build parent / OUTPUT_DATA_PATH / model_name / today's date 
    // Subpath = build parent / OUTPUT_DATA_PATH / model_name / today's date / HISTOGRAM_SUBFOLDER
    FS::path sub_path = output_path / HISTOGRAM_SUBFOLDER;
    create_directory( sub_path.string() );
#if JOB_ARRAYS
    sub_path /= std::string("ID_") + job_id_string;
    create_directory( sub_path.string() );
#endif

    // Path = build parent / OUTPUT_DATA_PATH / model_name / today's date 
    // Subpath = build parent / OUTPUT_DATA_PATH / model_name / today's date / HISTOGRAM_SUBFOLDER / size_string_subfolder
    sub_path /= size_string;
    create_directory( sub_path.string() );

#if AT_DENSITIES 
    // Add a folder to house the density plots.
    // Path = build parent / OUTPUT_DATA_PATH / model_name / today's date 
    // density_path = Path / Density_Plots
    FS::path density_path = output_path / std::string("Density_Plots");
    create_directory( density_path.string() );
#endif
 
    return output_path;
}
