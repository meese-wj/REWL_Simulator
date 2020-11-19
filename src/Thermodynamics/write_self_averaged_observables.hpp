#ifndef WRITE_SELF_AVERAGED_OBSERVABLES
#define WRITE_SELF_AVERAGED_OBSERVABLES
/* This file will write the self-averaged observables
 * and the energy observables to a specified file.
 * The more complicated cumulants must be written 
 * elsewhere since they are model-dependent. */

#include <vector>
#include <file_manager.hpp>
#include <fstream>

// This string must be defined somewhere!
// const std::string DELIMITER = "  ";

template<typename energy_t, typename obs_t>
void write_observables_to_file( const size_t num_temps, const size_t num_obs,
                                const std::string & file_string,
                                const std::string & file_header,
                                const std::vector<std::string> obs_names,
                                const std::filesystem::path & file_path,
                                const energy_t * const temperature_array,
                                const obs_t * const observables_array )
{
    std::string filename = "self_averaged_observables-" + file_string + ".txt";
    std::ofstream output_file; 

    output_file.open(file_path / filename);
    output_file << file_header;
    
    size_t counter = 1;
    output_file << "# Intensive Observable Names by Column";
    output_file << "\n#    1: Temperature";
    for ( size_t obs_idx = 0, num_obs = obs_names.size(); obs_idx != num_obs; ++obs_idx )
    {
        if ( obs_names[ obs_idx ].compare("NUM OBS") != 0 )
            output_file << "\n#    " << ++counter << ": " << obs_names[ obs_idx ];
    }
    output_file << "\n#\n";

    for ( size_t Tidx = 0; Tidx != num_temps; ++Tidx )
    {
        output_file << std::setprecision(std::numeric_limits<energy_t>::digits10) << std::scientific;
        output_file << temperature_array[ Tidx ] << DELIMITER; 
        output_file << std::setprecision(std::numeric_limits<obs_t>::digits10) << std::scientific;
        
        for ( size_t ob = 0; ob != num_obs; ++ob )
        {
            output_file << observables_array[ Tidx * num_obs + ob ];

            if ( ob != num_obs - 1 ) output_file << DELIMITER; 
        }

        if ( Tidx != num_temps - 1 )
            output_file << "\n";
    }


    output_file.close();
}

#endif
