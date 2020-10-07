#ifndef WRITE_NONLINEAR_OBSERVABLES
#define WRITE_NONLINEAR_OBSERVABLES
/* This file will write the nonlinear observables
 * to a specified file.  */

#include <vector>
#include <file_manager.hpp>
#include <fstream>

// This string needs to be defined somewhere!
// const std::string DELIMITER = "  ";

template<typename energy_t, typename obs_t>
void write_nonlinear_obs_to_file( const size_t num_temps, const size_t num_obs,
                                  const std::string & file_string,
                                  const std::string & file_header,
                                  const std::vector<std::string> nonlinear_obs_names,
                                  const std::filesystem::path & file_path,
                                  const energy_t * const temperature_array,
                                  const obs_t * const nonlinear_obs_array )
{
    std::string filename = "nonlinear_observables-" + file_string + ".txt";
    std::ofstream output_file; 

    output_file.open(file_path / filename);
    output_file << file_header;
    
    size_t counter = 1;
    output_file << "# Extensive Nonlinear Observable Names by Column";
    output_file << "\n#    1: Temperature";
    for ( size_t obs_idx = 0, num_obs = nonlinear_obs_names.size(); obs_idx != num_obs; ++obs_idx )
    {
        if ( nonlinear_obs_names[ obs_idx ].compare("NUM OBS") != 0 )
            output_file << "\n#    " << ++counter << ": " << nonlinear_obs_names[ obs_idx ];
    }
    output_file << "\n#\n";

    for ( size_t Tidx = 0; Tidx != num_temps; ++Tidx )
    {
        output_file << temperature_array[ Tidx ] << DELIMITER;
        
        for ( size_t ob = 0; ob != num_obs; ++ob )
        {
            output_file << nonlinear_obs_array[ Tidx * num_obs + ob ];

            if ( ob != num_obs - 1 ) output_file << DELIMITER; 
        }

        if ( Tidx != num_temps - 1 )
            output_file << "\n";
    }


    output_file.close();
}

#endif
