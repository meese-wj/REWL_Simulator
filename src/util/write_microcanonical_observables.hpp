#ifndef WRITE_MICROCANONICAL_OBSERVABLES
#define WRITE_MICROCANONICAL_OBSERVABLES
/* This file will write out observables as a
 * function of the energy. These are not 
 * thermally averaged. These observables are
 * intensive, meaning they are normalized by
 * system size (once). */

#include <vector>
#include <file_manager.hpp>
#include <fstream>

// This string must be defined somewhere!
// const std::string DELIMITER = "  ";

template<typename energy_t, typename logdos_t, typename obs_t>
void write_microcanonical_observables( const size_t system_size, 
                                       const size_t num_bins, const size_t num_obs,
                                       const size_t counts_per_bin,
                                       const std::string & file_string,
                                       const std::string & file_header,
                                       const std::vector<std::string> & obs_names,
                                       const std::filesystem::path & file_path,
                                       const energy_t * const energy_array,
                                       const logdos_t * const logdos_array,
                                       const obs_t * const system_observables )
{
    const std::string file_name = "microcanonical_observables-" + file_string + ".txt"; 
    std::ofstream output_file;
    output_file.open( file_path / file_name );

    output_file << file_header;
    size_t counter = 2;
    output_file << "# Intensive Observable Names by Column";
    output_file << "\n#    1: Energy\n#    2: logDoS";
    for ( size_t obs_idx = 0, num_obs = obs_names.size(); obs_idx != num_obs; ++obs_idx )
    {
        if ( obs_names[ obs_idx ].compare("NUM OBS") != 0 )
            output_file << "\n#    " << ++counter << ": " << obs_names[ obs_idx ];
    }
    output_file << "\n#\n";

    for ( size_t bin = 0; bin != num_bins; ++bin )
    {
        output_file << std::setprecision(std::numeric_limits<energy_t>::digits10) << std::scientific;
        output_file << energy_array[ bin ] / system_size << DELIMITER;
        output_file << std::setprecision(std::numeric_limits<logdos_t>::digits10) << std::scientific;
        output_file << logdos_array[ bin ] / system_size << DELIMITER;
        output_file << std::setprecision(std::numeric_limits<obs_t>::digits10) << std::scientific;

        for ( size_t ob = 0; ob != num_obs; ++ob )
        {
            if ( ob == counts_per_bin )
                output_file << system_observables[ bin * num_obs + ob ];
            else
                output_file << system_observables[ bin * num_obs + ob ] / system_size; 

            if ( ob != num_obs - 1 ) output_file << DELIMITER;
        }

        if ( bin != num_bins - 1 ) output_file << "\n";
    }
}

#endif
