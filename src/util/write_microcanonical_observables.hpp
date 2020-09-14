#ifndef WRITE_MICROCANONICAL_OBSERVABLES
/* This file will write out observables as a
 * function of the energy. These are not 
 * thermally averaged. */

#include <string>
#include <fstream>

// This string must be defined somewhere!
// const std::string DELIMITER = "  ";

template<typename energy_t, typename logdos_t, typename obs_t>
void write_microcanonical_observables( const size_t num_bins, const size_t num_obs,
                                       const std::string & file_path,
                                       const std::string & file_header,
                                       const energy_t * const energy_array,
                                       const logdos_t * const logdos_array,
                                       const obs_t * const system_observables )
{
    const std::string file_name = "microcanonical_observables.txt"; 
    std::ofstream output_file;
    output_file.open( file_path + "/" + file_name );

    output_file << file_header;
    for ( size_t bin = 0; bin != num_bins; ++bin )
    {
        output_file << energy_array[ bin ] << DELIMITER << logdos_array[ bin ] << DELIMITER;

        for ( size_t ob = 0; ob != num_bins; ++ob )
        {
            output_file << system_observables[ bin * num_obs + ob ];
            if ( ob != num_obs - 1 ) output_file << DELIMITER;
        }

        if ( bin != num_bins - 1 ) output_file << "\n";
    }
}

#endif
