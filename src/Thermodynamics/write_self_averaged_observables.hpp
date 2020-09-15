#ifndef SELF_AVERAGED_OBSERVABLES_WRITER
/* This file will write the self-averaged observables
 * and the energy observables to a specified file.
 * The more complicated cumulants must be written 
 * elsewhere since they are model-dependent. */

#include <string>
#include <fstream>

// This string must be defined somewhere!
// const std::string DELIMITER = "  ";

template<typename energy_t, typename obs_t>
void write_observables_to_file(                             const size_t num_temps, const size_t num_obs,
                                const std::string & file_path,
                                const std::string & file_header,
                                const energy_t * const temperature_array,
                                const obs_t * const observables_array )
{
    std::string filename = "self_averaged_observables.txt";
    std::ofstream output_file; 

    output_file.open(file_path + "/" + filename);
    output_file << file_header;

    for ( size_t Tidx = 0; Tidx != num_temps; ++Tidx )
    {
        output_file << temperature_array[ Tidx ] << DELIMITER;
        
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
