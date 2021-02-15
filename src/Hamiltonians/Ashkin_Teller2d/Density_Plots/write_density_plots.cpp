#ifndef AT_DENSITY_PLOT_WRITER
#define AT_DENSITY_PLOT_WRITER
/* This file will control the export of
 * the microcanonical Ashkin-Teller 
 * density plots per energy bin. */

#include <string>
#include <fstream>
#include <filesystem>
#include "ashkin_teller_densities_parameters.cxx"
#include "ashkin_teller_densities.hpp"

void write_density_plots( const density_int nbins,
                          const std::filesystem::path & density_path,
                          const density_float * const float_data,
                          const std::string & file_string,
                          const std::string & density_header )
{
    for ( density_int energy_bin = 0; energy_bin != nbins; ++energy_bin )
    {
        const std::string file_name = "density_plots_bin-" + std::to_string( energy_bin ) + file_string + ".txt";
        std::ofstream output_file;
        output_file.open( density_path / file_name );

        output_file << density_header;
        output_file << "\n";
        
        for ( density_int idx = 0; idx != AT_Density_Parameters::total_bins; ++idx )
        {
            output_file << std::setprecision(std::numeric_limits<density_float>::digits10) << std::scientific;
            output_file << float_data[ energy_bin_density_pointer( energy_bin ) + idx ];
            
            if ( (idx + 1) % AT_Density_Parameters::axis_bins == 0 )
            {
                if ( idx != AT_Density_Parameters::total_bins - 1 )
                    output_file << "\n";
            }
            else
                output_file << "  ";
        }

        output_file.close();
    }
}

#endif
