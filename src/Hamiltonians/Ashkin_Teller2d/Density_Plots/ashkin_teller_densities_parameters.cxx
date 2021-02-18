#ifndef AT_DENSITIES_PARAMETERS
#define AT_DENSITIES_PARAMETERS
/* This file will store the parameters
 * for the 2d density histograms. */

#include <string>
#include <cstdint>
using density_int = std::uint32_t;  // This maxes out at 2^31 - 1 = 4.29e9
                                    // If more counts are needed, use 64 bits
using density_float = double;

namespace AT_Density_Parameters
{
    // Number of bins along the sigma or tau axis
    // For memory reasons it would be best to keep this
    // a power of 2.
    constexpr density_int axis_bins = 32;

    // Total number of bins in the 2d plot
    constexpr density_int total_bins = axis_bins * axis_bins;

    // Define the maximum and minimum values
    // in the 2d plane. This should probably
    // just be +/- 1. 
    constexpr density_float sigma_tau_max = 1.05;
    constexpr density_float sigma_tau_min = -sigma_tau_max;
    constexpr density_float density_binwidth = ( sigma_tau_max - sigma_tau_min ) / static_cast<density_float>( axis_bins );


    /* Store the density parameters in a 
     * set of strings to be used in the header
     * of the output files. */
    struct Strings
    {
        const std::string axis_bins = std::to_string( AT_Density_Parameters::axis_bins );
        const std::string total_bins = std::to_string( AT_Density_Parameters::total_bins );
        const std::string sigma_tau_max = std::to_string( AT_Density_Parameters::sigma_tau_max );
        const std::string sigma_tau_min = std::to_string( AT_Density_Parameters::sigma_tau_min );
        const std::string density_binwidth = std::to_string( AT_Density_Parameters::density_binwidth );

        std::string header = "";

        Strings();
        ~Strings()
        {}
    };
}

#endif
