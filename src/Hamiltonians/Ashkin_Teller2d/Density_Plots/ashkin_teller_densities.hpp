#ifndef AT_DENSITY_FUNCTIONS
#define AT_DENSITY_FUNCTIONS
/* This file will set up the 2d histograms 
 * in the sigma - tau plane.  It is a 
 * microcanonical observable, and so there
 * will be a 2d histogram for EACH energy 
 * bin. */

#include "ashkin_teller_densities_parameters.cxx"

constexpr density_int energy_bin_density_pointer( const density_int bin )
{
    /* Return the pointer to the first element of the density
     * plots corresponding to the energy bin. */
    return bin * AT_Density_Parameters::total_bins;
}

constexpr density_int density_axis_index( const density_float observable )
{
    /* Return the index along one axis for the given observable. */
    return static_cast<density_int>( ( observable - AT_Density_Parameters::sigma_tau_min ) / AT_Density_Parameters::density_binwidth ); 
}

constexpr density_int density_index( const density_float sigma, const density_float tau )
{
    /* Return the index in the density array corresponding to sigma and tau. */
   return density_axis_index( tau ) * AT_Density_Parameters::axis_bins + density_axis_index( sigma );
}

inline void update_densities( density_int * const density_data, const density_int bin, const density_float sigma, const density_float tau )
{
    /* Update the density plot for a given energy bin. */
    ++density_data[ energy_bin_density_pointer( bin ) + density_index( sigma, tau ) ];    
}

void transfer_and_normalize( const density_int nbins, const density_int * const density_data, density_float * const float_data  )
{
    /* Populate float_data with the normalized density_data.
     * The normalization is with respect to the maximum value
     * of the density data. */
    density_int * max_values = new density_int [nbins]; 
    for ( density_int bin = 0; bin != nbins; ++bin )
    {
        max_values[ bin ] = 0.;
        // Scan the density plot for the maximum value
        for ( density_int density_bin = 0; density_bin != AT_Density_Parameters::total_bins; ++density_bin  )
        {
            max_values[ bin ] = ( max_values[ bin ] < density_data[ energy_bin_density_pointer( bin ) + density_bin ] ? density_data[ energy_bin_density_pointer( bin ) + density_bin ] : max_values[ bin ] );
        }
        // Transfer to the float_data for this specific energy bin
        for ( density_int density_bin = 0; density_bin != AT_Density_Parameters::total_bins; ++density_bin  )
        {
            float_data[ energy_bin_density_pointer( bin ) + density_bin ] = static_cast<density_float>( density_data[ energy_bin_density_pointer( bin ) + density_bin ] ) / static_cast<density_float>( max_values[ bin ] );
        }
    }
    delete [] max_values;
}

#endif
