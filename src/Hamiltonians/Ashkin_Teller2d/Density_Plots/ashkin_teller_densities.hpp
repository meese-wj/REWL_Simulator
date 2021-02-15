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

constexpr density_int density_axis_index( const float observable )
{
    /* Return the index along one axis for the given observable. */
    return static_cast<density_int>( ( observable - AT_Density_Parameters::sigma_tau_min ) / AT_Density_Parameters::density_binwidth ); 
}

constexpr density_int density_index( const float sigma, const float tau )
{
    /* Return the index in the density array corresponding to sigma and tau. */
   return density_axis_index( tau ) * AT_Density_Parameters::axis_bins + density_axis_index( sigma );
}

inline void update_densities( density_int * const density_data, const density_int bin, const float sigma, const float tau )
{
    /* Update the density plot for a given energy bin. */
    ++density_data[ energy_bin_density_pointer( bin ) + density_index( sigma, tau ) ];    
}

#endif
