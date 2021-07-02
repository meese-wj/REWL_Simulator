#ifndef ARRAY_SHIFT
#define ARRAY_SHIFT
/* This function shifts an array by
 * a common scalar value. */
#include <cmath>

template<typename data_t>
void array_shift_by_value( const data_t value, const size_t size, data_t * const data_array )
{
   for ( size_t idx = 0; idx != size; ++idx )
       data_array[ idx ] += value;
}

template<typename data_t, data_t energy_min, data_t energy_max = 0.>
void normalize_logDoS( const data_t states_per_site, const data_t system_size, const size_t size, data_t * const logDoS_array )
{
    static_assert( energy_max == 0. || energy_max == -energy_min, "Either normalize by half of a symmetric energy spectrum or the whole asymmetric spectrum." );

    // First determine the maximum logDoS value
    data_t logDoS_max = 0.;
    for ( size_t bin = 0; bin != size; ++bin )
        logDoS_max = logDoS_max * ( logDoS_max > logDoS_array[bin] ) + logDoS_array[bin] * ( logDoS_max <= logDoS_array[bin] );

    // Next evaluate the sum of exponential
    // differences between the logDoS and the
    // max. Add smallest to largest for better
    // error.
    data_t logDoS_sum = 0.;
    for ( size_t bin = 0; bin != size; ++bin )
        logDoS_sum += exp( -( logDoS_max - logDoS_array[bin] ) );
    // Add the max logDoS
    logDoS_sum += logDoS_max;

    // Now shift the logDoS array according
    // to this normalizer.
    data_t logDoS_shifter = static_cast<data_t>(system_size) * log( states_per_site ) - logDoS_sum - log(2.) * ( energy_max == 0. ); 
    for ( size_t bin = 0; bin != size; ++bin )
        logDoS_array[bin] += logDoS_shifter;
}


#endif
