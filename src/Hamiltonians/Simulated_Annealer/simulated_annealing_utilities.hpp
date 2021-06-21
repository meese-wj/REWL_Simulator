#ifndef SIMULATED_ANNEALING_UTILITIES
#define SIMULATED_ANNEALING_UTILITIES

#include <cmath>
#include <cstdint>
#include "simulated_annealing_parameters.cxx"

template<typename energy_t>
inline bool not_frozen( const energy_t value )
{
    return value != SA_Parameters::frozen_criterion;
}

template<typename energy_t>
inline bool not_within_error( const energy_t value, const energy_t avg, const energy_t stdev )
{
    return !( avg - stdev <= value && value <= avg + stdev );
}

template<typename energy_t>
inline bool not_within_tolerance( const energy_t new_stdev, const energy_t old_stdev, const energy_t tolerance )
{
    return ( abs(new_stdev - old_stdev) / old_stdev > tolerance );
}

template<typename energy_t>
void array_statistics( const uint32_t arr_size, const energy_t * const array, 
                       energy_t * const average, energy_t * const stdev )
{
    *average = 0.;
    *stdev = 0.;
    // First calculate the average of the array
    for ( uint32_t idx = 0; idx != arr_size; ++idx )
    {
        *average += array[idx];
    }
    *average /= static_cast<energy_t>( arr_size );

    // Now calculate the standard deviation
    for ( uint32_t idx = 0; idx != arr_size; ++idx )
    {
        *stdev += (array[idx] - *average) * (array[idx] - *average);
    }
    *stdev = sqrt( *stdev / static_cast<energy_t>( arr_size ) );
}


#endif
