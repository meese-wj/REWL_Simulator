#ifndef ORDER_PARAMETER_CUMULANTS
#define ORDER_PARAMETER_CUMULANTS
/* This file will have templated functions 
 * for higher-order combinations of the self-
 * averaged observables. The self-averaged 
 * ones are the moments of the order parameter
 * and the other ones we care about are not
 * self-averaging, for example the suscep-
 * tibility. */

#include <cmath>
#include <cfloat>

// Calculate the EXTENSIVE susceptibility from
// an intensive m and m2. If m and m2 are extensive
// then set size = 1
template<typename obs_t, typename energy_t>
inline obs_t calculate_susceptibility( const obs_t m2, const obs_t m, const energy_t Tvalue, const size_t size )
{
    return (size * m2 - size * size * m * m) / static_cast<obs_t>(Tvalue);
}

// Calculate the Binder parameter (kurtosis) 
// from an intensive m2 and m4. If they are 
// extensive, set size = 1.
template<typename obs_t>
inline obs_t calculate_Binder_parameter( const obs_t m4, const obs_t m2, const size_t size )
{
    return m4 / ( size * m2 * m2 );
}

// Calculate the Binder cumulant of a scalar order parameter
// from an intensive m2 and m4. If they are extensive, set
// size = 1.
template<typename obs_t>
inline obs_t calculate_Binder_cumulant( const obs_t m4, const obs_t m2, const size_t size )
{
    return 1. - calculate_Binder_parameter(m4, m2, size) / 3.;
}

// Calculate the Binder cumulant of a two-component order parameter
// from an intensive m2 and m4. If they are extensive, set
// size = 1.
template<typename obs_t>
inline obs_t calculate_two_component_Binder_cumulant( const obs_t m4, const obs_t m2, const size_t size )
{
    return 2. * ( 1. - 0.5 * ( m4 / ( size * m2 * m2 ) ) );
}

// Calculate the correlation length proxy from the Fourier
// components of the field correlator
template<typename obs_t>
inline obs_t calculate_correlation_length( const obs_t G0, const obs_t Gq, const size_t L )
{
    if ( Gq == 0. )
        return DBL_MAX;

    return ( static_cast<obs_t>(L) / (2. * acos(-1.)) ) * sqrt( (G0 - Gq) / Gq );
}

#endif
