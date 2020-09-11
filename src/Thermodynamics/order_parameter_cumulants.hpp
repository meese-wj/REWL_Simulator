#ifndef ORDER_PARAMETER_CUMULANTS
/* This file will have templated functions 
 * for higher-order combinations of the self-
 * averaged observables. The self-averaged 
 * ones are the moments of the order parameter
 * and the other ones we care about are not
 * self-averaging, for example the suscep-
 * tibility. */

// Calculate the EXTENSIVE susceptibility
template<typename obs_t>
inline obs_t calculate_susceptibility( const obs_t m2, const obs_t m, const energy_t Tvalue )
{
    return (m2 - m * m) / static_cast<obs_t>(Tvalue);
}

// Calculate the Binder parameter (kurtosis)
template<template obs_t>
inline obs_t calculate_Binder_parameter( const obs_t m4, const obs_t m2 )
{
    return m4 / ( m2 * m2 );
}

// Calculate the Binder cumulant of a scalar order parameter
template<template obs_t>
inline obs_t calculate_Binder_parameter( const obs_t m4, const obs_t m2 )
{
    return 1. - calculate_Binder_parameter(m4, m2) / 3.;
}


#endif
