#ifndef THERMODYNAMICS
/* Create an object to calculate the 
 * canonical thermodynamics from a density
 * of states using trapezoidal rule. */

#include <cmath>

template<class enum_t>
constexpr size_t convert( const enum_t ob )
{
    return static_cast<size_t>(ob);
}

enum class Energy_Obs
{
    free_energy, internal_energy, internal_energy2, entropy, specific_heat, NUM_OBS
};

template<typename energy_t, typename logdos_t, typename obs_t, class Obs_enum_t>
struct Thermodynamics
{
    const size_t total_observables = convert(Energy_Obs::NUM_OBS) + convert(Obs_enum_t::NUM_OBS);

    const energy_t energy_min;
    const energy_t energy_max;
    const energy_t energy_bin_size;
    const size_t num_energy_bins = static_cast<size_t>( (energy_max - energy_min) / energy_bin_size );

    const energy_t Tmin;
    const energy_t Tmax; 
    const size_t num_T;
    const energy_t dT = (Tmax - Tmin) / static_cast<energy_t>(num_T);
    energy_t * temperatures = nullptr;

    obs_t * canonical_observables = nullptr;

    Thermodynamics(const energy_t _emin, const energy_t _emax, const energy_t _ebsize,
                   const energy_t _Tmin, const energy_t _Tmax, const size_t _nT) : energy_min(_emin), energy_max(_emax), energy_bin_size(_ebsize),
                                                                                   Tmin(_Tmin), Tmax(_Tmax), num_T(_nT)
    {
        temperatures = new energy_t [ num_T ];
        for ( size_t Tidx = 0; Tidx != num_T; ++Tidx )
            temperatures = Tmin + dT * static_cast<energy_t>(Tidx);

        canonical_observables = new obs_t [ num_T * total_observables ];
        for ( size_t idx = 0; idx != num_T * total_observables; ++idx )
            canonical_observables[ idx ] = 0.;
    }

    ~Thermodynamics()
    { 
        if (temperatures != nullptr) delete [] temperatures;
        if (canonical_observables != nullptr) delete [] canonical_observables;
    }
    
    // Calculate the exponent in the partition function
    obs_t get_exponent( const size_t bin, const energy_t Tvalue, const logdos_t * const logdos_array ) const
    {
        return static_cast<obs_t> ( logdos_array[ bin ] - static_cast<logdos_t>( (energy_min + static_cast<energy_t>(bin) * energy_bin_size) / Tvalue ) );
    }
    
    // Calculate the exponent in the partition function shifted by its maximum
    obs_t get_reduced_exponent( const size_t bin, const energy_t Tvalue, const obs_t maximum, const logdos_t * const logdos_array ) const
    {
        return ( static_cast<obs_t> ( logdos_array[ bin ] - static_cast<logdos_t>( (energy_min + static_cast<energy_t>(bin) * energy_bin_size) / Tvalue ) ) - maximum );
    }

    obs_t reduced_exponential( const size_t bin, const energy_t Tvalue, const obs_t maximum, const logdos_t * const logdos_array ) const
    {
        return exp( get_reduced_exponent( bin, Tvalue, maximum, logdos_array ) );
    }
    
    // Calculate the maximum exponent for a given temperature
    obs_t find_maximum_exponent( const energy_t Tvalue, const logdos_t * const logdos_array ) const;

    // Calculate the thermodynamics (both energy and self-averaging observables)
    void calculate_thermodynamics( const size_t system_size,
                                   const energy_t * const energy_array,
                                   const logdos_t * const logdos_array,
                                   const obs_t * const observables ) const;  
};

// Calculate the maximum microcanonical exponent as a function
// of input temperature.
template<typename energy_t, typename logdos_t, typename obs_t, class Obs_enum_t>
obs_t Thermodynamics<energy_t, logdos_t, obs_t, Obs_enum_t>::find_maximum_exponent( const energy_t Tvalue, const logdos_t * const logdos_array ) const
{
    obs_t max = get_exponent(0, Tvalue, logdos_array);
    obs_t value = max;
    for ( size_t bin = 1; bin != num_energy_bins; ++bin )
    {
        value = get_exponent(bin, Tvalue, logdos_array);
        max = ( max < value ? value : max );
    }
    return max;
}

// Calculate all the thermodynamics.
template<typename energy_t, typename logdos_t, typename obs_t, class Obs_enum_t>
void Thermodynamics<energy_t, logdos_t, obs_t,
                    Obs_enum_t>::calculate_thermodynamics( const size_t system_size,
                                                           const energy_t * const energy_array,
                                                           const logdos_t * const logdos_array,
                                                           const obs_t * const observables_array ) const
{
    obs_t partition = 0.;
    obs_t max_exponent = 0.;
    for ( size_t Tidx = 0; Tidx != num_T; ++Tidx )
    {
        energy_t Tvalue = temperatures[ Tidx ];

        // First find the maximum exponent
        max_exponent = find_maximum_exponent( Tvalue, logdos_array );
        
        // Second compute the partition function with the maximum exponent scaled out
        partition = 0.;
        for ( size_t bin = 0; bin != num_energy_bins; ++bin )
        {
            partition += reduced_exponential( bin, Tvalue, max_exponent, logdos_array);
        }

        // Third compute the energy observables
        canonical_observables[ Tidx * total_observables + convert(Energy_Obs::free_energy) ] = -Tvalue * log( partition ) / static_cast<obs_t> (system_size);
        
        for ( size_t bin = 0; bin != num_energy_bins; ++bin )
        {
            obs_t energy_value = energy_array[ bin ]; 
            obs_t weight = reduced_exponential( bin, Tvalue, max_exponent, logdos_array );

            canonical_observables[ Tidx * total_observables + convert(Energy_Obs::internal_energy) ] += energy_value * weight;
            canonical_observables[ Tidx * total_observables + convert(Energy_Obs::internal_energy2) ] += energy_value * energy_value * weight;
            canonical_observables[ Tidx * total_observables + convert(Energy_Obs::entropy) ] += logdos_array[bin] * weight;
            
            // Compute the extra observables
            for ( size_t ob = 0; ob != convert(Obs_enum_t::NUM_OBS); ++ob  )
            {
                if ( ob != convert(Obs_enum_t::counts_per_bin) )
                    canonical_observables[ Tidx * total_observables + convert(Energy_Obs::NUM_OBS) + ob ] += ( observables_array[ bin * convert(Obs_enum_t::NUM_OBS) + ob ] ) * weight;
                else
                    canonical_observables[ Tidx * total_observables + convert(Energy_Obs::NUM_OBS) + ob ] += ( observables_array[ bin * convert(Obs_enum_t::NUM_OBS) + ob ] );
            }
        }
        // Normalize the result
        partition *= static_cast<obs_t> (system_size);
        canonical_observables[ Tidx * total_observables + convert(Energy_Obs::internal_energy) ] /= partition;
        canonical_observables[ Tidx * total_observables + convert(Energy_Obs::internal_energy2) ] /= partition;
        canonical_observables[ Tidx * total_observables + convert(Energy_Obs::entropy) ] /= partition;

        // Compute the specific heat
        canonical_observables[ Tidx * total_observables + convert(Energy_Obs::specific_heat) ] = ( canonical_observables[ Tidx * total_observables + convert(Energy_Obs::internal_energy2) ] 
                                                                                                 - canonical_observables[ Tidx * total_observables + convert(Energy_Obs::internal_energy) ] 
                                                                                                 * canonical_observables[ Tidx * total_observables + convert(Energy_Obs::internal_energy) ] ) 
                                                                                                 / ( static_cast<obs_t>( system_size * Tvalue * Tvalue ) ); 

        for ( size_t ob = 0; ob != convert(Energy_Obs::NUM_OBS); ++ob )
        {
            canonical_observables[ Tidx * total_observables + convert(Energy_Obs::NUM_OBS) + ob ] /= partition;
            if ( ob != convert(Obs_enum_t::counts_per_bin) )
                canonical_observables[ Tidx * total_observables + convert(Energy_Obs::NUM_OBS) + ob ] /= partition;
            else
                canonical_observables[ Tidx * total_observables + convert(Energy_Obs::NUM_OBS) + ob ] /= num_energy_bins;
        }
        
    }
}



#endif
