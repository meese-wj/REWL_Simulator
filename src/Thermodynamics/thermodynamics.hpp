#ifndef THERMODYNAMICS
#define THERMODYNAMICS
/* Create an object to calculate the 
 * canonical thermodynamics from a density
 * of states using trapezoidal rule. */

#include <cmath>
#include <string>
#include <vector>

template<class enum_t>
constexpr size_t convert( const enum_t ob )
{
    return static_cast<size_t>(ob);
}

namespace Energy_Obs
{
    enum class enum_names
    {
        free_energy, internal_energy, internal_energy2, entropy, specific_heat, NUM_OBS
    };

    const std::vector<std::string> string_names = { "Free Energy", "Energy", "Energy2", "Entropy", "Specific Heat", "NUM OBS" };
}

template<typename energy_t, typename logdos_t, typename obs_t, class Obs_enum_t>
struct Thermodynamics
{
    const size_t total_observables = convert(Energy_Obs::enum_names::NUM_OBS) + convert(Obs_enum_t::NUM_OBS);

    const size_t num_energy_bins;

    const energy_t Tmin;
    const energy_t Tmax; 
    const size_t num_T;
    const energy_t dT = (Tmax - Tmin) / static_cast<energy_t>(num_T);
    energy_t * temperatures = nullptr;

    obs_t * canonical_observables = nullptr;

    Thermodynamics(const size_t _numEbins,
                   const energy_t _Tmin, const energy_t _Tmax, const size_t _nT) : num_energy_bins(_numEbins),
                                                                                   Tmin(_Tmin), Tmax(_Tmax), num_T(_nT)
    {
        temperatures = new energy_t [ num_T ];
        for ( size_t Tidx = 0; Tidx != num_T; ++Tidx )
            temperatures[ Tidx ] = Tmin + dT * static_cast<energy_t>(Tidx);

        canonical_observables = new obs_t [ num_T * total_observables ];
        for ( size_t idx = 0; idx != num_T * total_observables; ++idx )
            canonical_observables[ idx ] = 0.;
    }

    ~Thermodynamics()
    { 
        delete [] temperatures;
        delete [] canonical_observables;
    }

    // Create a getter function to copy out the thermally-averaged observables
    obs_t get_energy_obs( const size_t Tidx, const Energy_Obs::enum_names ob ) const { return canonical_observables[ Tidx * total_observables + convert(ob) ]; }
    obs_t get_system_obs( const size_t Tidx, const size_t ob ) const { return canonical_observables[ Tidx * total_observables + convert(Energy_Obs::enum_names::NUM_OBS) + ob ]; }

    // Create an alias function to the thermally-averaged energy observables
    obs_t * energy_obs( const size_t Tidx, const Energy_Obs::enum_names ob ) const
    { 
        return &canonical_observables[ Tidx * total_observables + convert(ob) ]; 
    }
  
    // Create an alias function to the thermally-averaged system observables
    obs_t * system_obs( const size_t Tidx, const size_t ob ) const
    { 
        return &canonical_observables[ Tidx * total_observables + convert(Energy_Obs::enum_names::NUM_OBS) + ob ]; 
    }
 
    // Calculate the exponent in the partition function
    obs_t get_exponent( const size_t bin, const energy_t Tvalue, const energy_t * const energy_array, const logdos_t * const logdos_array ) const
    {
        return static_cast<obs_t> ( logdos_array[ bin ] - static_cast<logdos_t>( energy_array[ bin ] / Tvalue ) );
    }
    
    // Calculate the exponent in the partition function shifted by its maximum
    obs_t get_reduced_exponent( const size_t bin, const energy_t Tvalue, const obs_t maximum, 
                                const energy_t * const energy_array, const logdos_t * const logdos_array ) const
    {
        return ( static_cast<obs_t> ( get_exponent(bin, Tvalue, energy_array, logdos_array) - maximum ) );
    }

    obs_t reduced_exponential( const size_t bin, const energy_t Tvalue, const obs_t maximum, 
                               const energy_t * const energy_array, const logdos_t * const logdos_array ) const
    {
        return exp( get_reduced_exponent( bin, Tvalue, maximum, energy_array, logdos_array ) );
    }
   
    // Calculate the partition function reduced by the 
    // exponential of the maximum exponent.
#if TRAPEZOIDAL_RULE
    obs_t reduced_partition_function( const energy_t Tvalue, const obs_t maximum, 
                                      const energy_t * const energy_array, const logdos_t * const logdos_array ) const
    {
        obs_t part = 0.;
        for ( size_t bin = 1; bin <= num_energy_bins; ++bin )
        {
            const size_t current_bin = num_energy_bins - bin;
            if ( current_bin != 0 && current_bin != num_energy_bins - 1 )
                part += 2. * reduced_exponential( current_bin, Tvalue, maximum, energy_array, logdos_array );
            else
                part += reduced_exponential( current_bin, Tvalue, maximum, energy_array, logdos_array );
        }
        return 0.5 * part;
    }
#else
    obs_t reduced_partition_function( const energy_t Tvalue, const obs_t maximum, 
                                      const energy_t * const energy_array, const logdos_t * const logdos_array ) const
    {
        obs_t part = 0.;
        for ( size_t bin = 1; bin <= num_energy_bins; ++bin )
        {
            const size_t current_bin = num_energy_bins - bin;
            part += reduced_exponential( current_bin, Tvalue, maximum, energy_array, logdos_array );
        }
        return part;
    }
#endif
    
    // Calculate the maximum exponent for a given temperature
    obs_t find_maximum_exponent( const energy_t Tvalue, 
                                 const energy_t * const energy_array, const logdos_t * const logdos_array ) const;

    // Calculate the thermodynamics (both energy and self-averaging observables)
    void calculate_thermodynamics( const size_t system_size,
                                   const energy_t * const energy_array,
                                   const logdos_t * const logdos_array,
                                   const obs_t * const observables_array ) const;  
};

// Calculate the maximum microcanonical exponent as a function
// of input temperature.
template<typename energy_t, typename logdos_t, typename obs_t, class Obs_enum_t>
obs_t Thermodynamics<energy_t, logdos_t, obs_t, Obs_enum_t>::find_maximum_exponent( const energy_t Tvalue, 
                                                                                    const energy_t * const energy_array,
                                                                                    const logdos_t * const logdos_array ) const
{
    obs_t max = get_exponent(num_energy_bins - 1, Tvalue, energy_array, logdos_array);
    obs_t value = max;
    for ( size_t bin = 2; bin <= num_energy_bins; ++bin )
    {
        const size_t current_bin = num_energy_bins - bin;
        value = get_exponent(current_bin, Tvalue, energy_array, logdos_array);
        max = ( max < value ? value : max );
    }
    return max;
}

// Calculate all the thermodynamics.
#if TRAPEZOIDAL_RULE
template<typename energy_t, typename logdos_t, typename obs_t, class Obs_enum_t>
void Thermodynamics<energy_t, logdos_t, obs_t,
                    Obs_enum_t>::calculate_thermodynamics( const size_t system_size,
                                                           const energy_t * const energy_array,
                                                           const logdos_t * const logdos_array,
                                                           const obs_t * const observables_array ) const
{
    for ( size_t Tidx = 1; Tidx <= num_T; ++Tidx )
    {
        const size_t current_Tidx = num_T - Tidx;
        energy_t Tvalue = temperatures[ current_Tidx ];
        
        // First find the maximum exponent
        const obs_t max_exponent = find_maximum_exponent( Tvalue, energy_array, logdos_array );
 
        // Second compute the partition function with the maximum exponent scaled out
        const obs_t partition = reduced_partition_function( Tvalue, max_exponent, energy_array, logdos_array );
        
        // Third compute the energy observables
        *energy_obs( current_Tidx, Energy_Obs::enum_names::free_energy ) = -Tvalue * ( max_exponent + log(partition) ) / static_cast<obs_t>(system_size);

        for ( size_t bin = 1; bin <= num_energy_bins; ++bin )
        {
            const size_t current_bin = num_energy_bins - bin;
            obs_t energy_value = energy_array[ current_bin ]; 
            obs_t weight = reduced_exponential( current_bin, Tvalue, max_exponent, energy_array, logdos_array );
            
            if ( current_bin == 0 || current_bin == num_energy_bins - 1 )
            {
                *energy_obs( current_Tidx, Energy_Obs::enum_names::internal_energy ) += 0.5 * energy_value * weight;
                *energy_obs( current_Tidx, Energy_Obs::enum_names::internal_energy2 ) += 0.5 * energy_value * energy_value * weight;
                *energy_obs( current_Tidx, Energy_Obs::enum_names::entropy ) += 0.5 * logdos_array[ current_bin ] * weight;
            }
            else 
            {
                *energy_obs( current_Tidx, Energy_Obs::enum_names::internal_energy ) += energy_value * weight;
                *energy_obs( current_Tidx, Energy_Obs::enum_names::internal_energy2 ) += energy_value * energy_value * weight;
                *energy_obs( current_Tidx, Energy_Obs::enum_names::entropy ) += logdos_array[ current_bin ] * weight;
            }

            // Compute the extra observables
            for ( size_t ob = 0; ob != convert(Obs_enum_t::NUM_OBS); ++ob  )
            {
                if ( ob != convert(Obs_enum_t::counts_per_bin) )
                { 
                    if ( current_bin == 0 || current_bin == num_energy_bins - 1 )
                    {
                        *system_obs( current_Tidx, ob ) += 0.5 * ( observables_array[ current_bin * convert(Obs_enum_t::NUM_OBS) + ob ] ) * weight;     
                    }
                    else
                    {
                        *system_obs( current_Tidx, ob ) += ( observables_array[ current_bin * convert(Obs_enum_t::NUM_OBS) + ob ] ) * weight;     
                    }
                }
                else
                    *system_obs( current_Tidx, ob ) += observables_array[ current_bin * convert(Obs_enum_t::NUM_OBS) + ob ];
            }
        }
        
        // Normalize the result
        *energy_obs( current_Tidx, Energy_Obs::enum_names::internal_energy ) /= ( system_size * partition );
        *energy_obs( current_Tidx, Energy_Obs::enum_names::internal_energy2 ) /= ( system_size * system_size * partition );
        *energy_obs( current_Tidx, Energy_Obs::enum_names::entropy ) /= ( system_size * partition );

        // Compute the specific heat
        obs_t En = get_energy_obs( current_Tidx, Energy_Obs::enum_names::internal_energy );
        obs_t En2 = get_energy_obs( current_Tidx, Energy_Obs::enum_names::internal_energy2 );
        *energy_obs( current_Tidx, Energy_Obs::enum_names::specific_heat ) = system_size * (En2 - En * En) / static_cast<obs_t>( Tvalue * Tvalue );


        for ( size_t ob = 0; ob != convert(Obs_enum_t::NUM_OBS); ++ob )
        {
            if ( ob != convert(Obs_enum_t::counts_per_bin) )
                *system_obs( current_Tidx, ob ) /= ( partition );
            else
                *system_obs( current_Tidx, ob ) /= num_energy_bins;
        }        
    }
}
#else
template<typename energy_t, typename logdos_t, typename obs_t, class Obs_enum_t>
void Thermodynamics<energy_t, logdos_t, obs_t,
                    Obs_enum_t>::calculate_thermodynamics( const size_t system_size,
                                                           const energy_t * const energy_array,
                                                           const logdos_t * const logdos_array,
                                                           const obs_t * const observables_array ) const
{
    for ( size_t Tidx = 1; Tidx <= num_T; ++Tidx )
    {
        const size_t current_Tidx = num_T - Tidx;
        energy_t Tvalue = temperatures[ current_Tidx ];
        
        // First find the maximum exponent
        const obs_t max_exponent = find_maximum_exponent( Tvalue, energy_array, logdos_array );
 
        // Second compute the partition function with the maximum exponent scaled out
        const obs_t partition = reduced_partition_function( Tvalue, max_exponent, energy_array, logdos_array );
        
        // Third compute the energy observables
        *energy_obs( current_Tidx, Energy_Obs::enum_names::free_energy ) = -Tvalue * ( max_exponent + log(partition) ) / static_cast<obs_t>(system_size);

        for ( size_t bin = 1; bin <= num_energy_bins; ++bin )
        {
            const size_t current_bin = num_energy_bins - bin;
            obs_t energy_value = energy_array[ current_bin ]; 
            obs_t weight = reduced_exponential( current_bin, Tvalue, max_exponent, energy_array, logdos_array );
            
            *energy_obs( current_Tidx, Energy_Obs::enum_names::internal_energy ) += energy_value * weight;
            *energy_obs( current_Tidx, Energy_Obs::enum_names::internal_energy2 ) += energy_value * energy_value * weight;
            *energy_obs( current_Tidx, Energy_Obs::enum_names::entropy ) += logdos_array[ current_bin ] * weight;
            
            // Compute the extra observables
            for ( size_t ob = 0; ob != convert(Obs_enum_t::NUM_OBS); ++ob  )
            {
                if ( ob != convert(Obs_enum_t::counts_per_bin) )
                    *system_obs( current_Tidx, ob ) += ( observables_array[ current_bin * convert(Obs_enum_t::NUM_OBS) + ob ] ) * weight;
                else
                    *system_obs( current_Tidx, ob ) += observables_array[ current_bin * convert(Obs_enum_t::NUM_OBS) + ob ];
            }
        }
        
        // Normalize the result
        *energy_obs( current_Tidx, Energy_Obs::enum_names::internal_energy ) /= ( system_size * partition );
        *energy_obs( current_Tidx, Energy_Obs::enum_names::internal_energy2 ) /= ( system_size * system_size * partition );
        *energy_obs( current_Tidx, Energy_Obs::enum_names::entropy ) /= ( system_size * partition );

        // Compute the specific heat
        obs_t En = get_energy_obs( current_Tidx, Energy_Obs::enum_names::internal_energy );
        obs_t En2 = get_energy_obs( current_Tidx, Energy_Obs::enum_names::internal_energy2 );
        *energy_obs( current_Tidx, Energy_Obs::enum_names::specific_heat ) = system_size * (En2 - En * En) / static_cast<obs_t>( Tvalue * Tvalue );

        for ( size_t ob = 0; ob != convert(Obs_enum_t::NUM_OBS); ++ob )
        {
            if ( ob != convert(Obs_enum_t::counts_per_bin) )
                *system_obs( current_Tidx, ob ) /= ( partition );
            else
                *system_obs( current_Tidx, ob ) /= num_energy_bins;
        }        
    }
}
#endif


#endif
