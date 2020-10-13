#ifndef CONCATENATE_REWL_TABLES
#define CONCATENATE_REWL_TABLES
/* This file contains the concatenation 
 * procedures for the energy tables, the 
 * logdos tables, and the observable 
 * tables. 
 * Here, a table is a typedef of the 
 * std::vector<std::vector< ... > >    */

#include <numerical_derivatives.hpp>

template <typename data_t>
using table = std::vector<std::vector<data_t> >;

// Compute the overlap averages 
template<typename data_t>
inline data_t average_two_windows( const data_t obs1, const data_t obs2, const size_t num1, const size_t num2 )
{
    return ( static_cast<data_t>(num1) * obs1 + static_cast<data_t>(num2) * obs2 ) / static_cast<data_t>( num1 + num2 );
}

// Push back the final values simply
// as in there are no overlaps to deal with
template<typename energy_t, typename logdos_t, typename obs_t>
void simply_push_back_vectors( const size_t walker, 
                               const size_t bin,
                               const size_t num_obs,
                               const logdos_t logdos_shifter,
                               const table<energy_t> & energy_table,
                               const table<logdos_t> & logdos_table,
                               const table<obs_t>    & obs_table,
                               std::vector<energy_t> & final_energy_values,
                               std::vector<logdos_t> & final_logdos_values,
                               std::vector<obs_t>    & final_obs_values )
{
    final_energy_values.push_back( energy_table[walker][bin] );
    final_logdos_values.push_back( logdos_table[walker][bin] + logdos_shifter );

    const size_t obs_index = bin * num_obs;
    for ( size_t ob = 0; ob != num_obs; ++ob )
        final_obs_values.push_back( obs_table[walker][obs_index + ob] );
}

// Concatenate the energy, logdos, and 
// observable tables after the REWL process
template<typename energy_t, typename logdos_t, typename obs_t>
void concatenate_tables_single_overlap( const size_t num_obs, 
                                        const size_t counts_index,
                                        const table<energy_t> & energy_table,
                                        const table<logdos_t> & logdos_table,
                                        const table<obs_t>    & obs_table,
                                        std::vector<energy_t> & final_energy_values,
                                        std::vector<logdos_t> & final_logdos_values,
                                        std::vector<obs_t>    & final_obs_values )
{
    const size_t num_walkers = energy_table.size();
    logdos_t logdos_shifter = 0;

    // Add values for the first walker's first bin
    simply_push_back_vectors<energy_t, 
                             logdos_t, 
                             obs_t> ( 0, 0, num_obs, logdos_shifter,
                                      energy_table, logdos_table, obs_table, 
                                      final_energy_values, final_logdos_values,
                                      final_obs_values); 

    for ( size_t walker = 0; walker != num_walkers; ++walker )
    {
        const size_t num_bins = energy_table[walker].size();
        // Fill up the tables for each walker up until the last bin
        // (the first bin is used in the concatenation)
        for ( size_t bin = 1; bin != num_bins - 1; ++bin )
        {
           simply_push_back_vectors<energy_t, 
                                    logdos_t, 
                                    obs_t> ( walker, bin, num_obs, logdos_shifter,
                                             energy_table, logdos_table, obs_table, 
                                             final_energy_values, final_logdos_values,
                                             final_obs_values);             
        }

        // Now for all but the last walker, concatenate at the last bin
        if ( walker != num_walkers - 1 )
        {
            // Subtract out the bottom of the logdos in the right window
            // and add the value of the logdos in the left window.
            // This fixes the logdos to be equal at the concatenation 
            // point.
            logdos_shifter = -logdos_table[walker + 1][0] + logdos_table[walker][num_bins - 1];

            // The energies should be identical
            final_energy_values.push_back( energy_table[walker][num_bins - 1] );

            // Add the logdos
            final_logdos_values.push_back( logdos_table[walker + 1][0] + logdos_shifter );

            // Average observables in the overlap by their counts
            size_t left_obs_bin  = ( num_bins - 1 ) * num_obs;
            size_t right_obs_bin = ( 0 ) * num_obs;
            obs_t obs_value = 0;
            for ( size_t ob = 0; ob != num_obs; ++ob )
            {
                if ( ob != counts_index )
                    obs_value = average_two_windows<obs_t>( obs_table[walker    ][left_obs_bin  + ob],
                                                            obs_table[walker + 1][right_obs_bin + ob],
                                                            obs_table[walker    ][left_obs_bin + counts_index],
                                                            obs_table[walker + 1][right_obs_bin + counts_index]);
                else  // Save the total number of measurements in the overlap. Not sure if this is better or worse than averaging?
                    obs_value = obs_table[walker][left_obs_bin + counts_index] + obs_table[walker + 1][right_obs_bin + counts_index];

                final_obs_values.push_back( obs_value );
            }

        }        
        else // For the last walker, include the last bin
        {
            simply_push_back_vectors<energy_t, 
                                     logdos_t, 
                                     obs_t> ( walker, num_bins - 1,
                                              num_obs, logdos_shifter,
                                              energy_table, logdos_table, obs_table, 
                                              final_energy_values, final_logdos_values,
                                              final_obs_values); 
        }

    }
}

// Concatenate the energy, logdos, and 
// observable tables after the REWL process
// depending on the degree of overlap
template<typename energy_t, typename logdos_t, typename obs_t>
void concatenate_tables( const bool single_overlap,
                         const size_t num_obs, 
                         const size_t counts_index,
                         const table<energy_t> & energy_table,
                         const table<logdos_t> & logdos_table,
                         const table<obs_t>    & obs_table,
                         std::vector<energy_t> & final_energy_values,
                         std::vector<logdos_t> & final_logdos_values,
                         std::vector<obs_t>    & final_obs_values )
{

    switch( single_overlap )
    {
        case true:
        {
            concatenate_tables_single_overlap<energy_t,
                                              logdos_t,
                                              obs_t>( num_obs, counts_index, energy_table, logdos_table, obs_table,
                                                      final_energy_values, final_logdos_values, final_obs_values );
            break;
        }
        case false:
        {
            // TODO: Fill this in for non singular overlaps
            break;
        }
    }
    
    return;
}

#endif
