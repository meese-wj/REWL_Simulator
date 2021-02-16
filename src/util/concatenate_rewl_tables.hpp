#ifndef CONCATENATE_REWL_TABLES
#define CONCATENATE_REWL_TABLES
/* This file contains the concatenation 
 * procedures for the energy tables, the 
 * logdos tables, and the observable 
 * tables. 
 * Here, a table is a typedef of the 
 * std::vector<std::vector< ... > >    */

#include <vector_matching.hpp>

template <typename data_t>
using table = std::vector<std::vector<data_t> >;

// Make a density table that is organized by
// windows --> energy bins --> values
template <typename data_t>
using dens_table = std::vector<std::vector<std::vector<data_t> > >;

template<typename data_t>
std::vector<data_t> operator + ( const std::vector<data_t> & a, const std::vector<data_t> & b )
{
    // This assumes a and b have the same dimension
    const size_t sze = a.size();
    std::vector<data_t> output (sze);
    for ( size_t idx = 0; idx != sze; ++idx )
        output[idx] = a[idx] + b[idx];
    return output;
}

template<typename data_t>
std::vector<data_t> operator * ( const data_t scalar, const std::vector<data_t> & a )
{
    const size_t sze = a.size();
    std::vector<data_t> output (a);
    for ( size_t idx = 0; idx != sze; ++idx )
        output[idx] = scalar * a[idx];
    return output;
}

// Compute the overlap averages 
template<typename data_t>
inline data_t average_two_windows( const data_t obs1, const data_t obs2, const size_t num1, const size_t num2 )
{
    // TODO: Weighting by counts seems to systematically favor the lower
    // energy windows...
    // return 0.5 * ( obs1 + obs2 );
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
// for a single bin overlap
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
            // The energies should be identical
            final_energy_values.push_back( energy_table[walker][num_bins - 1] );

            // Add the logdos
            final_logdos_values.push_back( logdos_table[walker][num_bins - 1] + logdos_shifter );

            // Subtract out the bottom of the logdos in the right window
            // and add the value of the shifted logdos in the left window
            // now at the end of the final_logdos_values vector.
            // This fixes the logdos to be equal at the concatenation 
            // point.
            logdos_shifter = -logdos_table[walker + 1][0] + final_logdos_values.back();

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
// for a multiple bin overlap
//
// This should concatenate AFTER the replicas
// from a single window have had their 
// measurements averaged.
template<typename energy_t, typename logdos_t, typename obs_t>
#if AT_DENSITIES
void concatenate_tables_multiple_overlap( const size_t num_obs, 
                                          const size_t counts_index,
                                          const table<energy_t>      & energy_table,
                                          const table<logdos_t>      & logdos_table,
                                          const table<obs_t>         & obs_table,
                                          const dens_table<obs_t>    & density_table,
                                          std::vector<energy_t>      & final_energy_values,
                                          std::vector<logdos_t>      & final_logdos_values,
                                          std::vector<obs_t>         & final_obs_values,
                                          table<obs_t>               & final_density_table)
#else
void concatenate_tables_multiple_overlap( const size_t num_obs, 
                                          const size_t counts_index,
                                          const table<energy_t> & energy_table,
                                          const table<logdos_t> & logdos_table,
                                          const table<obs_t>    & obs_table,
                                          std::vector<energy_t> & final_energy_values,
                                          std::vector<logdos_t> & final_logdos_values,
                                          std::vector<obs_t>    & final_obs_values )
#endif
{
    const size_t num_windows = energy_table.size();
    logdos_t logdos_shifter = 0;

    // Record how many overlapping windows 
    // in each bin
    std::vector<size_t> overlapping_windows;
    
    // Add values for the first window to the final
    // values. Then find the best concatenation with
    // the final windows.
    for ( size_t bin = 0, size = energy_table[0].size(); bin != size; ++bin )
    {
        simply_push_back_vectors<energy_t,
                                 logdos_t,
                                 obs_t> ( 0, bin, num_obs, logdos_shifter, energy_table, logdos_table,
                                          obs_table, final_energy_values, final_logdos_values, 
                                          final_obs_values );
#if AT_DENSITIES
        std::cout << "\nPushing back density table\n";
        final_density_table.push_back( std::vector<obs_t> () );
        final_density_table.back() = density_table[0][bin];
        std::cout << "\nPushed back density table\n";
#endif
 
        overlapping_windows.push_back(1);
    }

    // Now go through and concatenate with the final
    // vectors.
    for ( size_t window = 1; window != num_windows; ++window )
    {
        // Find the start indices for the
        // concatenation and the concatenation
        // indices. The final vectors are always
        // the leftward ones.
        index_pair start_indices;
        index_pair concat_indices;
        find_left_concatenation_indices( &concat_indices, &start_indices, final_logdos_values, final_energy_values,
                                         logdos_table[window], energy_table[window] );
       
        // Concatenate first and then add the
        // leftovers in the right window

        // Concatenate:
        size_t left_idx  = 0;
        size_t right_idx = 0;
        for ( size_t idx = start_indices.left; idx != final_energy_values.size(); ++idx )
        {
            left_idx  = idx;
            right_idx = start_indices.right + (idx - start_indices.left);

            // The energies are already in the final
            // energy vector up until the final energy values

            // Average observables in the overlap
            size_t left_obs_bin  = left_idx  * num_obs;
            size_t right_obs_bin = right_idx * num_obs;
            obs_t obs_value = 0;
            for ( size_t ob = 0; ob != num_obs; ++ob )
            {
                obs_value = average_two_windows<obs_t>( final_obs_values[ left_obs_bin + ob ],
                                                        overlapping_windows[ left_idx ],
                                                        obs_table[window][ right_obs_bin + ob ], 1);
                
                final_obs_values[ left_obs_bin + ob ] = ( obs_value );
            }
            
#if AT_DENSITIES
                std::cout << "\nleft_obs_bin " << left_obs_bin << " density table size " << final_density_table.size() << " final obs size " << final_obs_values.size() << "\n";
                std::cout << "Updating density table in window " << window << " at bin " << idx << "\n";
                final_density_table[ left_idx ] = (1./( 1. + overlapping_windows[ left_idx ] )) * ( (1. * overlapping_windows[ left_idx ]) * final_density_table[ left_idx ] + density_table[window][right_idx] );
#endif
            // Increment the number of overlapping
            // windows at this bin.
            ++overlapping_windows[ left_idx ];
 
            // If left_idx < concat_indices.left, then
            // don't do anything to the final logdos
            if ( left_idx < concat_indices.left )
            {}
            // If left_idx == concat_indices.left, then
            // add the logdos and change the shifter
            else if ( left_idx == concat_indices.left )
            {
                logdos_shifter = final_logdos_values[ left_idx ] - logdos_table[window][right_idx];
            }
            // If left_idx > concat_indices.left, then
            // change the final logdos to the rightward 
            // one
            else 
            {
                final_logdos_values[ left_idx ] = logdos_table[window][right_idx] + logdos_shifter;
            }
        }

        // Push leftovers from right window into
        // the final vectors.
        for ( size_t idx = right_idx + 1, right_size = energy_table[window].size(); idx != right_size; ++idx )
        {
            simply_push_back_vectors<energy_t,
                                     logdos_t,
                                     obs_t> ( window, idx, num_obs, logdos_shifter,
                                              energy_table, logdos_table, obs_table, 
                                              final_energy_values, final_logdos_values, 
                                              final_obs_values );

#if AT_DENSITIES
            std::cout << "\nPushing back density table in the final window\n";
            final_density_table.push_back( std::vector<obs_t> () );
            final_density_table.back() = density_table[window][idx];
            std::cout << "\nPushed back density table in the final window\n";
#endif
 
            overlapping_windows.push_back(1);
        }
    }
}


/*
template<typename energy_t, typename logdos_t, typename obs_t>
void concatenate_tables_multiple_overlap( const size_t num_obs, 
                                          const size_t counts_index,
                                          const table<energy_t> & energy_table,
                                          const table<logdos_t> & logdos_table,
                                          const table<obs_t>    & obs_table,
                                          std::vector<energy_t> & final_energy_values,
                                          std::vector<logdos_t> & final_logdos_values,
                                          std::vector<obs_t>    & final_obs_values )
{
    const size_t num_windows = energy_table.size();
    logdos_t logdos_shifter = 0;

    // Create a vector for the indices in the left
    // window where the overlap first occurs
    std::vector<size_t> left_indices_at_overlap (num_windows - 1);
    find_left_indices<energy_t>(left_indices_at_overlap, energy_table);

    // Create a vector for the indices after the 
    // left starter index where the concatenation
    // will take place
    std::vector<size_t> left_concatenation_indices (num_windows - 1);
    find_left_concatenation_indices<energy_t, logdos_t>( left_concatenation_indices,
                                                         left_indices_at_overlap,
                                                         logdos_table, energy_table );

    // Add values for the first window's first bin
    simply_push_back_vectors<energy_t, 
                             logdos_t, 
                             obs_t> ( 0, 0, num_obs, logdos_shifter,
                                      energy_table, logdos_table, obs_table, 
                                      final_energy_values, final_logdos_values,
                                      final_obs_values); 

    // Now go through the overlapping regions. Save
    // the end of the last window for last (this one
    // presumably has no overlaps).
    size_t previous_bin_overlap = 0;
    size_t right_bin = 0;
    
    for ( size_t window = 0; window != num_windows - 1; ++window )
    {
        const size_t left_concat_bin = left_indices_at_overlap[window] + left_concatenation_indices[window];

        // Start from the bin above the last overlapping one
        for ( size_t bin = 0, num_bins = energy_table[window].size() - (previous_bin_overlap + 1); bin != num_bins; ++bin )
        {
            size_t left_bin = bin + previous_bin_overlap + 1;
            if ( previous_bin_overlap > left_concat_bin )
            {
                // This should never happen as long as there is the 
                // bin >= left_concatenate_indices[window - 1] conditional
                // in find_left_concatenation_indices(...)
                printf("\n\nERROR: MISSING DATA IN OVERLAPPING REGIONS FOR WINDOWS %ld AND %ld\n\n", window, window + 1);
            }
                
            // First if there is no overlap, then simply add
            if ( left_bin < left_indices_at_overlap[window] )
            {
                 simply_push_back_vectors<energy_t, 
                                          logdos_t, 
                                          obs_t> ( window, left_bin, num_obs, logdos_shifter,
                                                   energy_table, logdos_table, obs_table, 
                                                   final_energy_values, final_logdos_values,
                                                   final_obs_values);  
            }
            // Now add the overlapping region stuff
            else 
            {
                // TODO: This else guarantees that only left-to-right averaging occurs once.
                // I.e. if three windows overlap, then only the lower two are averaged
                // and the third is lost in the overlap region. This might be a problem
                // if the concatenation for the higher two windows occurs in this 
                // discarded region. If the overlapping regions are sufficiently small,
                // then this will not happen.
                
                // Calculate the right_bin 
                right_bin = left_bin - left_indices_at_overlap[window];

                // Energies should be identical
                final_energy_values.push_back( energy_table[window][left_bin] );

                // Take the leftward logdos up until
                // the concatenation point, and then
                // take the rightward logdos thereafter
                if ( left_bin <= left_concat_bin )
                    final_logdos_values.push_back( logdos_table[window][left_bin] + logdos_shifter );
                else
                    final_logdos_values.push_back( logdos_table[window + 1][right_bin] + logdos_shifter );

                // At the concatenation point, change the logdos_shifter
                // after the leftward logdos was added to the final logdos
                if ( left_bin == left_concat_bin )
                {
                    logdos_shifter = -logdos_table[window + 1][left_concatenation_indices[window]] + final_logdos_values.back();
                }


                // Average observables throughout overlap
                size_t left_obs_bin  = left_bin  * num_obs;
                size_t right_obs_bin = right_bin * num_obs;
                obs_t obs_value = 0;
                for ( size_t ob = 0; ob != num_obs; ++ob )
                {
                    if ( ob != counts_index )
                        obs_value = average_two_windows<obs_t>( obs_table[window    ][left_obs_bin  + ob],
                                                                obs_table[window + 1][right_obs_bin + ob],
                                                                obs_table[window    ][left_obs_bin + counts_index],
                                                                obs_table[window + 1][right_obs_bin + counts_index]);
                    else  // Save the total number of measurements in the overlap. Not sure if this is better or worse than averaging?
                        obs_value = 0.5 * ( obs_table[window][left_obs_bin + counts_index] + obs_table[window + 1][right_obs_bin + counts_index] );

                    final_obs_values.push_back( obs_value );
                } 
            }    
        }     
        // Finally, set the bin of the last overlapping window
        previous_bin_overlap = right_bin;
    }

    // After all of the above, we must now simply add
    // the remainder of the last window
    for ( size_t bin = previous_bin_overlap + 1, num_bins = energy_table[num_windows - 1].size(); bin != num_bins; ++bin )
    {
        simply_push_back_vectors<energy_t, 
                                 logdos_t, 
                                 obs_t> ( num_windows - 1, bin, num_obs, logdos_shifter,
                                          energy_table, logdos_table, obs_table, 
                                          final_energy_values, final_logdos_values,
                                          final_obs_values); 
    }

}
*/

// Concatenate the energy, logdos, and 
// observable tables after the REWL process
// depending on the degree of overlap
template<typename energy_t, typename logdos_t, typename obs_t>
#if AT_DENSITIES
void concatenate_tables( const bool single_overlap,
                         const size_t num_obs, 
                         const size_t counts_index,
                         const table<energy_t>      & energy_table,
                         const table<logdos_t>      & logdos_table,
                         const table<obs_t>         & obs_table,
                         const dens_table<obs_t>    & density_table,
                         std::vector<energy_t>      & final_energy_values,
                         std::vector<logdos_t>      & final_logdos_values,
                         std::vector<obs_t>         & final_obs_values,
                         table<obs_t>               & final_density_table )

#else
void concatenate_tables( const bool single_overlap,
                         const size_t num_obs, 
                         const size_t counts_index,
                         const table<energy_t> & energy_table,
                         const table<logdos_t> & logdos_table,
                         const table<obs_t>    & obs_table,
                         std::vector<energy_t> & final_energy_values,
                         std::vector<logdos_t> & final_logdos_values,
                         std::vector<obs_t>    & final_obs_values )
#endif
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
#if AT_DENSITIES
            concatenate_tables_multiple_overlap<energy_t,
                                                logdos_t,
                                                obs_t>( num_obs, counts_index, energy_table, logdos_table, obs_table, density_table,
                                                        final_energy_values, final_logdos_values, final_obs_values, final_density_table );
#else            
            concatenate_tables_multiple_overlap<energy_t,
                                                logdos_t,
                                                obs_t>( num_obs, counts_index, energy_table, logdos_table, obs_table,
                                                        final_energy_values, final_logdos_values, final_obs_values );
#endif
            break;
        }
    }
    
    return;
}

#endif
