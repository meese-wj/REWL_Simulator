#ifndef VECTOR_MATCHING
#define VECTOR_MATCHING
/* These functions are to help with
 * the concatenation functions. */

#include <numerical_derivatives.hpp>

// Find the indices that starts the overlap 
// in the left energy windows
template<typename data_t>
void find_left_indices( std::vector<size_t> & left_indices, const std::vector<std::vector<data_t> > & energy_table )
{
    const size_t num_windows = energy_table.size();

    for( size_t window = 0; window != num_windows - 1; ++window )
    {
        data_t right_energy = energy_table[window + 1].front();

        for ( size_t left_bin = energy_table[window].size() - 1, num_bins = energy_table[window].size(); left_bin != num_bins; --left_bin )
        {
            if ( energy_table[window][left_bin] == right_energy )
            {
                left_indices[window] = left_bin;
                break;
            }
        }
    }
}

// Find the indices of the concatenation
// in the leftward energy windows via
// numerical differentiation of the 
// logDoS
template<typename energy_t, typename logdos_t>
void find_left_concatenation_indices( std::vector<size_t> & left_concatenate_indices,
                                      const std::vector<size_t> & left_start_indices,
                                      const std::vector<std::vector<logdos_t> > & logdos_table,
                                      const std::vector<std::vector<energy_t> > & energy_table )
{
    const size_t num_windows = energy_table.size();

    for( size_t window = 0; window != num_windows - 1; ++window )
    {
        const size_t left_start = left_start_indices[window];
        logdos_t left_beta = derivative_of_vector( left_start, logdos_table[window], energy_table[window] ); 
        logdos_t right_beta = derivative_of_vector( 0, logdos_table[window + 1], energy_table[window + 1] ); 
        logdos_t min_beta_diff = abs( left_beta - right_beta );
        left_concatenate_indices[window] = 0;

        // Go through the overlapping bins and 
        // find the minimum beta difference
        for ( size_t bin = 1, num_bins = energy_table[window].size(); bin != num_bins; ++bin )
        {
            left_beta = derivative_of_vector( left_start + bin, logdos_table[window], energy_table[window] );
            right_beta = derivative_of_vector( bin, logdos_table[window + 1], energy_table[window + 1] );
            logdos_t test_diff = abs( left_beta - right_beta );

            if ( test_diff < min_beta_diff )
            {
                // Store only the bin since the left_start_indices are
                // already known
                if ( window > 0 && bin >= left_concatenate_indices[window - 1] )
                {
                    // This if statement guarantees that only two windows
                    // will be contatenated moving from left to right. 
                    // In other words, if three windows overlap, the 
                    // closest concatenation allowed between the first
                    // and third windows is to completely sit on top of
                    // the second.
                    // TODO: This is ugly... Can I combine it in some way?
                    left_concatenate_indices[window] = bin;
                    min_beta_diff = test_diff; 
                }
                else if ( window == 0 )
                {
                    left_concatenate_indices[window] = bin;
                    min_beta_diff = test_diff; 
                }
            }
        }
    } 
}


#endif
