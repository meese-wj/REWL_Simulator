#ifndef REWL_HISTOGRAMS

/*  This header sets up an SOA for the relevant
    REWL histograms for a single walker.  */

#include <cmath>

template<typename data_type>
struct rewl_histogram
{
    data_type min_value;        // Greatest lower bound of the REWL histogram (included).
    data_type max_value;        // Least upper bound of the REWL histogram (not included).
    data_type bin_size;         // Size of the bins for histograms. (This must be the same
                                // across all walkers!).
    
    const size_t num_bins = static_cast<size_t> ((max_value - min_value) / bin_size);

    // Declare the arrays.
    
    ~rewl_histogram()
    {

    }
};


#endif
