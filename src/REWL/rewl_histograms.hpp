#ifndef REWL_HISTOGRAMS

/*  This header sets up an SOA for the relevant
    REWL histograms for a single walker.  

    The data_t should almost always be of type
    64-bit double for memory management reasons.
    It CANNOT be a 32-bit float which does not 
    have a large enough mantissa for the level of
    precision these simulations require.
*/

#include <cmath>
#include <cstdint>

static constexpr size_t num_pair = 2;
static constexpr uint64_t count_initializer = 0;
static constexpr uint64_t logdos_initializer = 1;

template<typename data_t>
struct value_pair
{
    uint64_t count;
    data_t logdos; 
};

template<typename data_t>
struct rewl_histogram
{
    data_t min_value;        // Greatest lower bound of the REWL histogram (included).
    data_t max_value;        // Least upper bound of the REWL histogram (not included).
    data_t bin_size;         // Size of the bins for histograms. (This must be the same
                             // across all walkers!).
    
    const size_t num_bins = static_cast<size_t> ((max_value - min_value) / bin_size);

    // Declare the arrays.
    value_pair<data_t> * histograms = nullptr;

    rewl_histogram(data_t min, data_t max, data_t b_size) : min_value(min), max_value(max), bin_size(b_size)
    {
        histograms = new value_pair [bin_size];
        for ( size_t idx = 0; idx != num_bins; ++idx )
        {
            histograms[ num_pair * idx ] = count_initializer;         // Initialize the counts to zero
            histograms[ num_pair * idx + 1 ] = logdos_initializer;    // Initialize the logdos values to 1
        }
    }
    
    ~rewl_histogram()
    {
        if (histograms != nullptr){ delete [] histograms; }
    }
};


#endif
