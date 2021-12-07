#ifndef REWL_HISTOGRAMS
#define REWL_HISTOGRAMS

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
#include <cfloat>
#include <iostream>
#include <string>

#if PRINT_HISTOGRAM
#include <string>
#include <filesystem>
#include <fstream>
#endif

/* Define some constexpr to be used later on. */
static constexpr uint64_t count_initializer = 0;
static constexpr uint64_t logdos_initializer = 0;
constexpr float ratio_failure = FLT_MAX;
#if PRINT_HISTOGRAM
const std::string histogram_file = "histogram_";
#endif
/* ****************************************** */

// Use these pairs for faster memory access
// later on in the simulation.
template<typename data_t>
struct value_pair
{
    uint64_t count;
    data_t logdos; 
};

// Define the templated rewl_histogram
// struct to be used in the REWL simulation.
template<typename data_t>
struct rewl_histograms
{
    const data_t min_value;        // Greatest lower bound of the REWL histogram (included).
    const data_t max_value;        // Least upper bound of the REWL histogram (not included).
    const data_t bin_size;         // Size of the bins for histograms. (This must be the same
                                   // across all walkers!).
    
    const size_t num_bins;

    // Declare the arrays.
    value_pair<data_t> * histograms = nullptr;

    rewl_histograms(const data_t min, const data_t max, const data_t b_size, const size_t nbins) : min_value(min), max_value(max), bin_size(b_size), num_bins(nbins)
    {
        histograms = new value_pair<data_t> [ num_bins ];
        for ( size_t idx = 0; idx != num_bins; ++idx )
        {
            histograms[ idx ].count = count_initializer;                               // Initialize the counts to zero
            histograms[ idx ].logdos = static_cast<data_t> (logdos_initializer);       // Initialize the logdos values to 1
        }
    }
    
    ~rewl_histograms()
    {
        delete [] histograms;
    }

    // Write the getter/setter for the histogram counts
    uint64_t get_count(const size_t bin) const { return histograms[ bin ].count; }
    void increment_count(const size_t bin){ ++( histograms[ bin ].count ); }
   
    // Write the getter/setter for the histogram logdos values
    data_t get_logdos( const size_t bin ) const { return histograms[ bin ].logdos; }
    void increment_logdos(const size_t bin, const data_t incrementer){ histograms[ bin ].logdos += incrementer; }

    // Reset the counts
    void reset_counts() const
    {
#if REDUCE_LOGDOS
        data_t logdos_reducer = -histograms[ 0 ].logdos + logdos_initializer;
#endif
        for ( size_t idx = 0; idx != num_bins; ++idx )
        {
            histograms[ idx ].count = count_initializer;
#if REDUCE_LOGDOS
            histograms[ idx ].logdos += logdos_reducer;
#endif
        }
    }

    // Return flatness in the counts
    float count_flatness() const;

    // Export out the density of states as a deep copy
    void export_logdos( data_t *& data_array ) const
    {
        data_array = new data_t [ num_bins ];
        for ( size_t bin = 0; bin != num_bins; ++bin )
            data_array[ bin ] = histograms[ bin ].logdos;
    }

#if PRINT_HISTOGRAM
    // Print counts for the histogram. This function overwrites
    // per each iteration.
    void print_histogram_counts( const size_t iteration, const std::filesystem::path & histogram_path ) const;
#endif
};

#if ONE_OVER_T_ALGORITHM
template<typename data_t>
float rewl_histograms<data_t>::count_flatness() const
{
    // This algorithm is based on that given in
    // DOI: 10.1103/PhysRevE.75.046701
    // Here "flatness" corresponds to every state
    // being sampled at least once.
    uint64_t sum_bins = ( get_count(0) > 0 );
    for ( size_t idx = 1; idx != num_bins; ++idx )
        sum_bins *= ( get_count(idx) > 0 );

    // If every bin has been sampled, then 
    // sum_bins == 1 otherwise sum_bins == 0.
    return sum_bins;
}
#else
template<typename data_t>
float rewl_histograms<data_t>::count_flatness() const
{
    uint64_t minimum = get_count(0);
    uint64_t maximum = minimum;    

    for ( size_t idx = 1; idx != num_bins; ++idx )
    {
        if ( get_count(idx) < minimum )
            minimum = get_count(idx);
        if ( get_count(idx) > maximum )
            maximum = get_count(idx);
    }

    if ( minimum != 0 )
        return (static_cast<float> (maximum - minimum) / static_cast<float> (minimum));

    // If the minimum is zero, then return a 
    // silly number to check against.
    return ratio_failure;
}
#endif

#if PRINT_HISTOGRAM
template<typename data_t>
void rewl_histograms<data_t>::print_histogram_counts( const size_t iteration, const std::filesystem::path & histogram_path ) const
{
    std::ofstream printer;
    std::filesystem::path file_name = histogram_path / ( histogram_file + std::to_string(iteration) + ".txt" );
    printer.open(file_name);
    for( size_t bin = 0; bin != num_bins; ++bin )
    {
        printer << min_value + bin_size * bin << "    " << histograms[bin].count << "    " << histograms[bin].logdos;
        if ( bin != num_bins - 1 )
            printer << "\n";
    }
    printer.close();
}
#endif

#endif
