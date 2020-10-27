#ifndef GLAZIER
#define GLAZIER

/* This header contains templated glazier 
 * functionality -- that is, window making
 * functions. */

#include <cmath>
#include <stdio.h>

/* ****************************************************** */
static constexpr int single_bin_overlap = -1;
/* ****************************************************** */

template<typename data_t>
struct window_data 
{
    // Relevant data for the energy windows
    // defined on the interval: [minimum, maximum)
    // with num_bins bins of size bin_size.
    data_t minimum = 0.;
    data_t maximum = 0.;
    data_t bin_size = 0.;
    size_t num_bins = 0;
};

template<typename data_t, class histogram_index_functor>
struct glazier 
{
    const data_t global_min;
    const data_t global_max;
    const data_t global_bin_size;     // The global bin size MUST be the same across ALL windows
    const size_t num_windows;
    const size_t replicas_per_window;
    const data_t window_overlap;

    window_data<data_t> * all_windows = nullptr;

    glazier(const data_t gmin, const data_t gmax, 
            const data_t gbs, const size_t nw,
            const size_t rpw, const data_t ow) : global_min(gmin), global_max(gmax),
                                                 global_bin_size(gbs), num_windows(nw), 
                                                 replicas_per_window(rpw), window_overlap(ow)
    {
        all_windows = new window_data<data_t> [ num_windows *  replicas_per_window ];
    }
    
    ~glazier()
    {
        delete [] all_windows; 
    }
    
    // TODO: This by default will only increase the
    // size of windows from smallest to largest. It
    // would be better to include a "symmetric" 
    // version which increases in size and then 
    // decreases.
    void construct_windows(); 
};

template<typename data_t, class histogram_index_functor>
void glazier<data_t, histogram_index_functor>::construct_windows()
{
    // Divide global windows into num_windows windows.
    // This measurement will be used to set up the windows initially.
    const data_t initial_window_size = (global_max - global_min) / static_cast<data_t> (num_windows);

    data_t window_min = global_min;
    data_t window_max = global_min + initial_window_size;
    size_t window_bins = static_cast<size_t> ( (window_max - window_min) / global_bin_size );
         
    // Set up the lowest energy window
    for ( size_t replica = 0; replica != replicas_per_window; ++replica )
    {
        all_windows[ replica ].minimum = window_min;
        all_windows[ replica ].maximum = window_max;
        all_windows[ replica ].bin_size = global_bin_size;
        all_windows[ replica ].num_bins = window_bins;
    }
    
    // Set up all other bins
    for ( size_t wdx = 1; wdx != num_windows; ++wdx )
    {
        window_min = global_min + initial_window_size * static_cast<data_t> (wdx);
        window_max = global_min + initial_window_size * static_cast<data_t> (wdx + 1);

        histogram_index_functor indexer (all_windows[ (wdx - 1) * replicas_per_window ].minimum,
                                         all_windows[ (wdx - 1) * replicas_per_window ].maximum,
                                         all_windows[ (wdx - 1) * replicas_per_window ].bin_size );
        
        if ( window_overlap == static_cast<data_t>(single_bin_overlap) )
        {
            window_min -= global_bin_size;
        }
        else
        {
            // Find the nearest bin to grab onto
            window_min -= window_overlap * ( all_windows[ (wdx - 1) * replicas_per_window ].maximum - all_windows[ (wdx - 1) * replicas_per_window ].minimum );
            // Get the index from the previous window and then scale
            // adjust appropriately
            window_min = all_windows[ (wdx - 1) * replicas_per_window ].minimum + global_bin_size * static_cast<data_t> (indexer( window_min )); 
        }

        window_bins = static_cast<size_t> ( (window_max - window_min) / global_bin_size );

        for ( size_t replica = 0; replica != replicas_per_window; ++replica )
        {
            all_windows[ wdx * replicas_per_window + replica ].minimum = window_min;
            all_windows[ wdx * replicas_per_window + replica ].maximum = window_max;
            all_windows[ wdx * replicas_per_window + replica ].bin_size = global_bin_size;
            all_windows[ wdx * replicas_per_window + replica ].num_bins = window_bins;
        }
    }
}

#endif 
