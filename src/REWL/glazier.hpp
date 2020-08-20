#ifndef GLAZIER

/* This header contains templated glazier 
 * functionality -- that is, window making
 * functions. */

#include <cmath>

template<typename data_t>
struct window_data 
{
    // Relevant data for the energy windows
    // defined on the interval: [minimum, maximum)
    // with num_bins bins of size bin_size.
    data_t minimum;
    data_t maximum;
    data_t bin_size;
    size_t num_bins;
};

template<typename data_t>
struct glazier 
{
    const data_t global_min;
    const data_t global_max;
    const size_t num_windows;
    const size_t replicas_per_window;

    window_data<data_t> * all_windows = nullptr;

    glazier(const data_t gmin, const data_t gmax, 
            const size_t nw, const size_t rpw) : global_min(gmin), global_max(gmax),
                                                 num_windows(nw), replicas_per_window(rpw)
    {}
    
    ~glazier()
    {
        if (all_windows != nullptr){ delete [] all_windows; }
    }
        
};

#endif 
