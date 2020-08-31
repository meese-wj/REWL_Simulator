#ifndef HISTOGRAM_INDEX

#include <cmath>

template<typename data_t>
size_t histogram_index(const data_t value, const data_t min_value, const data_t bin_size)
{
    return static_cast<size_t> ( (value - min_value) / bin_size );
}

#endif
