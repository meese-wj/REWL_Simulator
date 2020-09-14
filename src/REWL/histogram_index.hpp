#ifndef HISTOGRAM_INDEX

#include <cmath>

// Define a functor to return the index value of
// either histogram given the energy value.
template<typename data_t>
struct histogram_index 
{
    const data_t min_value;
    const data_t max_value;
    const data_t bin_size;
    
    histogram_index( const data_t min, const data_t max, const data_t b_size ) : min_value(min), max_value(max), bin_size(b_size)
    {}

    size_t operator () (const data_t value) const { return static_cast<size_t> ( (value - min_value) / bin_size ); }

    data_t get_bin_min( const size_t bin ) const { return min_value + bin_size * static_cast<data_t>(bin); }
    
    bool energy_too_low( const data_t value ) const { return ( static_cast<int>(value + 0.5) < static_cast<int>(min_value + 0.5) ); }
    bool energy_too_high( const data_t value ) const { return ( static_cast<int>(value + 0.5) >= static_cast<int>(max_value + 0.5) ); }

    bool energy_in_range( const data_t value ) const
    {
        // Add 0.5 and then truncate to get around floating point error
        return ( static_cast<int>(value + 0.5) >= static_cast<int>(min_value + 0.5) && 
                 static_cast<int>(value + 0.5) <  static_cast<int>(max_value + 0.5)   );
    }
};

#endif
