#ifndef ARRAY_SHIFT
/* This function shifts an array by
 * a common scalar value. */
#include <cmath>

template<typename data_t>
void array_shift_by_value( const data_t value, const size_t size, data_t * const data_array )
{
   for ( size_t idx = 0; idx != size; ++idx )
       data_array[ idx ] += value;
}


#endif
