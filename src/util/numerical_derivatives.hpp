#ifndef NUMERICAL_DERIVATIVES
#define NUMERICAL_DERIVATIVES
/* This file takes derivatives of some vector
 * with respect to another.  */

#include <cmath>
#include <vector>

/* *************************************************************************************************** */
/* These functions apply for any arrays                                                                */
/* *************************************************************************************************** */

// O(h) derivative (Eulerian)
template<typename indep_t, typename dep_t>
inline indep_t order_one_forward_deriv( const size_t index, const indep_t * const fxn, const dep_t * const axis )
{
    return ( fxn[index + 1] - fxn[index] ) / static_cast<indep_t>( axis[index + 1] - axis[index] );
}

template<typename indep_t, typename dep_t>
inline indep_t order_one_backward_deriv( const size_t index, const indep_t * const fxn, const dep_t * const axis )
{
    return ( fxn[index] - fxn[index - 1] ) / static_cast<indep_t>( axis[index] - axis[index - 1] );
}

// O(h^2) derivative (central Eulerian)
// This assumes that the axis discretization is uniform
template<typename indep_t, typename dep_t>
inline indep_t order_two_central_deriv( const size_t index, const indep_t * const fxn, const dep_t * const axis )
{
    return ( fxn[index + 1] - fxn[index - 1] ) / static_cast<indep_t>( 2. * ( axis[index] - axis[index - 1] ) );
}

// O(h^4) derivative 
// This assumes that the axis discretization is uniform
template<typename indep_t, typename dep_t>
inline indep_t order_four_central_deriv( const size_t index, const indep_t * const fxn, const dep_t * const axis )
{
    return ( -fxn[index + 2] + 8. * fxn[index + 1] - 8. * fxn[index - 1] + fxn[index - 2] ) / static_cast<indep_t>( 12. * ( axis[index] - axis[index - 1] ) );
}

/* *************************************************************************************************** */
/* Vector functions are below                                                                          */
/* *************************************************************************************************** */

// Compute the derivative of a vector type
template<typename indep_t, typename dep_t>
indep_t derivative_of_vector( const size_t index, const std::vector<indep_t> & fxn, const std::vector<dep_t> & axis )
{
    indep_t deriv = 0.;
    const size_t sze = fxn.size();

    if      ( index == 0 )
    {
        deriv = order_one_forward_deriv( index, fxn.data(), axis.data() );
    }
    else if ( index == sze - 1 )
    {
        deriv = order_one_backward_deriv( index, fxn.data(), axis.data() );
    }
    else if ( index == 1 || index == sze - 2 )
    {
        deriv = order_two_central_deriv( index, fxn.data(), axis.data() );
    }
    else
    {
        deriv = order_four_central_deriv( index, fxn.data(), axis.data() );
    }
    return deriv;
}


#endif
