#ifndef PBC_2D_VECTORS
#define PBC_2D_VECTORS
/* This file defines 2d spatial
 * vector operations on a periodic 
 * lattice. 
 * 
 * There are notable differences 
 * between the nature of vector 
 * operations between the PBC
 * case and the non PBC case.
 * */

#include <cstdint>
#include <cmath>

// Define a datatype that will 
// be useful for the PBC functions.
template<typename T>
struct pbc_2d_vector
{
    T x = 0;
    T y = 0;

    pbc_2d_vector(){}
    pbc_2d_vector( T _x, T _y ) : x(_x), y(_y) {}
    pbc_2d_vector( const pbc_2d_vector<T> & a ) : x(a.x), y(a.y) {}

    ~pbc_2d_vector(){}
};


// Add ax1 and ax2 
// in a way that obeys 
// the periodic boundary 
// conditions. Their sum
// cannot exceed L_axis.
template<typename T>
T pbc_add_1d( const std::uint32_t L_axis, const T ax1, const T ax2 )
{
    T result = ax1 + ax2;

    while ( result >= static_cast<T>(L_axis) )
        result -= L_axis;
    while ( result < static_cast<T>(L_axis) )
        result += L_axis;

    return result;
}

// Subtract ax2 from ax1: ax1 - ax2
// The distance between two points
// can be AT MOST L_axis / 2. 
template<typename T>
T pbc_subtract_1d( const std::uint32_t L_axis, const T ax1, const T ax2 )
{
    T result = ax1 - ax2;

    while ( result >= static_cast<T>(L_axis) / 2. )
        result -= L_axis;
    while ( result < static_cast<T>(L_axis) / 2. )
        result += L_axis;

    return result;
}

// Add vectors A and B
template<typename T>
pbc_2d_vector<T> pbc_add( const std::uint32_t Lx, const std::uint32_t Ly, const pbc_2d_vector<T> & A, const pbc_2d_vector<T> & B )
{
    return pbc_2d_vector<T>( pbc_add_1d(Lx, A.x, B.x), pbc_add_1d(Ly, A.y, B.y) );
}

// Subtract B from A: A - B
template<typename T>
pbc_2d_vector<T> pbc_subtract( const std::uint32_t Lx, const std::uint32_t Ly, const pbc_2d_vector<T> & A, const pbc_2d_vector<T> & B )
{
    return pbc_2d_vector<T>( pbc_subtract_1d(Lx, A.x, B.x), pbc_subtract_1d(Ly, A.y, B.y) );
}

// Compute the magnitude^2 of the 
// given vector and possibly 
// cast it into a different type.
template<typename O, typename T>
O square_magnitude( const pbc_2d_vector<T> & A )
{
    return static_cast<O>( A.x * A.x + A.y * A.y );
}

// Compute the magnitude of the 
// given vector and possibly 
// cast it into a different type.
template<typename O, typename T>
O magnitude( const pbc_2d_vector<T> & A )
{
    return sqrt( square_magnitude<O, T>( A ) );
}

// Calculate the distance 
// between two vectors
template<typename O, typename T>
O pbc_distance( const std::uint32_t Lx, const std::uint32_t Ly, const pbc_2d_vector<T> & A, const pbc_2d_vector<T> & B )
{
    T resultant = pbc_subtract( Lx, Ly, A, B );
    return magnitude<O, T>( resultant );
}

// When computing the dot product
// between two vectors we need to 
// measure both relative to some 
// origin.
// This should be the last function
// that needs to be modified relative
// to the non-periodic case.
template<typename O, typename T>
O pbc_scalar_product( const std::uint32_t Lx, const std::uint32_t Ly, 
                      const pbc_2d_vector<T> & A, const pbc_2d_vector<T> & B,
                      const pbc_2d_vector<T> origin = pbc_2d_vector<T> (0, 0) )
{
    const pbc_2d_vector<T> temp_A = pbc_subtract( Lx, Ly, A, origin );
    const pbc_2d_vector<T> temp_B = pbc_subtract( Lx, Ly, B, origin );
    return static_cast<O>( temp_A.x * temp_B.x + temp_A.y * temp_B.y );
}

#endif
