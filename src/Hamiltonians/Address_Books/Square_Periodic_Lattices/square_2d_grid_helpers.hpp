#ifndef SQUARE_2D_GRID_HELPERS_H
#define SQUARE_2D_GRID_HELPERS_H

#include <cstdint>
#include <helpful_global_macros.hpp>   // included from util

/* All functions are based on the flattened mapping:
   
   site = y_index * Lx + x_index
*/

// Modify these macros only so the inline/template constexpr 
// functions are always identical!!!
#define SITE_FUNCTION(x_index, y_index, Lx) (y_index * Lx + x_index)
#define SITE_X_FUNCTION( site, Lx ) (site % Lx)
#define SITE_Y_FUNCTION( site, Lx ) (site / Lx)

// Functions for use when Lx is not known at compile time
inline std::uint32_t site_from_indices( std::uint32_t x_index, std::uint32_t y_index, std::uint32_t Lx ) { return SITE_FUNCTION(x_index, y_index, Lx); }
inline std::uint32_t site_x_index( std::uint32_t site, std::uint32_t Lx ) { return SITE_X_FUNCTION(site, Lx); }
inline std::uint32_t site_y_index( std::uint32_t site, std::uint32_t Lx ) { return SITE_Y_FUNCTION(site, Lx); }


// Functions for use when Lx is known at compile time
template<std::uint32_t Lx>
constexpr std::uint32_t site_from_indices( std::uint32_t x_index, std::uint32_t y_index ) { return SITE_FUNCTION(x_index, y_index, Lx); }

template<std::uint32_t Lx>
constexpr std::uint32_t site_x_index( std::uint32_t site ) { return SITE_X_FUNCTION(site, Lx); }

template<std::uint32_t Lx>
constexpr std::uint32_t site_y_index( std::uint32_t site ) { return SITE_Y_FUNCTION(site, Lx); }

// PBC functions for a 2D square grid when Lx and Ly are known at compile time
// The order is what I think is the most cache friendly with this site decomposition.

// Zero neighbor (x, y - 1)
template<std::uint32_t Lx, std::uint32_t Ly = Lx>
constexpr std::uint32_t rect_pbc_2d_neighbor_0( const std::uint32_t site_x, const std::uint32_t site_y )
{
   return site_from_indices<Lx>( site_x, _BRANCHLESS_TERNARY( site_y == 0, Ly - 1, site_y - 1 ) );
}

// One neighbor (x - 1, y)
template<std::uint32_t Lx, std::uint32_t Ly = Lx>
constexpr std::uint32_t rect_pbc_2d_neighbor_1( const std::uint32_t site_x, const std::uint32_t site_y )
{
   return site_from_indices<Lx>( _BRANCHLESS_TERNARY( site_x == 0, Lx - 1, site_x - 1 ), site_y );
}

// Two neighbor (x + 1, y)
template<std::uint32_t Lx, std::uint32_t Ly = Lx>
constexpr std::uint32_t rect_pbc_2d_neighbor_2( const std::uint32_t site_x, const std::uint32_t site_y )
{
   return site_from_indices<Lx>( _BRANCHLESS_TERNARY( site_x == Lx - 1, 0, site_x + 1 ), site_y );
}

// Three neighbor (x, y + 1)
template<std::uint32_t Lx, std::uint32_t Ly = Lx>
constexpr std::uint32_t rect_pbc_2d_neighbor_3( const std::uint32_t site_x, const std::uint32_t site_y )
{
   return site_from_indices<Lx>( site_x, _BRANCHLESS_TERNARY( site_y == Ly - 1, 0, site_y + 1 ) );
}

#undef SITE_FUNCTION
#undef SITE_X_FUNCTION
#undef SITE_Y_FUNCTION

#endif /* SQUARE_2D_GRID_HELPERS_H */