#ifndef SQUARE_2D_GRID_HELPERS_H
#define SQUARE_2D_GRID_HELPERS_H

#include <cstdint>

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

#undef SITE_FUNCTION
#undef SITE_X_FUNCTION
#undef SITE_Y_FUNCTION

#endif /* SQUARE_2D_GRID_HELPERS_H */