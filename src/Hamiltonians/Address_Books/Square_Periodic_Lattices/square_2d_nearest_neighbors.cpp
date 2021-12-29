// Implementation file for the square (periodic) 2D
// lattice with nearest neighbors only.

#include "square_2d_nearest_neighbors.hpp"
#include <helpful_global_macros.hpp>   // included from util
#include <iostream>

Square_2D_Nearest_Neighbors::Square_2D_Nearest_Neighbors( const SiteType Lx, const SiteType Ly ) : _Lx(Lx), _Ly(Ly)
{
    /* intentionally empty */
}

void Square_2D_Nearest_Neighbors::assign_neighbors( SiteType site ) const
{
    SiteType site_x = site_x_index( site, _Lx );
    SiteType site_y = site_y_index( site, _Lx );
    SiteType neighbor_x = 0, neighbor_y = 0;

    // Order of the neighbors is (x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)
    // Neighbor 0:
    neighbor_x = _BRANCHLESS_TERNARY( site_x == 0, _Lx - 1, site_x - 1 );
    neighbor_y = site_y;
    _site_neighbors[site][0] = site_from_indices( neighbor_x, neighbor_y, _Lx );
    
    // Neighbor 1:
    neighbor_x = _BRANCHLESS_TERNARY( site_x == _Lx - 1, 0, site_x + 1 );
    neighbor_y = site_y;
    _site_neighbors[site][1] = site_from_indices( neighbor_x, neighbor_y, _Lx );
    
    // Neighbor 2:
    neighbor_x = site_x;
    neighbor_y = _BRANCHLESS_TERNARY( site_y == 0, _Ly - 1, site_y - 1 );
    _site_neighbors[site][2] = site_from_indices( neighbor_x, neighbor_y, _Lx );
    
    // Neighbor 3:
    neighbor_x = site_x;
    neighbor_y = _BRANCHLESS_TERNARY( site_y == _Ly - 1, 0, site_y + 1 );
    _site_neighbors[site][3] = site_from_indices( neighbor_x, neighbor_y, _Lx );

    return;
}

void Square_2D_Nearest_Neighbors::write_address_book() const
{
    for (SiteType site = 0; site != _Nxy; ++site)
        assign_neighbors( site );
}

void Square_2D_Nearest_Neighbors::initialize()
{
    _site_neighbors = new NeighborList [ _Nxy ];
    write_address_book();
}

void Square_2D_Nearest_Neighbors::print() const
{
    std::cout << "\nNearest Neighbor Lists for a " << _Lx << " x " << _Ly << " periodic square grid\n";
    for ( SiteType site = 0; site != _Nxy; ++site )
    {
        std::cout << "\nSite " << site << ": ";
        for (NeighborIterator itr = neighbor_begin(site); itr != neighbor_end(site); ++itr)
            std::cout << *itr << " ";
        std::cout << "\n";
    }
    std::cout << "\n";
}

Lattice_Address_Book::SiteType Square_2D_Nearest_Neighbors::neighbor( const SiteType site, const SiteType neighbor ) const 
{
    return _site_neighbors[site][neighbor];
}

Square_2D_Nearest_Neighbors::NeighborIterator Square_2D_Nearest_Neighbors::neighbor_begin( const SiteType site ) const 
{
    return _site_neighbors[site].begin();
}

Square_2D_Nearest_Neighbors::NeighborIterator Square_2D_Nearest_Neighbors::neighbor_end( const SiteType site ) const 
{
    return _site_neighbors[site].end();
}
