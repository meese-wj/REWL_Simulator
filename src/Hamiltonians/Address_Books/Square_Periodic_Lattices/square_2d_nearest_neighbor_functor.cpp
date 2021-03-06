#ifndef SQUARE_2D_NEAREST_NEIGHBOR_FUNCTOR_S
#define SQUARE_2D_NEAREST_NEIGHBOR_FUNCTOR_S

// Implementation file for the square (periodic) 2D
// lattice with nearest neighbors only. 
// 
// THIS SOURCE FILE NEEDS TO BE INCLUDED TO AVOID
// LINKER ERRORS!!!!

#include "square_2d_nearest_neighbor_functor.hpp"
#include "square_2d_grid_helpers.hpp"
#include <helpful_global_macros.hpp>   // included from util
#include <iostream>

template<Lattice_Address_Book::SiteType Lx, Lattice_Address_Book::SiteType Ly>
Square_2D_Nearest_Neighbor_Functor<Lx, Ly>::Square_2D_Nearest_Neighbor_Functor( )
{
    /* intentionally empty */
}

template<Lattice_Address_Book::SiteType Lx, Lattice_Address_Book::SiteType Ly>
void Square_2D_Nearest_Neighbor_Functor<Lx, Ly>::assign_neighbors( SiteType site ) const
{
    SiteType site_x = site_x_index<Lx>( site );
    SiteType site_y = site_y_index<Lx>( site );

    // Order of the neighbors is can be found in "square_2d_grid_helpers.hpp"
    // Neighbor 0:
    _neighbors[0] = rect_pbc_2d_neighbor_0<Lx>(site_x, site_y);
    
    // Neighbor 1:
    _neighbors[1] = rect_pbc_2d_neighbor_1<Lx>(site_x, site_y);
    
    // Neighbor 2:
    _neighbors[2] = rect_pbc_2d_neighbor_2<Lx>(site_x, site_y);
    
    // Neighbor 3:
    _neighbors[3] = rect_pbc_2d_neighbor_3<Lx>(site_x, site_y);

    update_current_site(site);
    return;
}

template<Lattice_Address_Book::SiteType Lx, Lattice_Address_Book::SiteType Ly>
void Square_2D_Nearest_Neighbor_Functor<Lx, Ly>::write_address_book() const
{
    /* intentionally empty */
}

template<Lattice_Address_Book::SiteType Lx, Lattice_Address_Book::SiteType Ly>
void Square_2D_Nearest_Neighbor_Functor<Lx, Ly>::initialize()
{
    update_current_site( 0 );
    assign_neighbors( _current_site ); // Just a dummy assignment to clear out uninitialized values   
}

template<Lattice_Address_Book::SiteType Lx, Lattice_Address_Book::SiteType Ly>
void Square_2D_Nearest_Neighbor_Functor<Lx, Ly>::print() const
{
    std::cout << "\nNearest Neighbor Lists for a " << Lx << " x " << Ly << " periodic square grid\n";
    for ( SiteType site = 0; site != total_sites(); ++site )
    {
        std::cout << "\nSite " << site << ": ";
        for (NeighborIterator itr = neighbor_begin(site); itr != neighbor_end(site); ++itr)
            std::cout << *itr << " ";
        std::cout << "\n";
    }
    std::cout << "\n";
}

template<Lattice_Address_Book::SiteType Lx, Lattice_Address_Book::SiteType Ly>
Lattice_Address_Book::SiteType Square_2D_Nearest_Neighbor_Functor<Lx, Ly>::current_site() const 
{
    return _current_site;
}

template<Lattice_Address_Book::SiteType Lx, Lattice_Address_Book::SiteType Ly>
Lattice_Address_Book::SiteType Square_2D_Nearest_Neighbor_Functor<Lx, Ly>::total_sites() const 
{
    return Lx * Ly;
}

template<Lattice_Address_Book::SiteType Lx, Lattice_Address_Book::SiteType Ly>
void Square_2D_Nearest_Neighbor_Functor<Lx, Ly>::update_current_site( const SiteType site ) const 
{
    _current_site = site;
}

template<Lattice_Address_Book::SiteType Lx, Lattice_Address_Book::SiteType Ly>
Lattice_Address_Book::SiteType Square_2D_Nearest_Neighbor_Functor<Lx, Ly>::neighbor( const SiteType site, const SiteType neighbor ) const 
{
    if (site != current_site())
        assign_neighbors( site );
    return _neighbors[neighbor];
}

template<Lattice_Address_Book::SiteType Lx, Lattice_Address_Book::SiteType Ly>
typename Square_2D_Nearest_Neighbor_Functor<Lx, Ly>::NeighborIterator Square_2D_Nearest_Neighbor_Functor<Lx, Ly>::neighbor_begin( const SiteType site ) const
{
    assign_neighbors( site );
    return _neighbors;
}

template<Lattice_Address_Book::SiteType Lx, Lattice_Address_Book::SiteType Ly>
typename Square_2D_Nearest_Neighbor_Functor<Lx, Ly>::NeighborIterator Square_2D_Nearest_Neighbor_Functor<Lx, Ly>::neighbor_end( const SiteType site ) const
{
    return _neighbors + _num_neighbors;
}

#endif /* SQUARE_2D_NEAREST_NEIGHBOR_FUNCTOR_S */