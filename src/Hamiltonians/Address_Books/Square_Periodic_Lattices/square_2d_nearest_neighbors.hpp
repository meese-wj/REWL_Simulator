#ifndef SQUARE_2D_NEAREST_NEIGHBORS_H
#define SQUARE_2D_NEAREST_NEIGHBORS_H

/*
    This header defines the address book for a 
    2D square (periodic) lattice. Each site has
    4 neighbors. 
*/
#include "../lattice_address_book.hpp"
#include "square_2d_grid_helpers.hpp"
#include <array>

static constexpr Lattice_Address_Book::SiteType NUM_NEIGHBORS = 4;

class Square_2D_Nearest_Neighbors : public Lattice_Address_Book
{    
public:
    using NeighborList = std::array<SiteType, NUM_NEIGHBORS>;
    using NeighborIterator = NeighborList::iterator;

    Square_2D_Nearest_Neighbors( const SiteType Lx, const SiteType Ly );
    void assign_neighbors( const SiteType site ) const override;
    void write_address_book() const override; 
    void initialize() override;

    SiteType neighbor( const SiteType site, const SiteType neighbor ) const;
    NeighborIterator neighbor_begin( const SiteType site ) const;
    NeighborIterator neighbor_end( const SiteType site ) const;

    virtual ~Square_2D_Nearest_Neighbors() override 
    {
        delete [] _site_neighbors;
    }

    void print() const;

private:
    const SiteType _Lx;
    const SiteType _Ly;
    const SiteType _Nxy = _Lx * _Ly;
    const SiteType _num_neighbors = NUM_NEIGHBORS;
    NeighborList * _site_neighbors = nullptr;
};


#endif /* SQUARE_2D_NEAREST_NEIGHBORS_H */
