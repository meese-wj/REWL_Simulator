#ifndef SQUARE_2D_NEAREST_NEIGHBOR_FUNCTOR_H
#define SQUARE_2D_NEAREST_NEIGHBOR_FUNCTOR_H

/*
    This header defines the address book for a 
    2D square (periodic) lattice. Each site has
    4 neighbors. 

    This class is different than the non-functor
    option by the fact that it computes the 
    nearest-neighbors for each site in real-time, 
    rather than storing them. 

    It is expected that for small lattices, this 
    will be slower than a neighbor table. For large
    systems, however, this approach is probably 
    better as it will lead to fewer cache misses.
*/
#include "../lattice_address_book.hpp"
#include <array>

static constexpr Lattice_Address_Book::SiteType NUM_NEIGHBORS = 4;

template<Lattice_Address_Book::SiteType Lx, Lattice_Address_Book::SiteType Ly = Lx>
class Square_2D_Nearest_Neighbor_Functor : public Lattice_Address_Book
{    
public:
    using NeighborList = SiteType*;
    using NeighborIterator = const SiteType*;

    Square_2D_Nearest_Neighbor_Functor();
    void assign_neighbors( const SiteType site ) const override;
    void write_address_book() const override; 
    void initialize() override;

    SiteType current_site( ) const; 
    void update_current_site( const SiteType site ) const; 
    SiteType neighbor( const SiteType site, const SiteType neighbor ) const;
    NeighborIterator neighbor_begin( const SiteType site ) const;
    NeighborIterator neighbor_end( const SiteType site ) const;

    virtual ~Square_2D_Nearest_Neighbor_Functor() override {}

    void print() const;

private:
    const SiteType _Nxy = Lx * Ly;
    const SiteType _num_neighbors = NUM_NEIGHBORS;
    mutable SiteType _current_site;
    mutable SiteType _neighbors [NUM_NEIGHBORS];
};


#endif /* SQUARE_2D_NEAREST_NEIGHBOR_FUNCTOR_H */
