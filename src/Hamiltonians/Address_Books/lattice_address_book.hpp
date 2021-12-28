#ifndef LATTICE_ADDRESS_BOOK_H
#define LATTICE_ADDRESS_BOOK_H

/* 
    This header defines a base class that unites lattice
    points through some geometrical mapping. It should be
    used essentially as an iterator that defines which 
    lattice points at a particular site index are nearby 
    others in space. Additionally, it should define an 
    iteration scheme for each site to move between its 
    neighbors.
*/

#include <cstdint>

class Lattice_Address_Book
{
public:
    using SiteType = std::uint32_t;
    using iterator = SiteType*;

    Lattice_Address_Book(/* args */) = default;
    virtual void assign_neighbors( const SiteType site ) const = 0;   // Makes this class abstract. 
    virtual void write_address_book() const = 0;
    virtual void initialize() = 0;
    virtual ~Lattice_Address_Book() {}
};


#endif /* LATTICE_ADDRESS_BOOK_H */