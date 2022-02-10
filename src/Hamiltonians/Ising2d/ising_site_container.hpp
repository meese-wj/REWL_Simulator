#ifndef _ISING_SITE_CONTAINER_H
#define _ISING_SITE_CONTAINER_H
// Use this object to house a AoS
// data structure of the Ising Sites.

#include <cstdint>
#include <array>
#include "ising_site.hpp"

#define _ISING_INDEX_TYPE std::uint32_t

// Define a Container Base class with all common
// functions. Then specialize the members with 
// a set of derived classes.
template<typename data_t, _ISING_INDEX_TYPE Nsites>
class Ising_Site_Container_Base
{
public:
    Ising_Site_Container_Base() = default;

    virtual const data_t get_field_value( const _ISING_INDEX_TYPE site ) const;
    virtual const data_t get_spin_value( const _ISING_INDEX_TYPE site ) const;
    virtual void set_spin_value( const _ISING_INDEX_TYPE site, const data_t spin_val );
    virtual void set_field_value( const _ISING_INDEX_TYPE site, const data_t field_val );

    virtual ~Ising_Site_Container_Base(){};

    void export_contiguous_spins( data_t *& contiguous_export_array ) const
    {
        delete [] contiguous_export_array;
        contiguous_export_array = new data_t [Nsites];
        for ( _ISING_INDEX_TYPE site = 0; site != Nsites; ++site )
            contiguous_export_array[site] = this -> get_spin_value(site);
    }

    void import_fields( const Ising_Site_Container_Base<data_t, Nsites> & other_container )
    {
        for ( _ISING_INDEX_TYPE site = 0; site != Nsites; ++site )
            this -> set_field_value( site, other_container.get_field_value( site ) );
    }

    void import_contiguous_spins( const data_t * const contiguous_import_array )
    {
        for ( _ISING_INDEX_TYPE site = 0; site != Nsites; ++site )
            this -> set_spin_value( site, contiguous_import_array[site] );
    }

    void import_spins( const Ising_Site_Container_Base<data_t, Nsites> & other_container )
    {
        for ( _ISING_INDEX_TYPE site = 0; site != Nsites; ++site )
            this -> set_field_value( site, other_container.get_spin_value( site ) );
    }
};

// The general Ising site container
template<typename data_t, _ISING_INDEX_TYPE Nsites, bool constant_field>
class Ising_Site_Container : public Ising_Site_Container_Base<data_t, Nsites> 
{
public:
    
};

// Spatially-varying field container
template<typename data_t, _ISING_INDEX_TYPE Nsites>
class Ising_Site_Container<data_t, Nsites, false> : public Ising_Site_Container_Base<data_t, Nsites>
{
public:
    Ising_Site_Container( const data_t const_field_value ) : _const_field_value(const_field_value){}

    const data_t get_spin_value( const _ISING_INDEX_TYPE site )  const { return _all_sites[site][Ising_Site<data_t>::variable_enum::spin]; }
    const data_t get_field_value( const _ISING_INDEX_TYPE site ) const { return _all_sites[site][Ising_Site<data_t>::variable_enum::field]; }
    void set_spin_value(  const _ISING_INDEX_TYPE site, const data_t spin_val ) {  _all_sites[site].set_variable( Ising_Site<data_t>::variable_enum::spin, spin_val ); }
    void set_field_value( const _ISING_INDEX_TYPE site, const data_t field_val ) { _all_sites[site].set_variable( Ising_Site<data_t>::variable_enum::field, field_val ); }
    
    void import_fields( const data_t * const contiguous_field_array ) const
    {
        for ( _ISING_INDEX_TYPE site = 0; site != Nsites; ++site )
            set_field_value( site, contiguous_field_array[site] );
    }

    ~Ising_Site_Container(){}
 
private:
    const data_t _const_field_value;
    std::array<Ising_Site<data_t>, Nsites> _all_sites;
};

// Spatially-constant field container
template<typename data_t, _ISING_INDEX_TYPE Nsites>
class Ising_Site_Container<data_t, Nsites, true> : public Ising_Site_Container_Base<data_t, Nsites>   
{
public:
    Ising_Site_Container( const data_t const_field_value ) : _field_value(const_field_value){}

    const data_t get_spin_value( const _ISING_INDEX_TYPE site )  const { return _all_sites[site]; }
    const data_t get_field_value( const _ISING_INDEX_TYPE site ) const { return _field_value; }
    void set_spin_value(  const _ISING_INDEX_TYPE site, const data_t spin_val ) {  _all_sites[site] = spin_val; }
    void set_field_value( const _ISING_INDEX_TYPE site, const data_t field_val ) { /* intentionally empty */ }

    ~Ising_Site_Container(){}
 
private:
    const data_t _field_value;
    std::array<data_t, Nsites> _all_sites;
};

#undef _ISING_INDEX_TYPE

#endif /* _ISING_SITE_CONTAINER_H */