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
    using SiteIndexType = _ISING_INDEX_TYPE;
    Ising_Site_Container_Base();

    virtual const data_t & get_field_value( const SiteIndexType site ) const;
    virtual const data_t & get_spin_value( const SiteIndexType site ) const;
    virtual const data_t get_field_value( const SiteIndexType site ) const;
    virtual const data_t get_spin_value( const SiteIndexType site ) const;
    virtual void set_spin_value( const SiteIndexType site, const data_t spin_val ) const;
    virtual void set_field_value( const SiteIndexType site, const data_t field_val ) const;

    virtual ~Ising_Site_Container_Base();
};

// The general Ising site container
template<typename data_t, _ISING_INDEX_TYPE Nsites, bool constant_field>
class Ising_Site_Container<data_t, Nsites, constant_field> : public Ising_Site_Container_Base<data_t, Nsites> 
{
public:
    void import_fields( const Ising_Site_Container<data_t, Nsites, constant_field> & other_container ) const
    {
        for ( SiteIndexType site = 0; site != Nsites; ++site )
            set_field_value( site, other_container.get_field_value( site ) );
    }

    void import_spins( const Ising_Site_Container<data_t, Nsites, constant_field> & other_container ) const
    {
        for ( SiteIndexType site = 0; site != Nsites; ++site )
            set_field_value( site, other_container.get_spin_value( site ) );
    }
};

// Spatially-varying field container
template<typename data_t, _ISING_INDEX_TYPE Nsites>
class Ising_Site_Container<data_t, Nsites, false>
{
public:
    Ising_Site_Container(){}

    const data_t & get_spin_value( const SiteIndexType site )  const override { return _all_sites[site][Ising_Site::variable_enum::spin]; }
    const data_t & get_field_value( const SiteIndexType site ) const override { return _all_sites[site][Ising_Site::variable_enum::field]; }
    const data_t get_spin_value( const SiteIndexType site )  const override { return _all_sites[site][Ising_Site::variable_enum::spin]; }
    const data_t get_field_value( const SiteIndexType site ) const override { return _all_sites[site][Ising_Site::variable_enum::field]; }
    void set_spin_value(  const SiteIndexType site, const data_t spin_val ) const override {  _all_sites[site].set_variable( Ising_Site::variable_enum::spin, spin_val ); }
    void set_field_value( const SiteIndexType site, const data_t field_val ) const override { _all_sites[site].set_variable( Ising_Site::variable_enum::field, field_val ); }
    
    void import_fields( const data_t * const contiguous_field_array ) const
    {
        for ( SiteIndexType site = 0; site != Nsites; ++site )
            set_field_value( site, contiguous_field_array[site] );
    }

    ~Ising_Site_Container(){}
 
private:
    std::array<Ising_Site<data_t>, Nsites> _all_sites;
};

// Spatially-constant field container
template<typename data_t, _ISING_INDEX_TYPE Nsites>
class Ising_Site_Container<data_t, Nsites, true>
{
public:
    Ising_Site_Container(){}

    const data_t & get_spin_value( const SiteIndexType site )  const override { return _all_sites[site]; }
    const data_t & get_field_value( const SiteIndexType site ) const override { return _field_value; }
    const data_t get_spin_value( const SiteIndexType site )  const override { return _all_sites[site]; }
    const data_t get_field_value( const SiteIndexType site ) const override { return _field_value; }
    void set_spin_value(  const SiteIndexType site, const data_t spin_val ) const override {  _all_sites[site] = spin_val; }
    void set_field_value( const SiteIndexType site, const data_t field_val ) const override { _field_value = field_val; }

    ~Ising_Site_Container(){}
 
private:
    data_t _field_value;
    std::array<data_t, Nsites> _all_sites;
};

#undef _ISING_INDEX_TYPE

#endif /* _ISING_SITE_CONTAINER_H */