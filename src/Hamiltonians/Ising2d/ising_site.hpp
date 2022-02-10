#ifndef _ISING_SITE_H
#define _ISING_SITE_H
// Use this to define a site for the Ising model.
#include "../Site_Variables/site_variables.hpp"

template<typename data_t>
class Ising_Site : public Site_Variables<data_t>
{
public:
    enum variable_enum 
    {
        spin=0, field, NUM_VARIABLES
    };
    
    void set_variable( const unsigned which_var, const data_t var_val ) { _variables[which_var] = var_val; }
    data_t get_variable( const unsigned which_var ) const { return _variables[which_var]; }

    void set_variable( const variable_enum which_var, const data_t var_val ) { set_variable( static_cast<unsigned>(which_var), var_val ); }
    data_t get_variable( const variable_enum which_var ) const { return get_variable( static_cast<unsigned>(which_var) ); }

    Ising_Site( ) {}
    Ising_Site( const data_t spin_val, const data_t field_val )
    {
        set_variable(variable_enum::spin,  spin_val);
        set_variable(variable_enum::field, field_val);
    }

    virtual ~Ising_Site(){}
private:
    data_t _variables [ NUM_VARIABLES ] = { 0., 0. };
};

#endif /* _ISING_SITE_H */