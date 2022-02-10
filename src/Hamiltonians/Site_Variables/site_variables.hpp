#ifndef _SITE_VARIABLES_H
#define _SITE_VARIABLES_H
// This is an abstract class meant to be used 
// by Hamiltonians down the line. The idea is 
// to optimize spatial locality by putting all
// variables relevant for a Lattice site 
// together. This is particularly relevant for
// random field models.

template<typename data_t>
class Site_Variables
{
public:
    Site_Variables() = default;
    virtual ~Site_Variables(){}
    virtual data_t get_variable( unsigned which ) const = 0;
};

#endif /* _SITE_VARIABLES_H */

