#ifndef RANDOM_FIELDS
#define RANDOM_FIELDS
/* This file will build the random fields
 * to be used in the simulations.      */

#include <cmath>
#include <random_number_generators.hpp>
#include <chrono>
#include <stdio.h>

enum disorder_distribution
{
    uniform
};

// Produce a random field value on the 
// interval [-strength/2, strength/2] for a
// uniform random number rng() in [0,1].
template<typename data_t>
data_t uniform_random_field( random_number_generator<data_t> & rng, const data_t strength )
{
    return strength * ( -0.5 + rng() );
}

template<typename data_t>
void generate_random_field( const size_t num_sites, const data_t strength,
                            data_t *& field_array, const int type )
{
    delete [] field_array;
    field_array = new data_t [num_sites]();
    std::cout << "\n field pointer in function = " << field_array << "\n";
    // Build a random number generator 
    std::uint64_t seed = static_cast<std::uint64_t>( std::chrono::high_resolution_clock::now().time_since_epoch().count() );
    random_number_generator<data_t> rng (seed);

    // Choose the right random field distribution
    data_t (*field_generator)( random_number_generator<data_t> &, const data_t);
    if ( type == disorder_distribution::uniform )
    {
        field_generator = &uniform_random_field<data_t>;
    }

    // Populate the fields
    for ( size_t idx = 0; idx != num_sites; ++idx )
    {
        data_t h = field_generator(rng, strength);
        printf("h[%ld] = %e\n", idx, h);
        field_array[idx] = h;
    }

}

#endif
