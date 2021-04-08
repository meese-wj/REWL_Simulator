#ifndef RANDOM_FIELDS
#define RANDOM_FIELDS
/* This file will build the random fields
 * to be used in the simulations.      */

#include <cmath>
#include <random_number_generators.hpp>
#include <chrono>

// Generate a set of random values between the fields
// at each site to make sure they are truly uncorrelated.
static constexpr size_t WAIT_BETWEEN_FIELD_GENERATION = 20;

enum disorder_distribution
{
    uniform
};

// Produce a random field value on the 
// interval [-strength, strength] for a
// uniform random number rng() in [0,1].
// This means the typical field strength 
// at any site is "strength".
template<typename data_t>
data_t uniform_random_field( random_number_generator<data_t> & rng, const data_t strength )
{
    return strength * ( -1. + 2. * rng() );
}

template<typename data_t>
void generate_random_field( const size_t num_sites, const data_t strength,
                            data_t *& field_array, const int type )
{
    delete [] field_array;
    field_array = new data_t [num_sites]();
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
        // Generate a bunch of wasted fields to guarantee
        // the state of the random number generator is 
        // uncorrelated between sites.
        for ( size_t reject = 0; reject != WAIT_BETWEEN_FIELD_GENERATION; ++reject )
            data_t fake_h = field_generator(rng, strength);

        data_t h = field_generator(rng, strength);
        field_array[idx] = h;
    }
}

#endif
