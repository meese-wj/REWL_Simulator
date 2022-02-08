#ifndef RANDOM_NUMBER_GENERATORS
#define RANDOM_NUMBER_GENERATORS
/* This file is really just a wrapper around the Standard Library 
 * random number generators. This is to separate the random 
 * number generator from the Wang Landau functionality -- this 
 * will be more helpful for the replica exchange part of the 
 * simulation. It can also be used to try different generators
 * for performance. */
#include <cmath>
#include <random>
#include <limits>


template<typename data_t>
struct random_number_generator
{
    const std::uint64_t my_seed;
    std::mt19937 generator;
    // std::default_random_engine generator;
    std::uniform_real_distribution<data_t> u_real;
    // TODO: Should I add a uniform int distribution too? 
    
    random_number_generator(const std::uint64_t _seed) : my_seed(_seed), generator(my_seed), u_real(0., 1.) 
    {}
    ~random_number_generator() {}
   
    data_t operator () () { return u_real(generator); }
    data_t get_real () { return u_real(generator); }
    data_t get_real_modified(const data_t min, const data_t max){ return min + (max - min) * u_real(generator); }

    data_t get_gaussian( const data_t mean, const data_t stdev )
    {
        static constexpr data_t machine_epsilon = std::numeric_limits<data_t>::epsilon();
        static constexpr data_t two_pi = 2.0 * acos(-1.0);

        data_t u_var = get_real();
        data_t v_var = get_real();

        // Make sure the magnitude won't blow up.
        while ( u_var <= machine_epsilon )
        {
            u_var = get_real();
            v_var = get_real();
        }

        // Compute z0 and z1 according to the Box-Muller transform
        data_t magnitude = stdev * sqrt( -2.0 * log(u_var) );
        data_t z0 = magnitude * cos( two_pi * v_var ) + mean;
        data_t z1 = magnitude * sin( two_pi * v_var ) + mean;
       
        // Both z0 and z1 are sampled from independent Gaussians
        // so return either with equal probability.
        bool return_cos = get_real() < 0.5;
        return z0 * return_cos + z1 * (!return_cos); 
    }
    
};

#endif
