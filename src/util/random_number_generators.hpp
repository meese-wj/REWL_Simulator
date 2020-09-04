#ifndef RANDOM_NUMBER_GENERATORS
/* This file is really just a wrapper around the Standard Library 
 * random number generators. This is to separate the random 
 * number generator from the Wang Landau functionality -- this 
 * will be more helpful for the replica exchange part of the 
 * simulation. It can also be used to try different generators
 * for performance. */
#include <cmath>
#include <random>

template<typename data_t>
struct random_number_generator
{
    const std::uint32_t my_seed;
    std::mt19937 generator;
    std::uniform_real_distribution<data_t> u_real;
    // TODO: Should I add a uniform int distribution too? 
    
    random_number_generator(const std::uint32_t _seed) : my_seed(_seed), generator(my_seed), u_real(0., 1.) {}
    ~random_number_generator() {}
   
    data_t operator () () { return u_real(generator); }
    data_t get_real () { return u_real(generator); }
    data_t get_real_modified(const data_t min, const data_t max){ return min + (max - min) * u_real(generator); }
    
};

#endif
