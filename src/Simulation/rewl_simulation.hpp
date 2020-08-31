#ifndef REWL_SIMULATION
/* This file contains the templated simulation
 * struct that will initialize the simulation
 * and contain all the relevant types. */

#include <rewl_walker.hpp>

// Create an alias for the observable
// data type. 
using LOGDOS_TYPE = double;
using OBS_TYPE = float;

struct REWL_simulation
{
    int my_world_rank;

    // TODO: Make a glazier and start this stuff up.
    
    REWL_Walker<LOGDOS_TYPE, OBS_TYPE> * my_walker = nullptr;

    
    ~REWL_simulation() { if (my_walker != nullptr) delete my_walker; }
};

#endif
