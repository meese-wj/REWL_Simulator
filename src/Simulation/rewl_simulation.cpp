#include "rewl_simulation.hpp"

void REWL_simulation::simulate()
{
    bool simulation_incomplete = true;
    
    while (simulation_incomplete)
    {
        // First update the walker up until the 
        // sweeps_per_check
        wang_landau_walk(REWL_Parameters::sweeps_per_check); 

        // Now check to see if the histogram is flat
        if ( my_walker -> wl_walker.is_flat( REWL_Parameters::flatness_criterion ) )
        {
            // Reset only the energy histogram and leave
            // the logdos alone.
            my_walker -> wl_walker.reset_histogram();

            // TODO: Generalize to 1/t algorithm.
            my_walker -> incrementer *= 0.5;

            simulation_incomplete = ( my_walker -> incrementer < REWL_Parameters::final_increment );
        }
    }
}
