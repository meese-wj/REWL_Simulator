#include <stdio.h>

#include "rewl_simulation.hpp"
#include <thermodynamics.hpp>

using thermo_t = Thermodynamics<ENERGY_TYPE, LOGDOS_TYPE, OBS_TYPE, System_Obs_enum_t, Observables_t>;

int main(const int argc, const char * argv[])
{
    if (argc > 1)
    {
        printf("\nThis simulation %s takes no command line arguments.", argv[0]);
        printf("\nReturning error value.\n\n");
        return 1;
    }


    REWL_simulation * simulation = new REWL_simulation();
   
    printf("\nStarting simulation...");
    simulation -> simulate();
    printf("\nEnd of simulation. Exiting.\n\n");

    printf("\nNow calculating canonical thermodynamics.");
    
    LOGDOS_TYPE * final_logdos_array = nullptr;
    OBS_TYPE * final_observable_array = nullptr;

    // Copy out the final logdos and the observables
    simulation -> my_walker -> wl_walker.wl_histograms.export_logdos( final_logdos_array );
    simulation -> my_walker -> system_obs.export_observables( final_observable_array );

    thermo_t * thermo = new thermo_t ( System_Parameters::energy_min, System_Parameters::energy_max, System_Parameters::energy_bin_size, 0.1, 4.7, 1000 );
    
    thermo -> calculate_thermodynamics( System_Parameters::N, final_logdos_array, 
                                        final_observable_array ); 


    /* ************************************************************************************* */
    /* Now it is time to clean up the heap.                                                  */
    /* ************************************************************************************* */

    delete simulation;
    delete thermo;
    
    delete [] final_logdos_array;
    delete [] final_observable_array;

    return 0;
}
