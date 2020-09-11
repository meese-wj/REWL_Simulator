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

    thermo_t * thermo = new thermo_t ( System_Parameters::energy_min, System_Parameters::energy_max, System_Parameters::energy_bin_size, 0.1, 4.7, 1000 );
    
    LOGDOS_TYPE * final_logdos_array = nullptr;
    thermo -> calculate_thermodynamics( System_Parameters::N, final_logdos_array, 
                                        final_observable_array ); 

    delete simulation;
    delete [] final_logdos_array;
    delete [] final_observable_array;
    delete thermo;

    return 0;
}
