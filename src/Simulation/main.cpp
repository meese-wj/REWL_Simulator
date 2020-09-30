#include <stdio.h>
#include <string>
#include <mpi.h>

const std::string DELIMITER = "  ";

#include "rewl_simulation.hpp"

#include <file_manager.hpp>
#include <file_header.hpp>
#include <array_shift.hpp>
#include <write_microcanonical_observables.hpp>
#include <thermodynamics.hpp>
#include <write_self_averaged_observables.hpp>

using thermo_t = Thermodynamics<ENERGY_TYPE, LOGDOS_TYPE, OBS_TYPE, System_Obs_enum_t>;
constexpr ENERGY_TYPE Tmin = 0.01;
constexpr ENERGY_TYPE Tmax = 4.71;
constexpr size_t num_T = 1000;

int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);
    int world_rank;
    int world_size;
    MPI_Comm_rank( MPI_COMM_WORLD, &world_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &world_size );

    if (argc > 1)
    {
        if ( world_rank == REWL_MASTER_PROC )
        {
            printf("\nThis simulation %s takes no command line arguments.", argv[0]);
            printf("\nReturning error value.\n\n");
        }
        return 1;
    }
 
    REWL_simulation * simulation = new REWL_simulation();

    /* ****************************************************************************** */
    /* Set up the file system for this simulation.                                    */
    /* ****************************************************************************** */
    
    System_Strings sys_strings = System_Strings();
    REWL_Parameter_String rewl_strings = REWL_Parameter_String();
    std::filesystem::path data_path;
    std::string data_file_header;

    if ( world_rank == REWL_MASTER_PROC )
    {
        data_path = create_output_path( sys_strings.model_name, sys_strings.size_string ); 
        data_file_header = create_file_header( sys_strings.file_header, rewl_strings.file_header );
    }
    
    /* ****************************************************************************** */
    
    if ( world_rank == REWL_MASTER_PROC )
        printf("\nStarting simulation...");
    MPI_Barrier(MPI_COMM_WORLD);
#if PRINT_HISTOGRAM
    simulation -> simulate( data_path / "Histograms" / sys_strings.size_string );
#else
    simulation -> simulate();
#endif
    if ( world_rank == REWL_MASTER_PROC )
        printf("\nEnd of simulation. Exiting.\n\n");
    MPI_Barrier(MPI_COMM_WORLD);

    if ( world_rank == REWL_MASTER_PROC )
        printf("\nExporting energy array, logDoS array, and observables array.\n");
   
    ENERGY_TYPE * final_energy_array = nullptr;
    LOGDOS_TYPE * final_logdos_array = nullptr;
    OBS_TYPE * final_observable_array = nullptr;

    const size_t final_num_bins = simulation -> my_walker -> wl_walker.wl_histograms.num_bins;

    // Copy out the final logdos and the observables
    simulation -> my_walker -> export_energy_bins( final_energy_array );
    simulation -> my_walker -> wl_walker.wl_histograms.export_logdos( final_logdos_array );
    simulation -> my_walker -> system_obs.export_observables( final_observable_array );

    // Adjust the logdos by the ground state degeneracy
    array_shift_by_value( System_Parameters::ground_state_degeneracy - final_logdos_array[0], final_num_bins, final_logdos_array );
    MPI_Barrier(MPI_COMM_WORLD);

    // TODO: Send arrays to master for concatenation
    // TODO: Concatenate microcanonical observables

    if ( world_rank == REWL_MASTER_PROC )
    {
        // Print out the microcanonical observables before thermally averaging
        write_microcanonical_observables<ENERGY_TYPE, LOGDOS_TYPE, OBS_TYPE>( System_Parameters::N, final_num_bins, convert<System_Obs_enum_t>(System_Obs_enum_t::NUM_OBS),
                                                                              convert<System_Obs_enum_t>(System_Obs_enum_t::counts_per_bin),
                                                                              sys_strings.file_name_base, data_file_header, System_Obs::string_names,
                                                                              data_path, final_energy_array, final_logdos_array, final_observable_array ); 

        thermo_t * thermo = new thermo_t ( final_num_bins, Tmin, Tmax, num_T );
        
        printf("\nNow calculating canonical thermodynamics.\n");
#if COLLECT_TIMINGS
        printf("Starting new timer.");
        fflush(stdout);
        auto timer_start = std::chrono::high_resolution_clock::now();
        auto timer_end = timer_start;
#endif
        thermo -> calculate_thermodynamics( System_Parameters::N, final_energy_array, final_logdos_array, final_observable_array ); 

        printf("\nWriting thermodynamics to file.");

        constexpr size_t total_observables =   convert<Energy_Obs::enum_names>(Energy_Obs::enum_names::NUM_OBS) 
                                             + convert<System_Obs_enum_t>(System_Obs_enum_t::NUM_OBS);
        
        const std::vector<std::string> thermal_obs_names = concatenate_vector_string( Energy_Obs::string_names, System_Obs::string_names );

        write_observables_to_file<ENERGY_TYPE, OBS_TYPE>( num_T, total_observables, sys_strings.file_name_base, data_file_header,
                                                          thermal_obs_names, data_path, thermo -> temperatures, thermo -> canonical_observables ); 
        
#if COLLECT_TIMINGS
        timer_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<float> time_elapsed = timer_end - timer_start;
        printf("\nTime elapsed: %e seconds.\n", time_elapsed.count());
#endif
        
        delete thermo;
    }

    /* ************************************************************************************* */
    /* Now it is time to clean up the heap.                                                  */
    /* ************************************************************************************* */

    delete simulation;
   
    delete [] final_energy_array;
    delete [] final_logdos_array;
    delete [] final_observable_array;

    MPI_Barrier(MPI_COMM_WORLD);    
    MPI_Finalize();
    return 0;
}
