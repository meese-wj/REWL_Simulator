#include "main_headers.hpp"

using thermo_t = Thermodynamics<ENERGY_TYPE, LOGDOS_TYPE, OBS_TYPE, System_Obs_enum_t>;
constexpr ENERGY_TYPE Tmin = 0.01;
constexpr ENERGY_TYPE Tmax = 6.71;
constexpr size_t num_T = 3000;

#if MPI_ON
int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);
    int world_rank = 0;
    int world_size = 0;
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
    MPI_Barrier(MPI_COMM_WORLD);
 
    REWL_simulation * simulation = new REWL_simulation();

    if ( world_rank == REWL_MASTER_PROC )
        printf("\nStarting simulation...\n");
    MPI_Barrier(MPI_COMM_WORLD);
#if PRINT_HISTOGRAM
    simulation -> simulate( data_path / "Histograms" / sys_strings.size_string );
#else
    simulation -> simulate();
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    if ( world_rank == REWL_MASTER_PROC )
        printf("\nEnd of simulation. Exiting.\n\n");

    if ( world_rank == REWL_MASTER_PROC )
        printf("\nExporting energy array, logDoS array, and observables array.\n");
   
    ENERGY_TYPE * final_energy_array = nullptr;
    LOGDOS_TYPE * final_logdos_array = nullptr;
    OBS_TYPE * final_observable_array = nullptr;

    size_t final_num_bins = simulation -> my_walker -> wl_walker.wl_histograms.num_bins;
    size_t final_num_obs_values = convert<System_Obs_enum_t>(System_Obs_enum_t::NUM_OBS) * final_num_bins;

    // Copy out the final logdos and the observables
    simulation -> my_walker -> export_energy_bins( final_energy_array );
    simulation -> my_walker -> wl_walker.wl_histograms.export_logdos( final_logdos_array );
    simulation -> my_walker -> system_obs.export_observables( final_observable_array );

    // Adjust the logdos by the ground state degeneracy
    array_shift_by_value( System_Parameters::ground_state_degeneracy - final_logdos_array[0], final_num_bins, final_logdos_array );
    MPI_Barrier(MPI_COMM_WORLD);

    printf("\nBefore thermodynamics with process %d\n", world_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if ( world_rank == REWL_MASTER_PROC )
    {

        /* ****************************************************************************** */
        /* Set up the file system for this simulation.                                    */
        /* ****************************************************************************** */
        
        System_Strings sys_strings = System_Strings();
        REWL_Parameter_String rewl_strings = REWL_Parameter_String();
        std::filesystem::path data_path = create_output_path( sys_strings.model_name, sys_strings.size_string ); 
        std::string data_file_header = create_file_header( sys_strings.file_header, rewl_strings.file_header );
 
        /* ****************************************************************************** */

        // Send arrays to master for concatenation
        table<ENERGY_TYPE> energy_table  ( world_size );
        table<LOGDOS_TYPE> logdos_table  ( world_size );
        table<OBS_TYPE> observable_table ( world_size );

        energy_table[ world_rank ] = std::vector<ENERGY_TYPE> ( final_energy_array, final_energy_array + final_num_bins ); 
        logdos_table[ world_rank ] = std::vector<LOGDOS_TYPE> ( final_logdos_array, final_logdos_array + final_num_bins );
        observable_table[ world_rank ] = std::vector<OBS_TYPE> ( final_observable_array, final_observable_array + final_num_obs_values );

        for ( int proc = 0; proc != world_size; ++proc )
        {
            if ( proc != world_rank )
            {
                MPI_Status status; 
                mpi_recv_array_to_vector<ENERGY_TYPE>( proc, energy_table[proc], MPI_FLOAT, final_energy_tag, MPI_COMM_WORLD, &status );
                mpi_recv_array_to_vector<LOGDOS_TYPE>( proc, logdos_table[proc], MPI_LOGDOS_TYPE, final_logdos_tag, MPI_COMM_WORLD, &status );
                mpi_recv_array_to_vector<OBS_TYPE>( proc, observable_table[proc], MPI_OBS_TYPE, final_obs_tag, MPI_COMM_WORLD, &status );
            }
        }
        
        std::vector<ENERGY_TYPE> final_energy_vector;
        std::vector<LOGDOS_TYPE> final_logdos_vector;
        std::vector<OBS_TYPE>    final_observable_vector;
        // TODO: Concatenate microcanonical observables
        concatenate_tables<ENERGY_TYPE, LOGDOS_TYPE, OBS_TYPE>
            ( REWL_Parameters::window_overlap == static_cast<float>( single_bin_overlap ),
              convert<System_Obs_enum_t>(System_Obs_enum_t::NUM_OBS),
              convert<System_Obs_enum_t>(System_Obs_enum_t::counts_per_bin),
              energy_table, logdos_table, observable_table,
              final_energy_vector, final_logdos_vector, final_observable_vector );
 
        // Finally, after concatenation, reset the final number of bins
        final_num_bins = final_energy_vector.size();
        final_num_obs_values = final_observable_vector.size();

        // Print out the microcanonical observables before thermally averaging
        write_microcanonical_observables<ENERGY_TYPE, LOGDOS_TYPE, OBS_TYPE>( System_Parameters::N, final_num_bins, convert<System_Obs_enum_t>(System_Obs_enum_t::NUM_OBS),
                                                                              convert<System_Obs_enum_t>(System_Obs_enum_t::counts_per_bin),
                                                                              sys_strings.file_name_base, data_file_header, System_Obs::string_names,
                                                                              data_path, final_energy_vector.data(), final_logdos_vector.data(), final_observable_vector.data() ); 

        thermo_t * thermo = new thermo_t ( final_num_bins, Tmin, Tmax, num_T );
        
        printf("\nNow calculating canonical thermodynamics.\n");
#if COLLECT_TIMINGS
        printf("Starting new timer.");
        fflush(stdout);
        auto timer_start = std::chrono::high_resolution_clock::now();
        auto timer_end = timer_start;
#endif
        thermo -> calculate_thermodynamics( System_Parameters::N, final_energy_vector.data(), final_logdos_vector.data(), final_observable_vector.data() ); 

        OBS_TYPE * nonlinear_obs_array = nullptr;
        calculate_nonlinear_observables<OBS_TYPE, thermo_t>( num_T, System_Parameters::N, thermo, nonlinear_obs_array ); 


        printf("\nWriting thermodynamics to file.");

        constexpr size_t total_observables =   convert<Energy_Obs::enum_names>(Energy_Obs::enum_names::NUM_OBS) 
                                             + convert<System_Obs_enum_t>(System_Obs_enum_t::NUM_OBS);
        
        const std::vector<std::string> thermal_obs_names = concatenate_vector_string( Energy_Obs::string_names, System_Obs::string_names );

        write_observables_to_file<ENERGY_TYPE, OBS_TYPE>( num_T, total_observables, sys_strings.file_name_base, data_file_header,
                                                          thermal_obs_names, data_path, thermo -> temperatures, thermo -> canonical_observables ); 
    
        write_nonlinear_obs_to_file<ENERGY_TYPE, OBS_TYPE>( num_T, convert(System_Obs::nonlinear_obs_enum::NUM_OBS), sys_strings.file_name_base, data_file_header,
                                                            System_Obs::nonlinear_obs_strings, data_path, thermo -> temperatures, nonlinear_obs_array ); 
        
#if COLLECT_TIMINGS
        timer_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<float> time_elapsed = timer_end - timer_start;
        printf("\nTime elapsed: %e seconds.\n", time_elapsed.count());
#endif
        
        delete thermo;
        delete [] nonlinear_obs_array;
    }
    else
    {
        // Send arrays to the master process for analysis.
        mpi_send_array<ENERGY_TYPE>( REWL_MASTER_PROC, final_num_bins, final_energy_array, MPI_FLOAT, final_energy_tag, MPI_COMM_WORLD );  
        mpi_send_array<LOGDOS_TYPE>( REWL_MASTER_PROC, final_num_bins, final_logdos_array, MPI_LOGDOS_TYPE, final_logdos_tag, MPI_COMM_WORLD );  
        mpi_send_array<OBS_TYPE>( REWL_MASTER_PROC, final_num_obs_values, final_observable_array, MPI_OBS_TYPE, final_obs_tag, MPI_COMM_WORLD );  
    }
    MPI_Barrier(MPI_COMM_WORLD);

    printf("\nCleaning up the heap with process %d\n", world_rank);
    MPI_Barrier(MPI_COMM_WORLD);

    /* ************************************************************************************* */
    /* Now it is time to clean up the heap.                                                  */
    /* ************************************************************************************* */

    delete simulation;
   
    delete [] final_energy_array;
    delete [] final_logdos_array;
    delete [] final_observable_array;

    MPI_Barrier(MPI_COMM_WORLD);    
    printf("\nCompletely done with the simulation on process %d\n", world_rank);
    fflush(stdout);
    int ret_val = MPI_Finalize();
    printf("\nID %d: Finalize return = %d", world_rank, ret_val);
    fflush(stdout);
    return 0;
}
#else
int main(int argc, char * argv[])
{
    int world_rank = 0;

    if (argc > 1)
    {
        printf("\nThis simulation %s takes no command line arguments.", argv[0]);
        printf("\nReturning error value.\n\n");
        return 1;
    }

    /* ****************************************************************************** */
    /* Set up the file system for this simulation.                                    */
    /* ****************************************************************************** */
    
    System_Strings sys_strings = System_Strings();
    REWL_Parameter_String rewl_strings = REWL_Parameter_String();
    std::filesystem::path data_path = create_output_path( sys_strings.model_name, sys_strings.size_string ); 
    std::string data_file_header = create_file_header( sys_strings.file_header, rewl_strings.file_header );

    /* ****************************************************************************** */
 
    REWL_simulation * simulation = new REWL_simulation();

    printf("\nStarting simulation...");
#if PRINT_HISTOGRAM
    simulation -> simulate( data_path / "Histograms" / sys_strings.size_string );
#else
    simulation -> simulate();
#endif
    printf("\nEnd of simulation. Exiting.\n\n");
    
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

    printf("\nBefore thermodynamics with process %d\n", world_rank);
       
    // Print out the microcanonical observables before thermally averaging
    write_microcanonical_observables<ENERGY_TYPE, LOGDOS_TYPE, OBS_TYPE>( System_Parameters::N, final_num_bins,
                                                                          convert<System_Obs_enum_t>(System_Obs_enum_t::NUM_OBS),
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

    OBS_TYPE * nonlinear_obs_array = nullptr;
    calculate_nonlinear_observables<OBS_TYPE, thermo_t>( num_T, System_Parameters::N, thermo, nonlinear_obs_array ); 

    printf("\nWriting thermodynamics to file.");

    constexpr size_t total_observables =   convert<Energy_Obs::enum_names>(Energy_Obs::enum_names::NUM_OBS) 
                                         + convert<System_Obs_enum_t>(System_Obs_enum_t::NUM_OBS);
    
    const std::vector<std::string> thermal_obs_names = concatenate_vector_string( Energy_Obs::string_names, System_Obs::string_names );

    write_observables_to_file<ENERGY_TYPE, OBS_TYPE>( num_T, total_observables, sys_strings.file_name_base, data_file_header,
                                                      thermal_obs_names, data_path, thermo -> temperatures, thermo -> canonical_observables ); 
   
    write_nonlinear_obs_to_file<ENERGY_TYPE, OBS_TYPE>( num_T, convert(System_Obs::nonlinear_obs_enum::NUM_OBS), sys_strings.file_name_base, data_file_header,
                                                        System_Obs::nonlinear_obs_strings, data_path, thermo -> temperatures, nonlinear_obs_array ); 
    
#if COLLECT_TIMINGS
    timer_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> time_elapsed = timer_end - timer_start;
    printf("\nTime elapsed: %e seconds.\n", time_elapsed.count());
#endif
    
    printf("\nCleaning up the heap with process %d\n", world_rank);

    /* ************************************************************************************* */
    /* Now it is time to clean up the heap.                                                  */
    /* ************************************************************************************* */

    delete thermo;
    delete simulation;
  
    delete [] nonlinear_obs_array;
    delete [] final_energy_array;
    delete [] final_logdos_array;
    delete [] final_observable_array;

    return 0;
}
#endif
