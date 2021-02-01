#include "main_headers.hpp"

// TODO: This probably will break
// the histograms.
#if JOB_ARRAYS
constexpr int num_args = 2;
constexpr int job_id_index = 1;
#else
constexpr int num_args = 1;
#endif

using thermo_t = Thermodynamics<ENERGY_TYPE, LOGDOS_TYPE, OBS_TYPE, System_Obs_enum_t>;
constexpr ENERGY_TYPE Tmin = 0.01;
constexpr ENERGY_TYPE Tmax = 20.01;
constexpr size_t num_T = 20000;

#if MPI_ON
int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);
    int world_rank = 0;
    int world_size = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, &world_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &world_size );

    if ( argc != num_args )
    {
        if ( world_rank == REWL_MASTER_PROC )
        {
#if JOB_ARRAYS
            printf("\nThis simulation %s takes 1 command line argument as the Job ID.", argv[0]);
            printf("\nReturning error value.\n\n");
#else
            printf("\nThis simulation %s takes no command line arguments.", argv[0]);
            printf("\nReturning error value.\n\n");
#endif
        }
        MPI_Finalize();
        return 1;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    ENERGY_TYPE ground_state_energy = 0.;
    ENERGY_TYPE highest_energy = 0.;
    Hamiltonian_t<OBS_TYPE> * toy_model = nullptr;

    if ( world_rank == REWL_MASTER_PROC )
    {
        toy_model = new Hamiltonian_t<OBS_TYPE> ();

#if SIMULATED_ANNEALING
        Simulated_Annealer<ENERGY_TYPE, Hamiltonian_t<OBS_TYPE>, State_t<OBS_TYPE> > * annealer = new Simulated_Annealer<ENERGY_TYPE, Hamiltonian_t<OBS_TYPE>, State_t<OBS_TYPE> >( (size_t)1e3, 50., 0.05, 250 );
        printf("\nSimulated annealer started. Initial energy = %e", toy_model -> current_state.energy);        
        annealer -> simulate_annealing( System_Parameters::N, System_Parameters::num_DoF / System_Parameters::N, toy_model );

        delete annealer;
#endif

        toy_model -> print_lattice();

        ground_state_energy = toy_model -> current_state.energy;
        highest_energy = System_Parameters::energy_max;
    }

    MPI_Bcast(&ground_state_energy, 1, MPI_ENERGY_TYPE, REWL_MASTER_PROC, MPI_COMM_WORLD);
    MPI_Bcast(&highest_energy, 1, MPI_ENERGY_TYPE, REWL_MASTER_PROC, MPI_COMM_WORLD);
    OBS_TYPE * dof_field_array = new OBS_TYPE [System_Parameters::num_DoF]();
    if ( world_rank == REWL_MASTER_PROC )
    {
        for ( size_t idx = 0; idx != System_Parameters::num_DoF; ++idx )
            dof_field_array[idx] = toy_model -> spin_array[idx];
    }
    MPI_Bcast( dof_field_array, System_Parameters::num_DoF, MPI_OBS_TYPE, REWL_MASTER_PROC, MPI_COMM_WORLD );

#if RANDOM_DISORDER
    ENERGY_TYPE * disorder_array = new ENERGY_TYPE [System_Parameters::N]();
    if ( world_rank == REWL_MASTER_PROC )
    {
        for ( size_t idx = 0; idx != System_Parameters::N; ++idx )
            disorder_array[idx] = toy_model -> field_array[idx];
    }
    MPI_Bcast( disorder_array, System_Parameters::N, MPI_ENERGY_TYPE, REWL_MASTER_PROC, MPI_COMM_WORLD );
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    delete toy_model;

    /* ****************************************************************************** */
    /* Set up the file system for this simulation.                                    */
    /* ****************************************************************************** */
    
    // Grab today's date at the start of the simulation
    std::string todays_date = get_todays_date();

#if JOB_ARRAYS
    std::string job_id_string ( argv[job_id_index] );
    System_Strings sys_strings = System_Strings( job_id_string );
#else
    System_Strings sys_strings = System_Strings();
#endif
    sys_strings.energy_min = std::to_string(ground_state_energy);
    sys_strings.energy_max = std::to_string(highest_energy);
    sys_strings.num_bins = std::to_string(static_cast<size_t>((highest_energy - ground_state_energy) / System_Parameters::energy_bin_size));
    sys_strings.update_file_header();
    
    REWL_Parameter_String rewl_strings = REWL_Parameter_String();
    // Update the parameter string with the world size
    rewl_strings.num_walkers = std::to_string(world_size);
    rewl_strings.recreate_parameter_string();
    std::filesystem::path data_path;
    std::string data_file_header = create_file_header( sys_strings.file_header, rewl_strings.file_header );
    if ( world_rank == REWL_MASTER_PROC )
    {
        std::cout << "\n**************************************************************************************\n";
        std::cout << "\n" << todays_date << "\n\n" << data_file_header; 
        std::cout << "\n**************************************************************************************\n";
#if JOB_ARRAYS
        data_path = create_output_path( sys_strings.model_name, todays_date, sys_strings.size_string, job_id_string );
#else
        data_path = create_output_path( sys_strings.model_name, todays_date, sys_strings.size_string ); 
#endif
    }

    /* ****************************************************************************** */
    MPI_Barrier(MPI_COMM_WORLD);

#ifndef INDEPENDENT_WALKERS
    /* **************************************************************************** */
    /* Set up the communicators                                                     */
    /*     Each processor belongs to four communicators:                            */
    /*         * MPI_COMM_WORLD                                                     */
    /*         * Their own energy window                                            */
    /*         * Their even communicator                                            */
    /*         * Their odd communicator                                             */
    /* **************************************************************************** */

    int my_ids_per_comm [ Communicators::NUM_COMMS ];     // IDs of the walker in 
                                                          // each communicator
                                                          // except for WORLD.
    int my_comm_ids [ Communicators::NUM_COMMS ];         // IDs of each of the 
                                                          // communicator the walker
                                                          // belongs to.
    // Set up the global group
    MPI_Group world_group;
    MPI_Comm_group( MPI_COMM_WORLD, &world_group );         

    // Define the local communicators
    int num_local_comms = ( world_size / REWL_Parameters::replicas_per_window ) - 1;
    MPI_Group * local_groups  = new MPI_Group [ num_local_comms ];
    MPI_Comm  * local_comms   = new MPI_Comm  [ num_local_comms ];
    MPI_Group * window_groups = new MPI_Group [ world_size / REWL_Parameters::replicas_per_window ];
    MPI_Comm  * window_comms  = new MPI_Comm  [ world_size / REWL_Parameters::replicas_per_window ];

    // Create the local groups and communicators for a single window
    define_window_communicators( world_size, REWL_Parameters::replicas_per_window, world_group,
                                 window_groups, window_comms ); 

    // Create the local groups and communicators between windows
    create_local_groups_and_communicators( world_size, REWL_Parameters::replicas_per_window, world_group,
                                           local_groups, local_comms );

    // Populate the IDs for the walkers and communicators
    determine_my_local_IDs( world_rank, world_size, REWL_Parameters::replicas_per_window, 
                            my_ids_per_comm, local_comms, window_comms ); 

    determine_my_communicators( world_rank, world_size, REWL_Parameters::replicas_per_window, my_comm_ids );

    MPI_Barrier(MPI_COMM_WORLD);
    /* **************************************************************************** */
#endif

    REWL_simulation * simulation = new REWL_simulation(ground_state_energy, highest_energy);
    simulation -> my_walker -> system.import_DoFs( dof_field_array );
    delete dof_field_array;

#if RANDOM_DISORDER
    simulation -> my_walker -> system.import_disorder( disorder_array );
    delete [] disorder_array;
#endif

    simulation -> my_walker -> adjust_state_to_range();

    if ( world_rank == REWL_MASTER_PROC )
        printf("\nStarting simulation...\n");
    MPI_Barrier(MPI_COMM_WORLD);
#if PRINT_HISTOGRAM
    // TODO: This won't work for job arrays
    simulation -> simulate( data_path / "Histograms" / sys_strings.size_string );
#else
#ifndef INDEPENDENT_WALKERS
    simulation -> simulate( my_ids_per_comm, my_comm_ids, local_comms, window_comms );
#else
    simulation -> simulate();
#endif
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

    // Find out which ranks are 0-processors
    int * window_ids      = new int [ world_size ] ();
    int * window_comm_ids = new int [ world_size ] ();
    MPI_Gather( &( my_ids_per_comm[ Communicators::window_comm ] ), 1, MPI_INT, window_ids, 1, MPI_INT, REWL_MASTER_PROC, MPI_COMM_WORLD );
    MPI_Gather( &( my_comm_ids[ Communicators::window_comm ] ), 1, MPI_INT, window_comm_ids, 1, MPI_INT, REWL_MASTER_PROC, MPI_COMM_WORLD );
 
    if ( world_rank == REWL_MASTER_PROC )
    {

                // Send arrays to master for concatenation
        table<ENERGY_TYPE> energy_table  ( world_size / REWL_Parameters::replicas_per_window );
        table<LOGDOS_TYPE> logdos_table  ( world_size / REWL_Parameters::replicas_per_window );
        table<OBS_TYPE> observable_table ( world_size / REWL_Parameters::replicas_per_window );

        energy_table[ world_rank ] = std::vector<ENERGY_TYPE> ( final_energy_array, final_energy_array + final_num_bins ); 
        logdos_table[ world_rank ] = std::vector<LOGDOS_TYPE> ( final_logdos_array, final_logdos_array + final_num_bins );
        observable_table[ world_rank ] = std::vector<OBS_TYPE> ( final_observable_array, final_observable_array + final_num_obs_values );

        for ( int proc = 0; proc != world_size; ++proc )
        {
            // Only communicate with 0-processors per window
            if ( proc != world_rank && window_ids[proc] == 0 )
            {
                int window_idx = window_comm_ids[ proc ];
                MPI_Status status; 
                mpi_recv_array_to_vector<ENERGY_TYPE>( proc,  energy_table[window_idx], MPI_FLOAT, final_energy_tag, MPI_COMM_WORLD, &status );
                mpi_recv_array_to_vector<LOGDOS_TYPE>( proc,  logdos_table[window_idx], MPI_LOGDOS_TYPE, final_logdos_tag, MPI_COMM_WORLD, &status );
                mpi_recv_array_to_vector<OBS_TYPE>( proc, observable_table[window_idx], MPI_OBS_TYPE, final_obs_tag, MPI_COMM_WORLD, &status );
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
    else if ( my_ids_per_comm[ Communicators::window_comm ] == 0 ) 
    {
        // Send arrays to the master process for analysis.
        mpi_send_array<ENERGY_TYPE>( REWL_MASTER_PROC, final_num_bins, final_energy_array, MPI_FLOAT, final_energy_tag, MPI_COMM_WORLD );  
        mpi_send_array<LOGDOS_TYPE>( REWL_MASTER_PROC, final_num_bins, final_logdos_array, MPI_LOGDOS_TYPE, final_logdos_tag, MPI_COMM_WORLD );  
        mpi_send_array<OBS_TYPE>( REWL_MASTER_PROC, final_num_obs_values, final_observable_array, MPI_OBS_TYPE, final_obs_tag, MPI_COMM_WORLD );  
    }
    // else { I'm a processor that's now bored... }
    MPI_Barrier(MPI_COMM_WORLD);

    printf("\nCleaning up the heap with process %d\n", world_rank);
    MPI_Barrier(MPI_COMM_WORLD);

    /* ************************************************************************************* */
    /* Now it is time to clean up the heap.                                                  */
    /* ************************************************************************************* */

    delete simulation;
   
    delete [] window_ids;
    delete [] final_energy_array;
    delete [] final_logdos_array;
    delete [] final_observable_array;

#ifndef INDEPENDENT_WALKERS
    delete [] local_groups;
    delete [] local_comms;
    delete [] window_groups;
    delete [] window_comms;
#endif

    MPI_Barrier(MPI_COMM_WORLD);    
    //printf("\nCompletely done with the simulation on process %d\n", world_rank);
    //fflush(stdout);
    MPI_Finalize();
    //int ret_val = MPI_Finalize();
    //printf("\nID %d: Finalize return = %d", world_rank, ret_val);
    //fflush(stdout);
    return 0;
}
#else
int main(int argc, char * argv[])
{
    int world_rank = 0;

    if ( argc != num_args )
    {
#if JOB_ARRAYS
            printf("\nThis simulation %s takes 1 command line argument as the Job ID.", argv[0]);
            printf("\nReturning error value.\n\n");
#else
            printf("\nThis simulation %s takes no command line arguments.", argv[0]);
            printf("\nReturning error value.\n\n");
#endif
        return 1;
    }
    
    ENERGY_TYPE ground_state_energy = 0.;
    ENERGY_TYPE highest_energy = 0.;
    Hamiltonian_t<OBS_TYPE> * toy_model = new Hamiltonian_t<OBS_TYPE> ();

#if SIMULATED_ANNEALING
        Simulated_Annealer<ENERGY_TYPE, Hamiltonian_t<OBS_TYPE>, State_t<OBS_TYPE> > * annealer = new Simulated_Annealer<ENERGY_TYPE, Hamiltonian_t<OBS_TYPE>, State_t<OBS_TYPE> >( (size_t)1e5, 50., 0.05, 20 );
        
        annealer -> simulate_annealing( System_Parameters::N, System_Parameters::num_DoF / System_Parameters::N, toy_model );

        delete annealer;
#endif

    ground_state_energy = toy_model -> current_state.energy;
    highest_energy = System_Parameters::energy_max;
    OBS_TYPE * dof_field_array = new OBS_TYPE [System_Parameters::num_DoF]();
    for ( size_t idx = 0; idx != System_Parameters::num_DoF; ++idx )
        dof_field_array[idx] = toy_model -> spin_array[idx];

#if RANDOM_DISORDER
    ENERGY_TYPE * disorder_array = new ENERGY_TYPE [System_Parameters::N];
    for ( size_t idx = 0; idx != System_Parameters::N; ++idx )
        disorder_array[idx] = toy_model -> field_array[idx];
#endif
    delete toy_model;

    /* ****************************************************************************** */
    /* Set up the file system for this simulation.                                    */
    /* ****************************************************************************** */
    
    System_Strings sys_strings = System_Strings();
    sys_strings.energy_min = std::to_string(ground_state_energy);
    sys_strings.energy_max = std::to_string(highest_energy);
    sys_strings.num_bins = std::to_string(static_cast<size_t>((highest_energy - ground_state_energy) / System_Parameters::energy_bin_size));
    sys_strings.update_file_header();

    REWL_Parameter_String rewl_strings = REWL_Parameter_String();
#if JOB_ARRAYS
    System_Strings sys_strings = System_Strings( job_id_string );
    std::string job_id_string ( argv[job_id_index] );
    std::filesystem::path data_path = create_output_path( sys_strings.model_name, sys_strings.size_string, job_id_string ); 
#else
    System_Strings sys_strings = System_Strings();
    std::filesystem::path data_path = create_output_path( sys_strings.model_name, sys_strings.size_string ); 
#endif
    std::string data_file_header = create_file_header( sys_strings.file_header, rewl_strings.file_header );

    /* ****************************************************************************** */
 
    REWL_simulation * simulation = new REWL_simulation(ground_state_energy, highest_energy);
    simulation -> my_walker -> system.import_DoFs( dof_field_array );
    delete dof_field_array;

#if RANDOM_DISORDER
    simulation -> my_walker -> system.import_disorder( disorder_array );
    delete [] disorder_array;
#endif

    simulation -> my_walker -> adjust_state_to_range();

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
