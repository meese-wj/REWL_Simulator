#ifndef REWL_SIMULATION
#define REWL_SIMULATION
/* This file contains the templated simulation
 * struct that will initialize the simulation
 * and contain all the relevant types. */
#include <string>
#include <chrono>
#include <algorithm> // For std::shuffle

#include <histogram_index.hpp>
#include <glazier.hpp>
#include <rewl_walker.hpp>
#include <rewl_parameter_string.hpp>
#include <mpi_rewl_comm_setup.hpp>

struct REWL_simulation
{
#if MPI_ON
    int my_world_rank = REWL_MASTER_PROC;
#else
    int my_world_rank = 0;
#endif
    bool i_am_the_master = true;

    // TODO: Make a glazier and start this stuff up.
    glazier<ENERGY_TYPE, histogram_index<ENERGY_TYPE> > * window_maker = nullptr;
    
    REWL_Walker<ENERGY_TYPE, LOGDOS_TYPE, OBS_TYPE, histogram_index<ENERGY_TYPE> > * my_walker = nullptr;

    REWL_simulation();

    ~REWL_simulation()
    { 
        if (window_maker != nullptr) delete window_maker;
        if (my_walker != nullptr) delete my_walker;
    }

#ifndef INDEPENDENT_WALKERS
#if SAMPLE_AFTER
    void replica_exchange_update( int & exchange_direction, const size_t iteration_counter, const bool sample_observables, 
                                  const int * const my_ids_per_comm, const int * const my_comm_ids, MPI_Comm * const local_communicators ) const;
#else
    void replica_exchange_update( int & exchange_direction, const size_t iteration_counter, 
                                  const int * const my_ids_per_comm, const int * const my_comm_ids, MPI_Comm * const local_communicators ) const;
#endif
    void average_and_redistribute_window( const bool simulation_incomplete, const int * my_comm_ids, MPI_Comm * const window_communicators ) const;
    bool check_if_window_is_flat( const int * my_comm_ids, MPI_Comm * const window_communicators ) const;
#endif

    void simulate(
#if PRINT_HISTOGRAM
                  const std::filesystem::path & histogram_path
#endif
#ifndef INDEPENDENT_WALKERS
                  const int * const my_ids_per_comm, const int * const my_comm_ids, MPI_Comm * const local_communicators, MPI_Comm * const window_communicators
#endif
            ) const;
};

// Simulation constructor
REWL_simulation::REWL_simulation()
{
#if MPI_ON
    MPI_Comm_rank( MPI_COMM_WORLD, &my_world_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &REWL_Parameters::num_walkers );
    i_am_the_master = ( my_world_rank == REWL_MASTER_PROC );
#endif
    
    // Construct the glazier
    window_maker = new glazier<ENERGY_TYPE, histogram_index<ENERGY_TYPE> >
                        (System_Parameters::energy_min, System_Parameters::energy_max,
                         System_Parameters::energy_bin_size, 
                         static_cast<size_t>(REWL_Parameters::num_walkers) / REWL_Parameters::replicas_per_window,
                         REWL_Parameters::replicas_per_window, 
                         static_cast<ENERGY_TYPE>(REWL_Parameters::window_overlap));

    window_maker -> construct_windows();

    // Construct the walker
    ENERGY_TYPE walker_min = window_maker -> all_windows[my_world_rank].minimum;
    ENERGY_TYPE walker_max = window_maker -> all_windows[my_world_rank].maximum;
    ENERGY_TYPE walker_bin_size = window_maker -> all_windows[my_world_rank].bin_size;
    size_t walker_num_bins = window_maker -> all_windows[my_world_rank].num_bins;

    // TODO: Set up the timer as the seed.
    std::uint32_t walker_seed = static_cast<std::uint32_t> (1);

    my_walker = new REWL_Walker<ENERGY_TYPE, LOGDOS_TYPE, OBS_TYPE, histogram_index<ENERGY_TYPE> >
                (walker_min, walker_max, walker_bin_size, walker_num_bins, walker_seed);
}

#if INDEPENDENT_WALKERS /* This will not turn on replica exchange */

// Main function for the simulation.
void REWL_simulation::simulate(
#if PRINT_HISTOGRAM
                                const std::filesystem::path & histogram_path
#endif
                               ) const
{
#if COLLECT_TIMINGS
    auto start = std::chrono::high_resolution_clock::now();
    auto iteration_start = start;
    auto timer = start;
#endif

    size_t iteration_counter = 1;
    size_t sweep_counter = 0;
    bool simulation_incomplete = true;
#if SAMPLE_AFTER
    bool sample_observables = false;
#endif
    
    while (simulation_incomplete)
    {
        // First update the walker up until the 
        // sweeps_per_check
#if SAMPLE_AFTER
        my_walker -> wang_landau_walk(REWL_Parameters::sweeps_per_check, sample_observables); 
#else
        my_walker -> wang_landau_walk(REWL_Parameters::sweeps_per_check); 
#endif
        sweep_counter += REWL_Parameters::sweeps_per_check;

#if PRINT_HISTOGRAM
        if ( sweep_counter % (REWL_Parameters::sweeps_per_check) == 0 )
            my_walker -> wl_walker.wl_histograms.print_histogram_counts(iteration_counter, histogram_path);
#endif

        // Now check to see if the histogram is flat
        if ( my_walker -> wl_walker.is_flat( REWL_Parameters::flatness_criterion ) )
        {
            // Reset only the energy histogram and leave
            // the logdos alone.
#if PRINT_HISTOGRAM
            my_walker -> wl_walker.wl_histograms.print_histogram_counts(iteration_counter, histogram_path);
#endif

            my_walker -> wl_walker.reset_histogram();

            printf("\nID %d: incrementer before = %e", my_world_rank, my_walker -> incrementer);
            // TODO: Generalize to 1/t algorithm.
            my_walker -> incrementer *= 0.5;

            printf("\nID %d: incrementer after = %e", my_world_rank, my_walker -> incrementer);
            
#if SAMPLE_AFTER
            if ( my_walker -> incrementer < REWL_Parameters::final_increment )
            {
                if (sample_observables == false) 
                {
                    sample_observables = true;
                    my_walker -> incrementer = 1.;   // Reset the incrementer and walk again
                    simulation_incomplete = true;
                }
                else simulation_incomplete = false;  // Kill the simulation if it is complete after sampling
            }
            else simulation_incomplete = true;
#else
            simulation_incomplete = ( my_walker -> incrementer >= REWL_Parameters::final_increment );
#endif

#if COLLECT_TIMINGS
            timer = std::chrono::high_resolution_clock::now();
            std::chrono::duration<float> time_elapsed = timer - iteration_start;
            printf("\n\n\nID %d: Iteration %ld complete after %e seconds.", my_world_rank, iteration_counter, time_elapsed.count());
            printf("\nID %d: Total Sweeps = %ld = %e Wang Landau Updates.", my_world_rank, sweep_counter, static_cast<float>(sweep_counter * System_Parameters::N) );
            printf("\nID %d: Sweep rate  = %e sweeps per second", my_world_rank, static_cast<float>(sweep_counter)/time_elapsed.count());
            printf("\nID %d: Update rate = %e updates per second", my_world_rank, static_cast<float>(sweep_counter * System_Parameters::N)/time_elapsed.count());
            fflush(stdout);
            iteration_start = std::chrono::high_resolution_clock::now();
#endif
            ++iteration_counter;
            sweep_counter = 0;
        }
    }

#if COLLECT_TIMINGS
            timer = std::chrono::high_resolution_clock::now();
            std::chrono::duration<float> simulation_time = timer - start;
            printf("\n\nID %d: Total Simulation Time: %e seconds.", my_world_rank, simulation_time.count());
#endif
}

#else /* This will turn on replica exchange. */

// Single replica exchange update
void REWL_simulation::replica_exchange_update( int & exchange_direction, const size_t iteration_counter, 
#if SAMPLE_AFTER
                                               const bool sample_observables,
#endif
                                               const int * const my_ids_per_comm, const int * const my_comm_ids, MPI_Comm * const local_communicators ) const
{
    int comm_id = my_comm_ids[ exchange_direction ];
    if ( comm_id != Communicators::NONE )
    {
        // TODO: Generalize this for multiple walkers
        int * partners = new int [ 2 * REWL_Parameters::replicas_per_window ];

        // Only have the communicator master assign partners
        if ( my_ids_per_comm[ exchange_direction ] == 0 )
        {
            // Get the partner indices
            std::vector<int> shuffled_partners ( REWL_Parameters::replicas_per_window );
            for ( int partner = 0; partner != static_cast<int>(REWL_Parameters::replicas_per_window); ++partner )
                shuffled_partners[ partner ] = static_cast<int>(REWL_Parameters::replicas_per_window) + partner;
            
            // Now shuffle the partners and pair by index
            std::shuffle( shuffled_partners.begin(), shuffled_partners.end(), my_walker -> random.generator );
            //printf("\n");
            for ( size_t replica = 0; replica != REWL_Parameters::replicas_per_window; ++replica )
            {
                // Assign the replica in the lower window the shuffled_partner at replica
                partners[ replica ] = shuffled_partners[ replica ];
                // Assign the shuffled_partner at replica the replica in the lower window
                partners[ shuffled_partners[ replica ] ] = replica;
                //printf("\nreplica = %ld, partners[%ld] = %d, partners[%d] = %d", replica, replica, partners[replica], shuffled_partners[replica], partners[shuffled_partners[replica]]);
            }
            //printf("\n");

            //partners[0] = REWL_Parameters::replicas_per_window;
            //partners[2 * REWL_Parameters::replicas_per_window - 1] = 0;
        }
         
        // Scatter the partner indices to everyone
        int partner_index = Communicators::NONE; 
        MPI_Scatter( partners, 1, MPI_INT, &partner_index, 1, MPI_INT, 0, local_communicators[ comm_id ] );

        // Delete the partners array because it is no longer needed
        delete [] partners;

        // Wait for everyone to get their partners
        MPI_Barrier( local_communicators[ comm_id ] );

        // Now proceed with the exchange
        if ( partner_index != Communicators::NONE )
        {
            MPI_Status status;

            ENERGY_TYPE current_energy = my_walker -> current_energy();
            ENERGY_TYPE new_energy = current_energy;      

            // Send the new energy to the partner index to use in the exchange
            MPI_Sendrecv_replace( &new_energy, 1, MPI_ENERGY_TYPE, partner_index, 1, partner_index, 1, local_communicators[ comm_id ], &status );

            float pexchange = 0.;
            if ( my_walker -> energy_in_range( new_energy ) )
            {
                pexchange = static_cast<float>( exp( std::min(0., my_walker -> get_logdos( current_energy ) - my_walker -> get_logdos( new_energy )) ) );
            }
            //printf("\nID %d = %d: pexchange before = %e", my_world_rank, my_ids_per_comm[exchange_direction], pexchange);

            bool we_do_exchange = false;
            if ( my_ids_per_comm[ exchange_direction ] < static_cast<int>(REWL_Parameters::replicas_per_window) )
            {
                // Have the lower ids be the calculator
                float other_pexchange = 0.;
                MPI_Recv( &other_pexchange, 1, MPI_FLOAT, partner_index, 2, local_communicators[ comm_id ], &status );

                // The exchange probability comes from the product of both
                // energy moves
                pexchange *= other_pexchange;

                we_do_exchange = ( pexchange != 0. && ( my_walker -> get_rand() < pexchange ) );

                // Send the whether the result is made to the partner
                MPI_Send( &we_do_exchange, 1, MPI_CXX_BOOL, partner_index, 3, local_communicators[ comm_id ] );
                MPI_Send( &pexchange, 1, MPI_FLOAT, partner_index, 4, local_communicators[ comm_id ] );
            }
            else
            {
                // Send the exchange probability to the calculator
                // and await a response
                MPI_Send( &pexchange, 1, MPI_FLOAT, partner_index, 2, local_communicators[ comm_id ] );
                MPI_Recv( &we_do_exchange, 1, MPI_CXX_BOOL, partner_index, 3, local_communicators[ comm_id ], &status );
                MPI_Recv( &pexchange, 1, MPI_FLOAT, partner_index, 4, local_communicators[ comm_id ], &status );
            }
            //printf("\nID %d = %d: pexchange after = %e\n", my_world_rank, my_ids_per_comm[exchange_direction], pexchange);

            if ( we_do_exchange )
            {
                //printf("\n\nID %d = %d EXCHANGED!!\n\n", my_world_rank, my_ids_per_comm[exchange_direction]);
                // Perform a MPI_Sendrecv_replace on the state and the degrees of freedom
                mpi_exchange_state<State_t<OBS_TYPE> >( my_walker -> current_state(), partner_index, comm_id, local_communicators, &status );
                mpi_exchange_DoFs<OBS_TYPE>( my_walker -> DoFs(), System_Parameters::num_DoF, partner_index, comm_id, local_communicators, &status );
            }

            // Finally, update the histograms after the exchanges
            my_walker -> update_histograms();
#if SAMPLE_AFTER
            my_walker -> update_observables( sample_observables ); 
#else
            my_walker -> update_observables();
#endif
            
        }
    }

    // Change the exchange direction
    exchange_direction = ( exchange_direction == Communicators::even_comm ? Communicators::odd_comm : Communicators::even_comm );
}

// Compute the weighted average in a window and redistribute
void REWL_simulation::average_and_redistribute_window( const bool simulation_incomplete, const int * const my_comm_ids, MPI_Comm * const window_communicators ) const
{
    if ( REWL_Parameters::replicas_per_window == 1 )
        return;

    int comm_id = my_comm_ids[ Communicators::window_comm ];

    // Store the sum of each walker's energy histogram in an array
    size_t num_bins = my_walker -> wl_walker.wl_histograms.num_bins;
    size_t * total_energy_histogram = new size_t [ num_bins ] ();
    OBS_TYPE * total_energy_histogram_obs_t = new OBS_TYPE [ num_bins ] ();

    MPI_Barrier( window_communicators[ comm_id ] );

    for ( size_t bin = 0; bin != num_bins; ++bin )
    {
        //printf("\nID %d: counts[%ld] = %ld, total_counts[%ld] = %ld", my_world_rank, bin, my_walker -> wl_walker.wl_histograms.histograms[bin].count, bin, total_energy_histogram[bin]);
        if ( simulation_incomplete )
            MPI_Allreduce( &( my_walker -> wl_walker.wl_histograms.histograms[bin].count ), &( total_energy_histogram[bin] ), 1, MPI_UNSIGNED, MPI_SUM, window_communicators[ comm_id ] );
        else
            MPI_Allreduce( &( my_walker -> system_obs.obs_array[ bin * convert(System_Obs_enum_t::NUM_OBS) + convert(System_Obs_enum_t::counts_per_bin) ] ), &( total_energy_histogram_obs_t[bin] ), 1, MPI_OBS_TYPE, MPI_SUM, window_communicators[ comm_id ] );
        if ( my_world_rank == 2 || my_world_rank == 3 )
            printf("\nID %d: counts[%ld] = %ld, total_counts[%ld] = %e", my_world_rank, bin, my_walker -> wl_walker.wl_histograms.histograms[bin].count, bin, total_energy_histogram_obs_t[bin]);
    }
    printf("\n");

    // Use this to compute this walker's weights
    LOGDOS_TYPE * logdos_weights = new LOGDOS_TYPE [ num_bins ] ();
    OBS_TYPE    * obs_weights    = new OBS_TYPE    [ num_bins * convert(System_Obs_enum_t::NUM_OBS) ] ();
    for ( size_t bin = 0; bin != num_bins; ++bin )
    {
        LOGDOS_TYPE weight = 0.;
        if (simulation_incomplete)
            weight = static_cast<LOGDOS_TYPE>( my_walker -> wl_walker.wl_histograms.histograms[bin].count ) / static_cast<LOGDOS_TYPE>( total_energy_histogram[bin] );
        else
        {
            weight = static_cast<LOGDOS_TYPE>( my_walker -> system_obs.obs_array[bin * convert(System_Obs_enum_t::NUM_OBS) + convert(System_Obs_enum_t::counts_per_bin)] ) / static_cast<LOGDOS_TYPE>( total_energy_histogram_obs_t[bin] );
            printf("\nID %d: my counts[%ld] = %e, total counts[%ld] = %e, weight = %e", my_world_rank, bin, 
                                                                my_walker -> system_obs.obs_array[bin * convert(System_Obs_enum_t::NUM_OBS) + convert(System_Obs_enum_t::counts_per_bin)],
                                                                bin, total_energy_histogram_obs_t[bin], weight );
        }
        
        printf("\nID = %d: bin %ld weight = %e", my_world_rank, bin, weight);
        logdos_weights[bin] = weight * ( my_walker -> wl_walker.wl_histograms.histograms[bin].logdos );
        if ( my_world_rank == 2 || my_world_rank == 3 )
            printf("\nID = %d: logdos_weights[%ld] = %e", my_world_rank, bin, logdos_weights[bin]);

        for ( size_t ob = 0; ob != convert(System_Obs_enum_t::NUM_OBS); ++ob )
        {
            obs_weights[ bin * convert(System_Obs_enum_t::NUM_OBS) + ob ] = my_walker -> system_obs.obs_array[ bin * convert(System_Obs_enum_t::NUM_OBS) + ob ];
            // Just sum the counts like normal
            if ( ob != convert(System_Obs_enum_t::counts_per_bin) )
               obs_weights[ bin * convert(System_Obs_enum_t::NUM_OBS) + ob ] *= static_cast<OBS_TYPE>(weight); 
        } 
    }

    // Break up the reductions for better memory management
    // Now use a reduction to sum the logdos and deposit it accordingly
    for ( size_t bin = 0; bin != num_bins; ++bin )
    {
        MPI_Allreduce( &( logdos_weights[bin] ), &( my_walker -> wl_walker.wl_histograms.histograms[bin].logdos ), 1, MPI_LOGDOS_TYPE, MPI_SUM, window_communicators[ comm_id ] );
    }
    // Now use a reduction to sum the observables and deposit them accordingly
    for ( size_t bin = 0; bin != num_bins; ++bin )
    {
        for ( size_t ob = 0; ob != convert(System_Obs_enum_t::NUM_OBS); ++ob )
        {
            const size_t ob_index = bin * convert(System_Obs_enum_t::NUM_OBS) + ob;
            MPI_Allreduce( &( obs_weights[ob_index] ), &( my_walker -> system_obs.obs_array[ob_index] ), 1, MPI_OBS_TYPE, MPI_SUM, window_communicators[ comm_id ] );
        }
    }

    delete [] total_energy_histogram;
    delete [] total_energy_histogram_obs_t;
    delete [] logdos_weights;
    delete [] obs_weights;
}

// Check for flatness within a window
// The window is only flat if every replica is flat.
bool REWL_simulation::check_if_window_is_flat( const int * const my_comm_ids, MPI_Comm * const window_communicators ) const
{
    int window_is_flat = false;
    int replica_is_flat = static_cast<int>( my_walker -> wl_walker.is_flat( REWL_Parameters::flatness_criterion ) );

    int comm_id = my_comm_ids[ Communicators::window_comm ];
    MPI_Barrier( window_communicators[comm_id] );

    // Reduce whether each replica is flat
    MPI_Allreduce( &replica_is_flat, &window_is_flat, 1, MPI_INT, MPI_PROD, window_communicators[ comm_id ] );
    
    return static_cast<bool>( window_is_flat );
}

// Main function for the simulation.
void REWL_simulation::simulate(
#if PRINT_HISTOGRAM
                                const std::filesystem::path & histogram_path
#endif
#ifndef INDEPENDENT_WALKERS
                                const int * const my_ids_per_comm, const int * const my_comm_ids, 
                                MPI_Comm * const local_communicators, MPI_Comm * const window_communicators
#endif
                               ) const
{
#if COLLECT_TIMINGS
    auto start = std::chrono::high_resolution_clock::now();
    auto iteration_start = start;
    auto timer = start;
#endif

    size_t iteration_counter = 1;
    size_t sweep_counter = 0;
    bool simulation_incomplete = true;
    int exchange_direction = Communicators::even_comm;
#if SAMPLE_AFTER
    bool sample_observables = false;
    int i_sample_observables = 0;
    //int we_all_sample_observables = 0;
    bool reset_incrementer = true;
#endif
    int i_am_done = 0;           // Integer to store whether this processor is finished
    constexpr size_t rewl_updates_per_check = REWL_Parameters::sweeps_per_check / REWL_Parameters::sweeps_per_exchange;
    
    while (simulation_incomplete)
    {
        for ( size_t rewl_update = 0; rewl_update != rewl_updates_per_check; ++rewl_update )
        {
            // First update the walker up until the sweeps_per_exchange
#if SAMPLE_AFTER
            my_walker -> wang_landau_walk(REWL_Parameters::sweeps_per_exchange, sample_observables); 
#else
            my_walker -> wang_landau_walk(REWL_Parameters::sweeps_per_exchange); 
#endif
            sweep_counter += REWL_Parameters::sweeps_per_exchange;
    
            // Then undergo the exchange update
#if SAMPLE_AFTER
            replica_exchange_update( exchange_direction, iteration_counter, sample_observables, my_ids_per_comm, my_comm_ids, local_communicators );
#else
            replica_exchange_update( exchange_direction, iteration_counter, my_ids_per_comm, my_comm_ids, local_communicators );
#endif
            //exchange_direction = ( exchange_direction == Communicators::even_comm ? Communicators::odd_comm : Communicators::even_comm ); 

        }
 
#if PRINT_HISTOGRAM
        if ( sweep_counter % (REWL_Parameters::sweeps_per_check) == 0 )
            my_walker -> wl_walker.wl_histograms.print_histogram_counts(iteration_counter, histogram_path);
#endif
        MPI_Barrier(MPI_COMM_WORLD);

        // Now check to see if the window is flat
        if ( check_if_window_is_flat( my_comm_ids, window_communicators ) )
        {
#if PRINT_HISTOGRAM
            my_walker -> wl_walker.wl_histograms.print_histogram_counts(iteration_counter, histogram_path);
#endif
            // If the window is flat, average the logdos and observables
            average_and_redistribute_window( simulation_incomplete, my_comm_ids, window_communicators );

            // Reset only the energy histogram and leave
            // the logdos alone.

            my_walker -> wl_walker.reset_histogram();

            printf("\nID %d: incrementer before = %e", my_world_rank, my_walker -> incrementer);
            // TODO: Generalize to 1/t algorithm.
            my_walker -> incrementer *= 0.5;

            printf("\nID %d: incrementer after = %e", my_world_rank, my_walker -> incrementer);
            
#if SAMPLE_AFTER
            if ( my_walker -> incrementer < REWL_Parameters::final_increment )
            {
                if (i_sample_observables == false) 
                {
                    i_sample_observables = 1;
                    //my_walker -> incrementer = 1.;   // Reset the incrementer and walk again
                    simulation_incomplete = true;
                }
                else
                {
#ifndef INDEPENDENT_WALKERS
                    i_am_done = 1;                  // This walker has finised complete
#else
                    simulation_incomplete = false;  // Kill the simulation if it is complete after sampling
#endif
                }
            }
            else simulation_incomplete = true;
#else
#ifndef INDEPENDENT_WALKERS
            i_am_done = ( my_walker -> incrementer >= REWL_Parameters::final_increment ? 0 : 1 );
#else
            simulation_incomplete = ( my_walker -> incrementer >= REWL_Parameters::final_increment );
#endif
#endif

#if COLLECT_TIMINGS
            timer = std::chrono::high_resolution_clock::now();
            std::chrono::duration<float> time_elapsed = timer - iteration_start;
            printf("\n\n\nID %d: Iteration %ld complete after %e seconds.", my_world_rank, iteration_counter, time_elapsed.count());
            printf("\nID %d: Total Sweeps = %ld = %e Wang Landau Updates.", my_world_rank, sweep_counter, static_cast<float>(sweep_counter * System_Parameters::N) );
            printf("\nID %d: Sweep rate  = %e sweeps per second", my_world_rank, static_cast<float>(sweep_counter)/time_elapsed.count());
            printf("\nID %d: Update rate = %e updates per second", my_world_rank, static_cast<float>(sweep_counter * System_Parameters::N)/time_elapsed.count());
            fflush(stdout);
            iteration_start = std::chrono::high_resolution_clock::now();
#endif

            ++iteration_counter;
            sweep_counter = 0;
        }

#if SAMPLE_AFTER
        // Throw up a barrier to check if we sample observables
        //MPI_Barrier(MPI_COMM_WORLD);
        //MPI_Allreduce( &i_sample_observables, &we_all_sample_observables, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD );
        
        //if ( we_all_sample_observables == 1 && reset_incrementer )
        if ( i_sample_observables == 1 && reset_incrementer )
        {
            sample_observables = true;
            my_walker -> incrementer = 1.;
            reset_incrementer = false;
        }        
#endif
        
        // Throw up a barrier to check if the simulation is over
        int completed_reduction = 0;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce( &i_am_done, &completed_reduction, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD );
        
        if ( completed_reduction == 1 )
            simulation_incomplete = false;
        else
            simulation_incomplete = true;

    }

#if COLLECT_TIMINGS
            timer = std::chrono::high_resolution_clock::now();
            std::chrono::duration<float> simulation_time = timer - start;
            printf("\n\nID %d: Total Simulation Time: %e seconds.", my_world_rank, simulation_time.count());
#endif

    // Finally average the results within a single window
    // Only the 0-processor in each window will communicate 
    // with the master processor at the end
    average_and_redistribute_window( simulation_incomplete, my_comm_ids, window_communicators ); 
}
#endif

#endif
