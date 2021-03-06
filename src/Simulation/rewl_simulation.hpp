#ifndef REWL_SIMULATION
#define REWL_SIMULATION
/* This file contains the templated simulation
 * struct that will initialize the simulation
 * and contain all the relevant types. */
#include <string>
#include <chrono>
#include <algorithm> // For std::shuffle
#include <unistd.h> // For usleep

#include <histogram_index.hpp>
#include <glazier.hpp>
#include <rewl_walker.hpp>
#include <rewl_parameter_string.hpp>
#include <mpi_rewl_comm_setup.hpp>

#if AT_DENSITIES
#include <Ashkin_Teller2d/Density_Plots/density_averaging.cpp>
#endif

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

    OBS_TYPE * old_counts = nullptr;

    REWL_simulation( const ENERGY_TYPE all_energy_min, const ENERGY_TYPE all_energy_max, const ENERGY_TYPE all_energy_binsize );

    ~REWL_simulation()
    {
        delete window_maker;
        delete my_walker;
        delete [] old_counts;
    }

#ifndef INDEPENDENT_WALKERS
#if SAMPLE_AFTER
    size_t replica_exchange_update( int & exchange_direction, const size_t iteration_counter, const bool sample_observables, 
                                    const int * const my_ids_per_comm, const int * const my_comm_ids, MPI_Comm * const local_communicators ) const;
#else
    size_t replica_exchange_update( int & exchange_direction, const size_t iteration_counter, 
                                    const int * const my_ids_per_comm, const int * const my_comm_ids, MPI_Comm * const local_communicators ) const;
#endif
    void average_and_redistribute_window( const bool simulation_incomplete, const int * const my_ids_per_comm, const int * const my_comm_ids, MPI_Comm * const window_communicators ) const;
    bool check_if_window_is_flat( const bool oot_engaged, const int * my_comm_ids, MPI_Comm * const window_communicators ) const;
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
REWL_simulation::REWL_simulation( const ENERGY_TYPE all_energy_min, const ENERGY_TYPE all_energy_max, const ENERGY_TYPE all_energy_binsize  )
{
#if MPI_ON
    MPI_Comm_rank( MPI_COMM_WORLD, &my_world_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &REWL_Parameters::num_walkers );
    i_am_the_master = ( my_world_rank == REWL_MASTER_PROC );
#endif
    
    // Construct the glazier
    window_maker = new glazier<ENERGY_TYPE, histogram_index<ENERGY_TYPE> >
                        (all_energy_min, all_energy_max,
                         all_energy_binsize, 
                         static_cast<size_t>(REWL_Parameters::num_walkers) / REWL_Parameters::replicas_per_window,
                         REWL_Parameters::replicas_per_window, 
                         static_cast<ENERGY_TYPE>(REWL_Parameters::window_overlap));

    window_maker -> construct_windows();

    // Construct the walker
    ENERGY_TYPE walker_min = window_maker -> all_windows[my_world_rank].minimum;
    ENERGY_TYPE walker_max = window_maker -> all_windows[my_world_rank].maximum;
    ENERGY_TYPE walker_bin_size = window_maker -> all_windows[my_world_rank].bin_size;
    size_t walker_num_bins = window_maker -> all_windows[my_world_rank].num_bins;

#if DIFFERENT_SEEDS
    // Sleep for my_world_rank * 5 milliseconds
    // and then check the high resolution clock
    // for my seed. This should be good enough to
    // keep all the seeds different.
    usleep( my_world_rank * 5000 );
    std::uint64_t walker_seed = static_cast<std::uint64_t>( std::chrono::high_resolution_clock::now().time_since_epoch().count() );
#else
    std::uint64_t walker_seed = static_cast<std::uint64_t> (1);
#endif

    my_walker = new REWL_Walker<ENERGY_TYPE, LOGDOS_TYPE, OBS_TYPE, histogram_index<ENERGY_TYPE> >
                (walker_min, walker_max, walker_bin_size, walker_num_bins, walker_seed);

    old_counts = new OBS_TYPE [ walker_num_bins ] ();
}

/* ==================================================================================================================== */
/* This will not turn on replica exchange */
#if INDEPENDENT_WALKERS 

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

    // Always have this defined for readability
    // although it is only used for the 1/t
    // algorithm.
    bool oot_engaged = false;
    
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
        if ( my_walker -> wl_walker.is_flat( REWL_Parameters::flatness_criterion, oot_engaged ) )
        {
            // Reset only the energy histogram and leave
            // the logdos alone.
#if PRINT_HISTOGRAM
            my_walker -> wl_walker.wl_histograms.print_histogram_counts(iteration_counter, histogram_path);
#endif

            my_walker -> wl_walker.reset_histogram();

            printf("\nID %d: incrementer before = %e", my_world_rank, my_walker -> incrementer);
#if ONE_OVER_T_ALGORITHM
            if (!oot_engaged)
            {
                oot_engaged = (my_walker -> incrementer <= 1./static_cast<double>(sweep_counter * iteration_counter));
                if (oot_engaged)
                    printf("\nID %d: 1/t algorithm engaged at %ld sweeps.", my_world_rank, sweep_counter * iteration_counter);
            }
#endif

            // This will work generally even without 1/t
            my_walker -> incrementer = ( 0.5 * my_walker -> incrementer ) * (!oot_engaged) + 1./static_cast<double>(sweep_counter * iteration_counter) * oot_engaged;

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
            if ( (!oot_engaged || iteration_counter % REWL_Parameters::iterations_per_stdout == 0 ) )
            {
                timer = std::chrono::high_resolution_clock::now();
                std::chrono::duration<float> time_elapsed = timer - iteration_start;
                printf("\n\n\nID %d: Iteration %ld complete after %e seconds.", my_world_rank, iteration_counter, time_elapsed.count());
                printf("\nID %d: Total Sweeps = %ld = %e Wang Landau Updates.", my_world_rank, sweep_counter, static_cast<float>(sweep_counter * System_Parameters::N) );
                printf("\nID %d: Sweep rate  = %e sweeps per second", my_world_rank, static_cast<float>(sweep_counter)/time_elapsed.count());
                printf("\nID %d: Update rate = %e updates per second", my_world_rank, static_cast<float>(sweep_counter * System_Parameters::N)/time_elapsed.count());
                fflush(stdout);
            }
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
/* ==================================================================================================================== */

/* ==================================================================================================================== */
/* This will turn on replica exchange. */
#else  // INDEPENDENT_WALKERS

// Single replica exchange update
size_t REWL_simulation::replica_exchange_update( int & exchange_direction, const size_t iteration_counter, 
#if SAMPLE_AFTER
                                               const bool sample_observables,
#endif
                                               const int * const my_ids_per_comm, const int * const my_comm_ids, MPI_Comm * const local_communicators ) const
{
    size_t ret_val = 0;
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
            //int wr = 0;
            //MPI_Comm_rank( MPI_COMM_WORLD, &wr );
            //printf("\n");
            for ( size_t replica = 0; replica != REWL_Parameters::replicas_per_window; ++replica )
            {
                // Assign the replica in the lower window the shuffled_partner at replica
                partners[ replica ] = shuffled_partners[ replica ];
                // Assign the shuffled_partner at replica the replica in the lower window
                partners[ shuffled_partners[ replica ] ] = replica;
                //printf("\nwr %d: replica = %ld, partners[%ld] = %d, partners[%d] = %d", wr, replica, replica, partners[replica], shuffled_partners[replica], partners[shuffled_partners[replica]]);
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

            double pexchange = 0.;
            if ( my_walker -> energy_in_range( new_energy ) )
            {
                pexchange = static_cast<double>( exp( std::min(0., my_walker -> get_logdos( current_energy ) - my_walker -> get_logdos( new_energy )) ) );
            }
            //printf("\nID %d = %d: pexchange before = %e", my_world_rank, my_ids_per_comm[exchange_direction], pexchange);

            bool we_do_exchange = false;
            if ( my_ids_per_comm[ exchange_direction ] < static_cast<int>(REWL_Parameters::replicas_per_window) )
            {
                // Have the lower ids be the calculator
                double other_pexchange = 0.;
                MPI_Recv( &other_pexchange, 1, MPI_DOUBLE, partner_index, 2, local_communicators[ comm_id ], &status );

                // The exchange probability comes from the product of both
                // energy moves
                pexchange *= other_pexchange;

                we_do_exchange = ( pexchange != 0. && ( my_walker -> get_rand() < pexchange ) );

                // Send the whether the result is made to the partner
                MPI_Send( &we_do_exchange, 1, MPI_CXX_BOOL, partner_index, 3, local_communicators[ comm_id ] );
                MPI_Send( &pexchange, 1, MPI_DOUBLE, partner_index, 4, local_communicators[ comm_id ] );
            }
            else
            {
                // Send the exchange probability to the calculator
                // and await a response
                MPI_Send( &pexchange, 1, MPI_DOUBLE, partner_index, 2, local_communicators[ comm_id ] );
                MPI_Recv( &we_do_exchange, 1, MPI_CXX_BOOL, partner_index, 3, local_communicators[ comm_id ], &status );
                MPI_Recv( &pexchange, 1, MPI_DOUBLE, partner_index, 4, local_communicators[ comm_id ], &status );
            }
            //printf("\nID %d = %d: pexchange after = %e\n", my_world_rank, my_ids_per_comm[exchange_direction], pexchange);

            if ( we_do_exchange )
            {
                ret_val = 1;
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
    return ret_val;
}

// Get the counts per bin index
inline size_t counts_index( const size_t bin )
{
    return bin * convert(System_Obs_enum_t::NUM_OBS) + convert(System_Obs_enum_t::counts_per_bin);
}

// Get an observable index at a bin
inline size_t obs_index( const size_t bin, const size_t ob )
{
    return bin * convert( System_Obs_enum_t::NUM_OBS ) + ob;
}

// Compute the weighted average in a window and redistribute
void REWL_simulation::average_and_redistribute_window( const bool simulation_incomplete, const int * const my_ids_per_comm, const int * const my_comm_ids, MPI_Comm * const window_communicators ) const
{
    if ( REWL_Parameters::replicas_per_window == 1 )
        return;

    int comm_id = my_comm_ids[ Communicators::window_comm ];

    // Store the sum of each walker's energy histogram in an array
    size_t num_bins = my_walker -> wl_walker.wl_histograms.num_bins;
    size_t * total_energy_histogram = nullptr; 
    OBS_TYPE * total_energy_histogram_obs_t = nullptr;

    if ( simulation_incomplete )
        total_energy_histogram = new size_t [ num_bins ] ();
    else
        total_energy_histogram_obs_t = new OBS_TYPE [ num_bins ] ();

    MPI_Barrier( window_communicators[ comm_id ] );

    for ( size_t bin = 0; bin != num_bins; ++bin )
    {
        if ( simulation_incomplete )
            MPI_Allreduce( &( my_walker -> wl_walker.wl_histograms.histograms[bin].count ), &( total_energy_histogram[bin] ), 1, MPI_UNSIGNED, MPI_SUM, window_communicators[ comm_id ] );
        else
            MPI_Allreduce( &( my_walker -> system_obs.obs_array[ counts_index(bin) ] ), &( total_energy_histogram_obs_t[bin] ), 1, MPI_OBS_TYPE, MPI_SUM, window_communicators[ comm_id ] ); 
    }

    // Use this to compute this walker's weights
    LOGDOS_TYPE * logdos_weights = new LOGDOS_TYPE [ num_bins ] ();
    OBS_TYPE    * obs_weights    = new OBS_TYPE    [ num_bins * convert(System_Obs_enum_t::NUM_OBS) ] ();
    for ( size_t bin = 0; bin != num_bins; ++bin )
    {
        LOGDOS_TYPE weight = 0.;
        
        if (simulation_incomplete)
            weight = static_cast<LOGDOS_TYPE>( my_walker -> wl_walker.wl_histograms.histograms[bin].count ) / static_cast<LOGDOS_TYPE>( total_energy_histogram[bin] );
        else
            weight = static_cast<LOGDOS_TYPE>( my_walker -> system_obs.obs_array[ counts_index(bin) ] ) / static_cast<LOGDOS_TYPE>( total_energy_histogram_obs_t[bin] );
        

        //weight = 1. / static_cast<LOGDOS_TYPE>(REWL_Parameters::replicas_per_window);
         
        logdos_weights[bin] = weight * ( my_walker -> wl_walker.wl_histograms.histograms[bin].logdos );
        
        if ( !simulation_incomplete )
        {
            for ( size_t ob = 0; ob != convert(System_Obs_enum_t::NUM_OBS); ++ob )
            {
                obs_weights[ obs_index(bin, ob) ] = my_walker -> system_obs.obs_array[ obs_index(bin, ob) ];
                // Sum the weighted observables in the next step
                if ( ob != convert(System_Obs_enum_t::counts_per_bin ) )
                   obs_weights[ obs_index(bin, ob) ] *= static_cast<OBS_TYPE>(weight); 
                // Just sum the counts like normal
                else
                {
                    // Subtract out the old counts so only the new counts 
                    // from the last iteration are counted
                    obs_weights[ counts_index(bin) ] -= old_counts[bin]; 
                }
            } 
        } 
    }

    // Break up the reductions for better memory management
    // Now use a reduction to sum the logdos and deposit it accordingly
    for ( size_t bin = 0; bin != num_bins; ++bin )
    {
        MPI_Allreduce( &( logdos_weights[bin] ), &( my_walker -> wl_walker.wl_histograms.histograms[bin].logdos ), 1, MPI_LOGDOS_TYPE, MPI_SUM, window_communicators[ comm_id ] );
    }
    // Now use a reduction to sum the observables and deposit them accordingly
    if (!simulation_incomplete)
    {
        for ( size_t bin = 0; bin != num_bins; ++bin )
        {
            for ( size_t ob = 0; ob != convert(System_Obs_enum_t::NUM_OBS); ++ob )
            {
                const size_t ob_index = obs_index(bin, ob);
                MPI_Allreduce( &( obs_weights[ob_index] ), &( my_walker -> system_obs.obs_array[ob_index] ), 1, MPI_OBS_TYPE, MPI_SUM, window_communicators[ comm_id ] );
                if ( ob == convert(System_Obs_enum_t::counts_per_bin) )
                {
                    // Add back the old counts
                    my_walker -> system_obs.obs_array[ob_index] += old_counts[bin];
                    // Record the current counts in the old counts now
                    old_counts[bin] = my_walker -> system_obs.obs_array[ob_index];
                }
            }
        }
    }

    delete [] total_energy_histogram;
    delete [] total_energy_histogram_obs_t;
    delete [] logdos_weights;
    delete [] obs_weights;
}

// Check for flatness within a window
// The window is only flat if every replica is flat.
bool REWL_simulation::check_if_window_is_flat( const bool oot_engaged, const int * const my_comm_ids, MPI_Comm * const window_communicators ) const
{
    int window_is_flat = false;
    int replica_is_flat = static_cast<int>( my_walker -> wl_walker.is_flat( REWL_Parameters::flatness_criterion, oot_engaged ) );

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
    size_t exchange_counter = 0;
    bool simulation_incomplete = true;
    int exchange_direction = Communicators::even_comm;
#if SAMPLE_AFTER
    bool sample_observables = false;
    int i_sample_observables = 0;
    //int we_all_sample_observables = 0;
    bool reset_incrementer = true;
#endif
    int i_am_done = 0;                         // Integer to store whether this processor is finished
    int * individuals_done = nullptr;          // Array to house which processes are done before broadcasting
    if ( my_world_rank == REWL_MASTER_PROC )  
        individuals_done = new int [REWL_Parameters::num_walkers] ();
    bool print_slowpokes = false;

    constexpr size_t rewl_updates_per_check = REWL_Parameters::sweeps_per_check / REWL_Parameters::sweeps_per_exchange;
    
    // Always have this defined for readability
    // although it is only used for the 1/t
    // algorithm.
    bool oot_engaged = false;

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
            exchange_counter += replica_exchange_update( exchange_direction, iteration_counter, sample_observables, my_ids_per_comm, my_comm_ids, local_communicators );
#else
            exchange_counter += replica_exchange_update( exchange_direction, iteration_counter, my_ids_per_comm, my_comm_ids, local_communicators );
#endif
            //exchange_direction = ( exchange_direction == Communicators::even_comm ? Communicators::odd_comm : Communicators::even_comm ); 

        }
 
#if PRINT_HISTOGRAM
        if ( sweep_counter % (REWL_Parameters::sweeps_per_check) == 0 )
            my_walker -> wl_walker.wl_histograms.print_histogram_counts(iteration_counter, histogram_path);
#endif
        MPI_Barrier(MPI_COMM_WORLD);

        // Now check to see if the window is flat
        if ( check_if_window_is_flat( oot_engaged, my_comm_ids, window_communicators ) )
        {
#if PRINT_HISTOGRAM
            my_walker -> wl_walker.wl_histograms.print_histogram_counts(iteration_counter, histogram_path);
#endif
            // If the window is flat, average the logdos and observables 
            average_and_redistribute_window( simulation_incomplete, my_ids_per_comm, my_comm_ids, window_communicators );
            
            // Reset only the energy histogram and leave
            // the logdos alone.
            my_walker -> wl_walker.reset_histogram();
            
            if ( my_world_rank == REWL_MASTER_PROC && (!oot_engaged || iteration_counter % REWL_Parameters::iterations_per_stdout == 0 ) )
                printf("\nID %d: incrementer before = %e", my_world_rank, my_walker -> incrementer);

#if ONE_OVER_T_ALGORITHM
            if (!oot_engaged)
            {
                oot_engaged = (my_walker -> incrementer <= 1./static_cast<double>(sweep_counter * iteration_counter));
                if (oot_engaged)
                    printf("\nID %d: 1/t algorithm engaged at %ld sweeps.", my_world_rank, sweep_counter * iteration_counter);
            }
#endif

            // This will work generally even without 1/t
            my_walker -> incrementer = ( 0.5 * my_walker -> incrementer ) * (!oot_engaged) + 1./static_cast<double>(iteration_counter * sweep_counter) * oot_engaged;


            if ( my_world_rank == REWL_MASTER_PROC && (!oot_engaged || iteration_counter % REWL_Parameters::iterations_per_stdout == 0 ) )
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
                    i_am_done = 1;                  // This walker has completed its run
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
            print_slowpokes = (!oot_engaged || iteration_counter % REWL_Parameters::iterations_per_stdout == 0 );
            timer = std::chrono::high_resolution_clock::now();
            std::chrono::duration<float> time_elapsed = timer - iteration_start;
            if ( my_world_rank == REWL_MASTER_PROC && ( !oot_engaged || iteration_counter % REWL_Parameters::iterations_per_stdout == 0 ) )
            {
                printf("\n\n\nID %d: Iteration %ld complete after %e seconds.", my_world_rank, iteration_counter, time_elapsed.count());
                printf("\nID %d: Total Sweeps = %ld = %e Wang Landau Updates.", my_world_rank, sweep_counter, static_cast<float>(sweep_counter * System_Parameters::N) );
                printf("\nID %d: Sweep rate  = %e sweeps per second", my_world_rank, static_cast<float>(sweep_counter)/time_elapsed.count());
                printf("\nID %d: Update rate = %e updates per second", my_world_rank, static_cast<float>(sweep_counter * System_Parameters::N)/time_elapsed.count());
                printf("\nID %d: Successful exchanges %ld", my_world_rank, exchange_counter);
                fflush(stdout);
            }
            iteration_start = std::chrono::high_resolution_clock::now();
#endif

            ++iteration_counter;
            sweep_counter = 0;
            exchange_counter = 0;
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
        int completed_reduction = 1;
        MPI_Barrier(MPI_COMM_WORLD);
        //MPI_Allreduce( &i_am_done, &completed_reduction, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD );

        // Gather all the i am dones from everyone
        MPI_Gather( &i_am_done, 1, MPI_INT, individuals_done, 1, MPI_INT, REWL_MASTER_PROC, MPI_COMM_WORLD );

        if ( my_world_rank == REWL_MASTER_PROC )
        {
            if ( print_slowpokes )
            {
                printf("\n\nWaiting on walkers:\n\t");
                for ( int wdx = 0; wdx != REWL_Parameters::num_walkers; ++wdx )
                {
                    completed_reduction *= individuals_done[wdx];
                    if ( individuals_done[wdx] == 0 )
                        printf("%d%s", wdx, wdx == REWL_Parameters::num_walkers - 1 ? "\n" : ", ");
                }
                printf("\n");
                print_slowpokes = false;
            }
            else
            {
                for ( int wdx = 0; wdx != REWL_Parameters::num_walkers; ++wdx )
                    completed_reduction *= individuals_done[wdx];
            }
        }

        // Distribute the result
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast( &completed_reduction, 1, MPI_INT, REWL_MASTER_PROC, MPI_COMM_WORLD );
 
        if ( completed_reduction == 1 )
            simulation_incomplete = false;
        else
            simulation_incomplete = true;

    }

    delete [] individuals_done;

#if COLLECT_TIMINGS
            timer = std::chrono::high_resolution_clock::now();
            std::chrono::duration<float> simulation_time = timer - start;
            if ( my_world_rank == REWL_MASTER_PROC )
                printf("\n\nID %d: Total Simulation Time: %e seconds.", my_world_rank, simulation_time.count());
#endif

    // Finally average the results within a single window
    // Only the 0-processor in each window will communicate 
    // with the master processor at the end
    average_and_redistribute_window( simulation_incomplete, my_ids_per_comm, my_comm_ids, window_communicators );
#if AT_DENSITIES
    average_density_in_window<Observables_t<OBS_TYPE> >( my_ids_per_comm, &( my_walker -> system_obs ), my_comm_ids, window_communicators );
#endif
}
#endif // INDEPENDENT_WALKERS

#endif // REWL_SIMULATION
/* ==================================================================================================================== */
