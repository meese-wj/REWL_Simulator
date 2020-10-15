#ifndef REWL_SIMULATION
#define REWL_SIMULATION
/* This file contains the templated simulation
 * struct that will initialize the simulation
 * and contain all the relevant types. */
#include <string>
#include <chrono>

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

    void simulate(
#if PRINT_HISTOGRAM
                  const std::filesystem::path & histogram_path
#endif
#ifndef INDEPENDENT_WALKERS
                  , const int * const my_ids_per_comm, const int * const my_comm_ids, const MPI_Comm * const local_communicators
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
                         static_cast<size_t>(REWL_Parameters::num_walkers),
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

// Main function for the simulation.
void REWL_simulation::simulate(
#if PRINT_HISTOGRAM
                                const std::filesystem::path & histogram_path
#endif
#ifndef INDEPENDENT_WALKERS
                                , const int * const my_ids_per_comm, const int * const my_comm_ids, const MPI_Comm * const local_communicators
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

        // TODO: After a set number of sweeps, undergo replica exchange.
        int comm_id = my_comm_ids[ exchange_direction ];
        if ( comm_id != Communicators::NONE )
        {
            int partners [ 2 * REWL_Parameters::replicas_per_window ];

            // Only have the communicator master assign partners
            if ( my_ids_per_comm[ exchange_direction ] == 0 )
            {
                // Get the partner indices
                partners[0] = REWL_Parameters::replicas_per_window;
                partners[2 * REWL_Parameters::replicas_per_window - 1] = 0;
            }
             
            // Scatter the partner indices to everyone
            int partner_index = Communicators::NONE; 
            MPI_Scatter( partners, 1, MPI_INT, &partner_index, MPI_INT, 0, local_communicators[ comm_id ] );

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
                    pexchange = static_cast<float>( exp( my_walker -> get_logdos( current_energy ) - my_walker -> get_logdos( new_energy ) ) );
                }

                bool we_do_exchange = false;
                if ( my_ids_per_comm[ exchange_direction ] < replicas_per_window )
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
                }
                else
                {
                    // Send the exchange probability to the calculator
                    // and await a response
                    MPI_Send( &pexchange, 1, MPI_FLOAT, partner_index, 2, local_communicators[ comm_id ] );
                    MPI_Recv( &we_do_exchange, 1, MPI_CXX_BOOL, partner_index, 3, local_communicators[ comm_id ], &status );
                }

                if ( we_do_exchange )
                {
                    // Perform a MPI_Sendrecv_replace on the state and the degrees of freedom
                    mpi_exchange_state<ENERGY_TYPE, OBS_TYPE>( my_walker -> current_state(), partner_index, comm_id, local_communicators, status );
                    mpi_exchange_DoFs<OBS_TYPE>( my_walker -> DoFs(), System_Parameters::num_dof, partner_index, comm_id, local_communicators, status );
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


#if PRINT_HISTOGRAM
        if ( sweep_counter % (REWL_Parameters::sweeps_per_check) == 0 )
            my_walker -> wl_walker.wl_histograms.print_histogram_counts(iteration_counter, histogram_path);
#endif

        // Now check to see if the histogram is flat
        // TODO: Generalize this for multiple walkers per window.
        if ( my_walker -> wl_walker.is_flat( REWL_Parameters::flatness_criterion ) )
        {
            // Reset only the energy histogram and leave
            // the logdos alone.
#if PRINT_HISTOGRAM
            my_walker -> wl_walker.wl_histograms.print_histogram_counts(iteration_counter, histogram_path);
#endif

            // TODO: Generalize for multiple walkers per window. Specifically average data here.

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
#endif

#endif
