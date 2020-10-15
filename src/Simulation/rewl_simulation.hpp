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

        // TODO: After a set number of sweeps, undergo replica exchange.

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
