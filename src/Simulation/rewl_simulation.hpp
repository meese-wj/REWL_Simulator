#ifndef REWL_SIMULATION
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
    
    while (simulation_incomplete)
    {
        // First update the walker up until the 
        // sweeps_per_check
        my_walker -> wang_landau_walk(REWL_Parameters::sweeps_per_check); 
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

            printf("\nincrementer before = %e", my_walker -> incrementer);
            // TODO: Generalize to 1/t algorithm.
            my_walker -> incrementer *= 0.5;

            printf("\nincrementer after = %e", my_walker -> incrementer);

            simulation_incomplete = ( my_walker -> incrementer >= REWL_Parameters::final_increment );

#if COLLECT_TIMINGS
            timer = std::chrono::high_resolution_clock::now();
            std::chrono::duration<float> time_elapsed = timer - iteration_start;
            printf("\n\n\nIteration %ld complete after %e seconds.", iteration_counter, time_elapsed.count());
            printf("\nTotal Sweeps = %ld = %e Wang Landau Updates.", sweep_counter, static_cast<float>(sweep_counter * System_Parameters::N) );
            printf("\nSweep rate  = %e sweeps per second", static_cast<float>(sweep_counter)/time_elapsed.count());
            printf("\nUpdate rate = %e updates per second", static_cast<float>(sweep_counter * System_Parameters::N)/time_elapsed.count());
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
            printf("\n\nTotal Simulation Time: %e seconds.", simulation_time.count());
#endif
}

#endif
