#ifndef SIMULATED_ANNEALER
#define SIMULATED_ANNEALER
/* Use this file to define a 
 * simulated annealing algorithm
 * to be used to find the system's 
 * ground state and ground state
 * energy.  */

#include <cmath>
#include <chrono>
#include <random_number_generators.hpp>
#include <stdio.h>
#include <limits>

#include "simulated_annealing_parameters.cxx"
#include "simulated_annealing_utilities.hpp"

template<typename energy_t>
inline energy_t calculate_temperature( const energy_t Ti, const energy_t Tf, const size_t Tidx, const size_t numT )
{
    return Ti * pow( Tf / Ti, static_cast<energy_t>(Tidx) / static_cast<energy_t>(numT) );
}

template<typename energy_t, class Hamiltonian_t, class State_t> 
struct Simulated_Annealer
{
    const size_t sweeps_per_temp;
    const energy_t initial_temp;
    const energy_t final_temp;
    energy_t equil_mean_energy = std::numeric_limits<energy_t>::max();
    energy_t equil_stdev_energy;
    energy_t min_energy_found;
    const size_t num_temps; 
    const std::uint64_t seed;
    random_number_generator<energy_t> rng;

    energy_t * equil_energies = nullptr;        // Use this to keep track of energies per block
    energy_t * temperature_energies = nullptr;

    Simulated_Annealer( const size_t _spt, 
                        const energy_t _iT, 
                        const energy_t _fT, 
                        const size_t _nT )
                        : sweeps_per_temp(_spt), 
                          initial_temp(_iT), 
                          final_temp(_fT), 
                          num_temps(_nT), 
                          seed(static_cast<std::uint64_t>( std::chrono::high_resolution_clock::now().time_since_epoch().count() )), rng(seed)
    {
        equil_energies = new energy_t [SA_Parameters::block_size];

        // Add 1 to num_temps to include final_temp
        temperature_energies = new energy_t [2*(num_temps + 1)];
    }

    ~Simulated_Annealer()
    { 
        delete [] equil_energies;
        delete [] temperature_energies;
    }

    size_t sa_update( const size_t idx, const energy_t temp, Hamiltonian_t * const system );
    void simulate_annealing( const size_t num_sites, const size_t num_flavors, Hamiltonian_t * const system );
};

template<typename energy_t, class Hamiltonian_t, class State_t>
size_t Simulated_Annealer<energy_t, Hamiltonian_t, State_t>::
       sa_update( const size_t idx, const energy_t beta, Hamiltonian_t * const system ) 
{
    State_t temporary_state;
    system -> change_state(idx, temporary_state);
    energy_t paccept = ( temporary_state.energy < system -> current_state.energy ? 1. : -1. );

    if (paccept == -1.)
        paccept = exp( beta * (system -> current_state.energy - temporary_state.energy) );

    if (rng() < paccept)
    {
        system -> set_state(idx, temporary_state);
        return 1;     
    }

    return 0;
}

// Simulated annealing function. Takes in a Hamiltonian object
// and finds its ground state configuration. Or it at least finds
// a good approximation for it.
template<typename energy_t, class Hamiltonian_t, class State_t>
void Simulated_Annealer<energy_t, Hamiltonian_t, State_t>::
     simulate_annealing( const size_t num_sites, const size_t num_flavors, Hamiltonian_t * const system ) 
{
    // First randomize the degrees of freedom
    system -> randomize_dofs();

    size_t total_sweeps = 0;
    size_t current_acceptances = 0;

    printf("\n");
#if COLLECT_TIMINGS
    auto start = std::chrono::high_resolution_clock::now();
    auto iteration_start = start;
    auto timer = start;
#endif

    for ( size_t Tidx = 0; Tidx <= num_temps; ++Tidx )
    {
        energy_t temperature = calculate_temperature(initial_temp, final_temp, Tidx, num_temps);
        energy_t beta = 1./temperature;
        size_t sweep = 0;
        
        equil_mean_energy = std::numeric_limits<energy_t>::max();
        equil_stdev_energy = 1.;
        bool sa_incomplete = true;
        while ( sa_incomplete )
        {
            for ( size_t block_idx = 0; block_idx != SA_Parameters::block_size; ++block_idx )
            {
                for ( size_t flavor = 0; flavor != num_flavors; ++flavor )
                {
                    for ( size_t site = 0; site != num_sites; ++site )
                    {
                        current_acceptances += sa_update( site, beta, system );
                        min_energy_found = ( system -> current_state.energy ) * ( system -> current_state.energy < min_energy_found )
                                           + min_energy_found * ( system -> current_state.energy >= min_energy_found );
                    }
                }
                equil_energies[block_idx] = system -> current_state.energy;
            }
            energy_t block_avg = 0.;
            energy_t block_stdev = 0.;
            array_statistics<energy_t>( SA_Parameters::block_size, equil_energies, &block_avg, &block_stdev );

            /*  First check if the new block average energy is within the previous
             *  average +/- standard deviation.
             *
             *  Then, make sure the energy variance has not changed significantly
             *  (otherwise the specific heat has a large error).
             *
             *  Finally, check if the system is frozen so as to not waste time.
             *
             *  If these three conditions are met, then this SA iteration is over.
             */
            sa_incomplete = not_within_error<energy_t>( block_avg, equil_mean_energy, equil_stdev_energy );
            sa_incomplete = sa_incomplete || not_within_tolerance<energy_t>( block_stdev * block_stdev, 
                                                                             equil_stdev_energy * equil_stdev_energy,
                                                                             SA_Parameters::energy_stdev_tolerance ); 

            sa_incomplete = not_frozen<energy_t>( block_stdev ) && sa_incomplete;
            
            equil_mean_energy  = block_avg;
            equil_stdev_energy = block_stdev;

            sweep += SA_Parameters::block_size;
        }

        total_sweeps += sweep;

        // Record the energy and the temperature
        temperature_energies[2 * Tidx + 0] = temperature; 
        temperature_energies[2 * Tidx + 1] = system -> current_state.energy; 

        printf("\nStep %ld / %ld: sweeps = %ld", Tidx, num_temps, sweep); 
        printf("\nEnergy = %e", system -> current_state.energy);
        printf("\nMin Energy Found = %e", min_energy_found);
        printf("\n|Estimate - Min| / |Estimate| = %e", (system -> current_state.energy - min_energy_found) / system -> current_state.energy );
        printf("\nTemperature = %.4e\nBeta = %.4e", temperature, beta);
#if COLLECT_TIMINGS
        timer = std::chrono::high_resolution_clock::now();
        std::chrono::duration<float> time_elapsed = timer - iteration_start;
        printf("\nTime to Complete: %e seconds", time_elapsed.count());
        printf("\nSweeps per second: %e", (sweep + 1) / time_elapsed.count());
        printf("\nUpdates per second: %e", num_flavors * num_sites * (sweep + 1) / time_elapsed.count());
        iteration_start = std::chrono::high_resolution_clock::now();
#endif
        printf("\n");
        system -> print_lattice();
        fflush(stdout);
    }

#if COLLECT_TIMINGS
    timer = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> time_elapsed = timer - start;
    printf("\n\nSimulated Annealing complete.");
    printf("\nGround State Estimate = %e", system -> current_state.energy);
    printf("\nMin Energy Found = %e", min_energy_found);
    printf("\n|Estimate - Min| / |Estimate| = %e", (system -> current_state.energy - min_energy_found) / system -> current_state.energy );
    printf("\nTotal Time Elapsed: %e seconds", time_elapsed.count());
    printf("\nAverage Sweeps per second: %e", (num_temps + 1) * total_sweeps / time_elapsed.count());
    printf("\nAverage Updates per second: %e", (num_temps + 1) * num_flavors * num_sites * total_sweeps / time_elapsed.count());
    printf("\n");
    fflush(stdout);
#endif
}

#endif
