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
    const size_t num_temps; 
    const std::uint64_t seed;
    random_number_generator<energy_t> rng;

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
        // Add 1 to num_temps to include final_temp
        temperature_energies = new energy_t [2*(num_temps + 1)];
    }

    ~Simulated_Annealer() { delete [] temperature_energies; }

    void sa_update( const size_t idx, const energy_t temp, Hamiltonian_t * const system );
    void simulate_annealing( const size_t num_sites, const size_t num_flavors, Hamiltonian_t * const system );
};

template<typename energy_t, class Hamiltonian_t, class State_t>
void Simulated_Annealer<energy_t, Hamiltonian_t, State_t>::
     sa_update( const size_t idx, const energy_t beta, Hamiltonian_t * const system ) 
{
    State_t temporary_state;
    system -> change_state(idx, temporary_state);
    energy_t paccept = ( temporary_state.energy < system -> current_state.energy ? 1. : -1. );

    if (paccept == -1.)
        paccept = exp( beta * (system -> current_state.energy - temporary_state.energy) );

    if (rng() < paccept)
        system -> set_state(idx, temporary_state);
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
        for ( size_t sweep = 0; sweep != sweeps_per_temp; ++sweep )
        {
           for ( size_t flavor = 0; flavor != num_flavors; ++flavor )
           {
               for ( size_t site = 0; site != num_sites; ++site )
               {
                   sa_update( site, beta, system );
               }
           }
        }

        // Record the energy and the temperature
        temperature_energies[2 * Tidx + 0] = temperature; 
        temperature_energies[2 * Tidx + 1] = system -> current_state.energy; 

        printf("\nStep %ld / %ld: Energy = %e", Tidx, num_temps, system -> current_state.energy);
        printf("\nTemperature = %.4e\nBeta = %.4e", temperature, beta);
#if COLLECT_TIMINGS
        timer = std::chrono::high_resolution_clock::now();
        std::chrono::duration<float> time_elapsed = timer - iteration_start;
        printf("\nTime to Complete: %e seconds", time_elapsed.count());
        printf("\nSweeps per second: %e", sweeps_per_temp / time_elapsed.count());
        printf("\nUpdates per second: %e", num_flavors * num_sites * sweeps_per_temp / time_elapsed.count());
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
    printf("\nTotal Time Elapsed: %e seconds", time_elapsed.count());
    printf("\nAverage Sweeps per second: %e", (num_temps + 1) * sweeps_per_temp / time_elapsed.count());
    printf("\nAverage Updates per second: %e", (num_temps + 1) * num_flavors * num_sites * sweeps_per_temp / time_elapsed.count());
    printf("\n");
    fflush(stdout);
#endif
}

#endif
