#ifndef SIMULATED_ANNEALING_PARAMETERS
#define SIMULATED_ANNEALING_PARAMETERS

#include <cstdint>
using SA_energy_t = double;

namespace SA_Parameters
{
    constexpr SA_energy_t initial_temperature = 25.;
    constexpr SA_energy_t final_temperature   = 0.25;
    constexpr uint32_t num_iterations = 250;

    constexpr uint32_t block_size = 1000;                   // Number of sweeps between equilibration
                                                            // checks

    constexpr SA_energy_t energy_stdev_tolerance = 0.01;    // Allowable percent difference for energy
                                                            // variance between blocks

    
    constexpr SA_energy_t frozen_variance = 0.;             // Defined frozen condition on energy variance
}

#endif
