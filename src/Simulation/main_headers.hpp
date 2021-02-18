#ifndef MAIN_HEADERS
#define MAIN_HEADERS
/* Include this header in main. */

#include <stdio.h>
#include <string>

#if MPI_ON
#include <mpi_rewl_helpers.hpp>
#include <concatenate_rewl_tables.hpp>
#ifndef INDEPENDENT_WALKERS
#include <mpi_rewl_comm_setup.hpp>
#endif
#endif

const std::string DELIMITER = "  ";

#include "rewl_simulation.hpp"

#include <file_manager.hpp>
#include <file_header.hpp>
#include <array_shift.hpp>
#include <write_microcanonical_observables.hpp>
#include <thermodynamics.hpp>
#include <write_self_averaged_observables.hpp>
#include <write_nonlinear_observables.hpp>

#if SIMULATED_ANNEALING
#include <Simulated_Annealer/simulated_annealer.hpp>
#endif

#if AT_DENSITIES
#include <Ashkin_Teller2d/Density_Plots/ashkin_teller_densities_parameters.cxx>
#include <Ashkin_Teller2d/Density_Plots/write_density_plots.cpp>
#endif



#endif
