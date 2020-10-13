#ifndef MAIN_HEADERS
#define MAIN_HEADERS
/* Include this header in main. */

#include <stdio.h>
#include <string>

#if MPI_ON
#include <mpi_rewl_helpers.hpp>
#include <concatenate_rewl_tables.hpp>
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




#endif
