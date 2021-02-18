#ifndef AT_DENSITY_AVERAGING
#define AT_DENSITY_AVERAGING

#include <mpi.h>
#include <mpi_rewl_comm_setup.hpp>
#include "ashkin_teller_densities_parameters.cxx"
#include "ashkin_teller_densities.hpp"

template<typename Observables_t>
void average_density_in_window( const int * const my_ids_per_comm, 
                                Observables_t * const observables,
                                const int * const my_comm_ids,
                                MPI_Comm * const window_communicators )
{
    const density_int nbins = static_cast<density_int>( observables -> num_bins );
    transfer_and_normalize( nbins, observables -> density_histograms, observables -> density_float_data );
    if ( REWL_Parameters::replicas_per_window == 1 )
        return;

    int comm_id = my_comm_ids[ Communicators::window_comm ];
   
    density_float multiplier = 1./static_cast<density_float>(REWL_Parameters::replicas_per_window);
    for ( density_int idx = 0; idx != nbins * AT_Density_Parameters::total_bins; ++idx )
    {
        MPI_Allreduce( MPI_IN_PLACE, &( observables -> density_float_data[ idx ] ), 1, MPI_DOUBLE, MPI_SUM, window_communicators[ comm_id ] );
        observables -> density_float_data[idx] *= multiplier;
    }
}

#endif
