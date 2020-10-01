#ifndef MPI_HELPERS

#include <mpi.h>

constexpr int REWL_MASTER_PROC = 0;

// This is a blocking send of an array
void mpi_send_float_array( int proc_send, int proc_recv, const float * const data );

#endif
