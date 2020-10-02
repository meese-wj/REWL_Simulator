#ifndef MPI_HELPERS

#include <cmath>
#include <mpi.h>

// This is a blocking send of an array
template<typename data_t>
void mpi_send_array( const int proc_recv, const size_t data_size, const data_t * const data,
                     const MPI_Datatype datatype, const int tag, const MPI_Comm comm )
{
    MPI_Send( data, static_cast<int>(data_size), datatype, proc_recv, tag, comm );
}

// This is a blocking recv of an array
template<typename data_t>
void mpi_recv_array( const int proc_send, data_t *& data, const MPI_Datatype datatype,
                     const int tag, const MPI_Comm comm, MPI_Status * status )
{
    // Destroy the data buffer 
    delete [] data;
    
    // Probe for a new size
    MPI_Probe( proc_send, tag, comm, status );
    int data_size = 0;
    MPI_Get_count( status, datatype, &data_size );

    // Resize the data buffer
    data = new data_t [ data_size ];

    // Receive the buffer
    MPI_Recv( data, data_size, datatype, proc_send, tag, comm, status );
}

#endif
