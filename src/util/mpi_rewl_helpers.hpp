#ifndef MPI_HELPERS
#define MPI_HELPERS

#include <cmath>
#include <vector>
#include <mpi.h>

/* *************************************************************************************** */
/* Array send/recv helpers                                                                 */
/* *************************************************************************************** */

// This is a blocking send of an array
template<typename data_t>
void mpi_send_array( const int proc_recv, const size_t data_size, const data_t * const data,
                     const MPI_Datatype datatype, const int tag, const MPI_Comm comm )
{
    MPI_Send( data, static_cast<int>(data_size), datatype, proc_recv, tag, comm );
}

// This is a blocking recv of an array
template<typename data_t>
void mpi_recv_array( const int proc_send, int * data_size, data_t *& data, const MPI_Datatype datatype,
                     const int tag, const MPI_Comm comm, MPI_Status * status )
{
    // Destroy the data buffer 
    delete [] data;
    
    // Probe for a new size
    MPI_Probe( proc_send, tag, comm, status );
    MPI_Get_count( status, datatype, data_size );

    // Resize the data buffer
    data = new data_t [ *data_size ];

    // Receive the buffer
    MPI_Recv( data, *data_size, datatype, proc_send, tag, comm, status );
}

// This is a blocking recv of an array 
// and conversion to a vector.
template<typename data_t>
void mpi_recv_array_to_vector( const int proc_send, std::vector<data_t> & data_vec, const MPI_Datatype datatype,
                               const int tag, const MPI_Comm comm, MPI_Status * status )
{  
    // Probe for a new size
    MPI_Probe( proc_send, tag, comm, status );
    int data_size = 0;
    MPI_Get_count( status, datatype, &data_size );

    // Resize the data buffer
    data_vec.resize( data_size );

    // Receive the buffer
    MPI_Recv( data_vec.data(), data_size, datatype, proc_send, tag, comm, status );
}

#endif
