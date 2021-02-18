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

// This is a blocking send of a 
// vector of vectors
template<typename data_t>
void mpi_send_table( const int proc_recv, const std::vector<std::vector<data_t> > & data,
                     const MPI_Datatype datatype, const int tag, const MPI_Comm comm )
{
    const unsigned nbins = data.size();
    // First send the outer size of the table
    // for proper resizing
    MPI_Send( &nbins, 1, MPI_UNSIGNED, proc_recv, tag, comm );

    // Now iterate through the data table and
    // send each vector to the recv processor
    int inner_tag = 1;
    for ( const auto &vec : data )
    {
        mpi_send_array( proc_recv, vec.size(), vec.data(), datatype, tag + inner_tag, comm );
        ++inner_tag;
    }
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

// This is a blocking recv of a vector
// of vectors into a vector of vectors.
template<typename data_t>
void mpi_recv_table_to_table( const int proc_send,  
                              std::vector<std::vector<data_t> > & data_vec,
                              const MPI_Datatype datatype,
                              const int tag, const MPI_Comm comm, MPI_Status * status )
{  
    // Recv the outer size of the 
    // vector of vectors
    unsigned outer_size = 0;
    MPI_Recv( &outer_size, 1, MPI_UNSIGNED, proc_send, tag, comm, status );
    
    // Resize the outer data buffer
    data_vec.resize( outer_size );

    int inner_tag = 1;
    for ( auto &vec : data_vec )
    {
        MPI_Status inner_status;
        mpi_recv_array_to_vector<data_t>( proc_send, vec, datatype, tag + inner_tag, comm, &inner_status );
        ++inner_tag;
    }
}
#endif
