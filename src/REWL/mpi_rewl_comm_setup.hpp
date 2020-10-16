#ifndef MPI_REWL_COMM_SETUP
#define MPI_REWL_COMM_SETUP
/* Use this file to setup communicators
 * for the Replica Exchange. */

#include <cmath>
#include <mpi.h>

enum Communicators
{
    NONE = -1, window_comm, odd_comm, even_comm, NUM_COMMS
};

int operator + ( const int idx, const Communicators enum_value )
{
    return idx + static_cast<int>(enum_value);
}

size_t operator + ( const size_t idx, const Communicators enum_value )
{
    return idx + static_cast<size_t>(enum_value);
}

// Get the communicator id for each ID.
// These will then go into an array
int get_local_comm_ID( const int my_world_rank, const int which_comm, const int num_procs, const int replicas_per_window )
{
    if ( which_comm == Communicators::window_comm )
        return my_world_rank / replicas_per_window;           // Return the window index

    int id = INT32_MAX;

    // Walker is in the lowest energy window
    if ( my_world_rank < replicas_per_window )
    {
        if      ( which_comm == Communicators::even_comm )
            id = 0;                                           // The lowest energy window is even
        else if ( which_comm == Communicators::odd_comm  )
            id = Communicators::NONE;                         // The lowest energy window has to odd window to the left
    }
    // Walker is in the highest window
    else if ( my_world_rank >= num_procs - replicas_per_window )
    {
        if ( ( num_procs / replicas_per_window ) % 2 == 0 )
        {
            // The number of processors is even so the last window
            // is odd and does not have an odd comm ID
            if      ( which_comm == Communicators::even_comm )
                id = 2 * ( my_world_rank / ( 2 * replicas_per_window ) );  // Get the id of the nearest even id from below
            else if ( which_comm == Communicators::odd_comm )
                id = Communicators::NONE;
        }
        else
        {
            // The number of processors is odd and so the last
            // window is even and does not have an even ID
            if      ( which_comm == Communicators::even_comm )
                id = Communicators::NONE;
            else if ( which_comm == Communicators::odd_comm )
                id = 2 * ( (my_world_rank - replicas_per_window) / (2 * replicas_per_window) ) + 1;
        }
    }
    // Walker is somewhere in between the second lowest and second highest window
    else
    {
        if      ( which_comm == Communicators::even_comm )
        {
            id = 2 * ( my_world_rank / ( 2 * replicas_per_window ) );
        }
        else if ( which_comm == Communicators::odd_comm )
        {
            id = 2 * ( (my_world_rank - replicas_per_window) / ( 2 * replicas_per_window ) ) + 1;
        }
    }

    return id;
}

// Set up the intra-window groups and communicators
void define_window_communicators( const int num_procs, const int replicas_per_window, const MPI_Group & World,
                                  MPI_Group * window_groups, MPI_Comm * window_communicators )
{
    int * ranks = new int [ replicas_per_window ];
    const int num_window_comms = num_procs / replicas_per_window;

    for ( int comm = 0; comm != num_window_comms; ++comm )
    {
        // Fill the ranks array with all the replicas in the energy window
        for ( int replica = 0; replica != replicas_per_window; ++replica )
            ranks[ replica ] = comm * replicas_per_window + replica;

        // Local group instantiation
        MPI_Group_incl( World, replicas_per_window, ranks, &window_groups[ comm ] );
        MPI_Comm_create( MPI_COMM_WORLD, window_groups[ comm ], &window_communicators[ comm ] );
    }
    delete [] ranks;
}

// Set up the inter-window groups and communicators
void create_local_groups_and_communicators( const int num_procs, const int replicas_per_window, const MPI_Group & World,
                                            MPI_Group * local_group, MPI_Comm * local_communicators )
{
    int * ranks = new int [ 2 * replicas_per_window ];
    int num_local_comms = ( num_procs / replicas_per_window ) - 1;

    for ( int comm = 0; comm != num_local_comms; ++comm )
    {
        // For each processor in a local communicator,
        // calculate its global rank in the MPI_COMM_WORLD
        // and cycle through two energy windows to establish
        // the even/odd communicators
        for ( int replica = 0; replica != 2 * replicas_per_window; ++replica )
            ranks[ replica ] = comm * replicas_per_window + replica;

        // Local group instantiation
        MPI_Group_incl( World,  2 * replicas_per_window, ranks, &local_group[ comm ] );
        // Local communicator instantiation
        MPI_Comm_create( MPI_COMM_WORLD, local_group[ comm ], &local_communicators[ comm ] );
    }

    delete [] ranks;
}

// Set up the IDs of a single walker within the local communicators
void determine_my_local_IDs( const int my_world_rank, const int num_procs, const int replicas_per_window,
                             int * const my_ids_per_comm, const MPI_Comm * const local_communicators, const MPI_Comm * const window_communicators )
{
    int comm_id = INT32_MAX;

    // First find the window communicator id
    comm_id = get_local_comm_ID( my_world_rank, Communicators::window_comm, num_procs, replicas_per_window );
    MPI_Comm_rank( window_communicators[ comm_id ], &my_ids_per_comm[ Communicators::window_comm ] );
   
    // Now populate the even/odd communicators
    // TODO: A picture will probably help...
    //
    // my_ids_per_comm[ Communicators::even_comm ] : even communicaotr ID defined for all walkers outside
    //                                            the highest energy window if the window number is odd
    // my_ids_per_comm[ Communicators::odd_comm  ] : odd communicator ID defined for all walkers outside 
    //                                            the highest energy window if the window number is even
    //                                            and for the lowest energy window

    for ( int which_comm = Communicators::window_comm + 1; which_comm != Communicators::NUM_COMMS; ++which_comm )
    {
        comm_id = get_local_comm_ID( my_world_rank, which_comm, num_procs, replicas_per_window );
        if ( comm_id != Communicators::NONE )
            MPI_Comm_rank( local_communicators[ comm_id ], &my_ids_per_comm[ which_comm ] );
        else
            my_ids_per_comm[ which_comm ] = Communicators::NONE;
    }
}

// Set up the communicator IDs for a single walker
void determine_my_communicators( const int my_world_rank, const int num_procs, const int replicas_per_window, int * my_comms_ids )
{
    for ( int which_comm = 0; which_comm != Communicators::NUM_COMMS; ++which_comm )
        my_comms_ids[ which_comm ] = get_local_comm_ID( my_world_rank, which_comm, num_procs, replicas_per_window );
}

#endif
