#ifndef GRID_SETUP

void define_2d_square_periodic_neighbors(const size_t Lx,
                                         const size_t Ly,
                                         const size_t num_neighbors,
                                         size_t * neighbor_array)
{
    neighbor_array = nullptr;
    neighbor_array = new size_t [num_neighbors * Lx * Ly];

    size_t neighbor_idx;
    size_t neighbor_jdx;
    for ( size_t idx = 0; idx != Lx; ++idx )
    {
        for ( size_t jdx = 0; jdx != Ly; ++jdx )
        {
            // Neighbor 0: idx - 1, jdx
            // Neighbor 1: idx + 1, jdx
            // Neighbor 2: idx    , jdx - 1
            // Neighbor 3: idx    , jdx + 1

            neighbor_idx = ( idx == 0 ? Lx - 1 : idx - 1 );
            neighbor_jdx = jdx;
            neighbor_array[ neighbor_idx * num_neighbors + 0 ] = neighbor_idx * Lx + neighbor_jdx;

            neighbor_idx = ( idx == Lx - 1 ? 0 : idx + 1 );
            neighbor_jdx = jdx;
            neighbor_array[ neighbor_idx * num_neighbors + 1 ] = neighbor_idx * Lx + neighbor_jdx;

            neighbor_idx = idx;
            neighbor_jdx = ( jdx == 0 ? Ly - 1 : jdx - 1 );
            neighbor_array[ neighbor_idx * num_neighbors + 2 ] = neighbor_idx * Lx + neighbor_jdx;

            neighbor_idx = idx;
            neighbor_jdx = ( jdx == Ly - 1 ? 0 : jdx + 1 );
            neighbor_array[ neighbor_idx * num_neighbors + 3 ] = neighbor_idx * Lx + neighbor_jdx;

        }
    }
}

#endif
