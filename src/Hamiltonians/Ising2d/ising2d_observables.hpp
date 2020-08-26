#if ISING2D

static constexpr float DATA_INITIALIZER = 0.;

// Set up an observables enum class.
// These variables are "non-energetic," meaning
// they are not directly obtained by the Wang 
// Landau sampling itself. They instead are 
// average quantities in each bin.
enum class Obs
{
    mag, mag2, mag4, counts_per_bin, NUM_OBS
};

constexpr size_t convert(const Obs obs_val)
{
    return static_cast<size_t>(obs_val);
}

template<typename data_t>
struct Ising2d_Obs
{
    const size_t num_bins; 

    data_t * obs_array = nullptr;

    Ising2d_Obs(const size_t nbins) : num_bins(nbins)
    {
        obs_array = new data_t [ num_samples * convert(NUM_OBS) ];

        // Initialize the observables to zero.
        // Store as AOS across structures since all
        // observables will be accessed in each bin.
        for ( size_t idx = 0; idx != num_bins; ++idx )
            for ( size_t ob = 0; ob != convert(NUM_OBS); ++ob )
                obs_array[ idx * convert(NUM_OBS) + ob ] = static_cast<data_t> (DATA_INITIALIZER);

    }

    ~Ising2d_Obs(){ if ( obs_array != nullptr ) delete [] obs_array; }

    // Set the data pointed to by the observable array
    void set_observable(const data_t value, const Obs ob, const size_t bin) const
    {
        obs_array[ bin * convert(NUM_OBS) + convert(ob) ] = value;
    }
    
    // Get the data pointed to by the observable array
    data_t get_observable(const Obs ob, const size_t bin) const
    {
        return obs_array[ bin * convert(NUM_OBS) + convert(ob) ];
    }

    // Update average observable with the given value
    void update_observable_average(const data_t value, const Obs ob, const size_t bin) const;

    // Increment the counter
    void increment_counts_per_bin(const size_t bin) const
    {
        ++obs_array[ bin * convert(NUM_OBS) + counts_per_bin ];
    }
};

template<typename data_t>
void Ising2d_Obs<data_t>::update_observable_average(const data_t value, 
                                                    const Obs ob, 
                                                    const size_t bin) const
{
    data_t current_avg = get_observable(ob, bin);
    data_t counts = get_observable(counts_per_bin, bin);

    current_avg = ( value + counts * current_avg ) / ( counts + 1 );

    set_observable(current_avg, ob, bin);
}

#endif
