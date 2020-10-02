#ifndef ISING2D_OBSERVABLES
#include <string>
#include <vector>

static constexpr float DATA_INITIALIZER = 0.;

// Set up an observables enum class.
// These variables are "non-energetic," meaning
// they are not directly obtained by the Wang 
// Landau sampling itself. They instead are 
// average quantities in each bin.
namespace Obs
{
    enum class enum_names
    {
        mag, mag2, mag4, counts_per_bin, NUM_OBS
    };

    const std::vector<std::string> string_names = { "Magnetization", "Magnetization2", "Magnetization4", "Counts per Bin", "NUM OBS" };
}

constexpr size_t convert(const Obs::enum_names obs_val)
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
        obs_array = new data_t [ num_bins * convert(Obs::enum_names::NUM_OBS) ];
        
        // Initialize the observables to zero.
        // Store as AOS across structures since all
        // observables will be accessed in each bin.
        for ( size_t idx = 0; idx != num_bins; ++idx )
        {
            for ( size_t ob = 0; ob != convert(Obs::enum_names::NUM_OBS); ++ob )
                obs_array[ idx * convert(Obs::enum_names::NUM_OBS) + ob ] = static_cast<data_t> (DATA_INITIALIZER);

        }
    }

    ~Ising2d_Obs(){ delete [] obs_array; }

    // Set the data pointed to by the observable array
    void set_observable(const data_t value, const Obs::enum_names ob, const size_t bin) const
    {
        obs_array[ bin * convert(Obs::enum_names::NUM_OBS) + convert(ob) ] = value;
    }
    
    // Get the data pointed to by the observable array
    data_t get_observable(const Obs::enum_names ob, const size_t bin) const
    {
        return obs_array[ bin * convert(Obs::enum_names::NUM_OBS) + convert(ob) ];
    }

    // Get the data pointed to by the observable array
    data_t get_observable(const size_t ob, const size_t bin) const
    {
        return obs_array[ bin * convert(Obs::enum_names::NUM_OBS) + ob ];
    }


    // Update average observable with the given value
    void update_observable_average(const data_t value, const Obs::enum_names ob, const size_t bin) const;

    // Increment the counter
    void increment_counts_per_bin(const size_t bin) const
    {
        ++obs_array[ bin * convert(Obs::enum_names::NUM_OBS) + convert(Obs::enum_names::counts_per_bin) ];
    }

    // Export the observables array as a deep copy.
    void export_observables( data_t *& data_array ) const
    {
        data_array = new data_t [ num_bins * convert(Obs::enum_names::NUM_OBS) ];
        for ( size_t bin = 0; bin != num_bins; ++bin )
        {
            for ( size_t ob = 0; ob != convert(Obs::enum_names::NUM_OBS); ++ob )
                data_array[ bin * convert(Obs::enum_names::NUM_OBS) + ob ] = get_observable(ob, bin);
        }
    }
};

template<typename data_t>
void Ising2d_Obs<data_t>::update_observable_average(const data_t value, 
                                                    const Obs::enum_names ob, 
                                                    const size_t bin) const
{
    data_t current_avg = get_observable(ob, bin);
    data_t counts = get_observable(Obs::enum_names::counts_per_bin, bin);

    current_avg = ( value + counts * current_avg ) / ( counts + 1 );

    set_observable(current_avg, ob, bin);
}

#endif
