#ifndef ASHKIN_TELLER2D_OBSERVABLES
#define ASHKIN_TELLER2D_OBSERVABLES
#include <string>
#include <vector>
#include <order_parameter_cumulants.hpp>

static constexpr float DATA_INITIALIZER = 0.;

// Set up an observables enum class.
// These variables are "non-energetic," meaning
// they are not directly obtained by the Wang 
// Landau sampling itself. They instead are 
// average quantities in each bin.
namespace Obs
{
    // The "nematicity" here is product sigma * tau.
    // The term comes from the J1-J2 antiferromagnetic
    // Heisenberg model.
    enum class enum_names
    {
        sigma_mag, sigma_mag2, sigma_mag4, 
        tau_mag,   tau_mag2,   tau_mag4, 
        order_param, order_param2, order_param4,
        nem_mag,   nem_mag2,   nem_mag4,
        counts_per_bin, NUM_OBS
    };

    enum class nonlinear_obs_enum
    {
        order_param, susc, binder_cumulant, NUM_OBS
    };

    const std::vector<std::string> string_names = { "Sigma Mag", "Sigma Mag2", "Sigma Mag4", 
                                                    "Tau Mag",   "Tau Mag2",   "Tau Mag4",
                                                    "Order Parameter", "Order Parameter2", "Order Parameter4",
                                                    "Nematicity", "Nematicity2", "Nematicity4",
                                                    "Counts per Bin", "NUM OBS" };

    const std::vector<std::string> nonlinear_obs_strings = { "Order Parameter", "Susceptibility", "Binder Cumulant" };
}

constexpr size_t convert(const Obs::enum_names obs_val)
{
    return static_cast<size_t>(obs_val);
}

template<typename data_t>
struct Ashkin_Teller2d_Obs
{
    const size_t num_bins; 

    data_t * obs_array = nullptr;

    Ashkin_Teller2d_Obs(const size_t nbins) : num_bins(nbins)
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

    ~Ashkin_Teller2d_Obs(){ delete [] obs_array; }

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
void Ashkin_Teller2d_Obs<data_t>::update_observable_average(const data_t value, 
                                                            const Obs::enum_names ob, 
                                                            const size_t bin) const
{
    data_t current_avg = get_observable(ob, bin);
    data_t counts = get_observable(Obs::enum_names::counts_per_bin, bin);

    current_avg = ( value + counts * current_avg ) / ( counts + 1 );

    set_observable(current_avg, ob, bin);
}


// Calculate the thermally-averaged nonlinear observables
// given a thermodynamics enum and the thermally-averaged
// linear observables.
template<typename data_t, typename thermo_enum, class thermo_t>
void calculate_nonlinear_observables( const size_t num_temps, const size_t system_size, const * thermo_t const thermo, data_t *& nonlinear_obs )
{
    // Wipe the nonlinear observables completely
    delete [] nonlinear_obs;
    // Create a new array (add an extra index for the temperature)
    const size_t num_nonlinear_obs = static_cast<size_t>(Obs::nonlinear_obs_enum::NUM_OBS);
    nonlinear_obs = new data_t [ num_nonlinear_obs * num_temps ];

    for ( size_t Tidx = 0; Tidx != num_temps; ++Tidx )
    {
        const data_t temperature = static_cast<data_t>( thermo -> temperatures[Tidx] );

        // Calculate the susceptibility
        nonlinear_obs[ Tidx * num_nonlinear_obs + Obs::nonlinear_obs_enum::susc ] = calculate_susceptibility<data_t>( thermo -> get_system_obs( Tidx, Obs::enum_names::order_param2 ), thermo -> get_system_obs( Tidx, Obs::enum_names::order_param ), temperature, system_size );

        // Calculate the two-component Binder cumulant
        nonlinear_obs[ Tidx * num_nonlinear_obs + Obs::nonlinear_obs_enum::binder_cumulant ] = calculate_two_component_Binder_cumulant( thermo -> get_system_obs( Tidx, Obs::enum_names::order_param4 ), thermo -> get_system_obs( Tidx, Obs::enum_names::order_param2 ), system_size );
    }
}

#endif
