#ifndef ASHKIN_TELLER2D_OBSERVABLES
#define ASHKIN_TELLER2D_OBSERVABLES
#include <string>
#include <vector>
#include "ashkin_teller2d_parameters.cxx"
#include <order_parameter_cumulants.hpp>

#if CORRELATION_LENGTHS
// Include the correlation functionality.
#include "../Correlations/fourier_correlator.hpp"
#endif

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
#if CORRELATION_LENGTHS
        sigma_mag, sigma_mag2, sigma_mag4, sigma_corr_qmin,
        tau_mag,   tau_mag2,   tau_mag4, tau_corr_qmin,
        order_param, order_param2, order_param4, order_corr_qmin,
        nem_mag,   nem_mag2,   nem_mag4, nem_corr_qmin,
        counts_per_bin, NUM_OBS
#else
        sigma_mag, sigma_mag2, sigma_mag4, 
        tau_mag,   tau_mag2,   tau_mag4, 
        order_param, order_param2, order_param4,
        nem_mag,   nem_mag2,   nem_mag4,
        counts_per_bin, NUM_OBS
#endif
    };

    enum class nonlinear_obs_enum
    {
#if CORRELATION_LENGTHS
        sigma_susc, sigma_binder, sigma_corr_length,
        tau_susc, tau_binder, tau_corr_length,
        nem_susc, nem_binder_cumulant, nem_corr_length,
        susc, binder_cumulant, corr_length, NUM_OBS
#else
        sigma_susc, sigma_binder, 
        tau_susc, tau_binder,
        nem_susc, nem_binder_cumulant,
        susc, binder_cumulant, NUM_OBS
#endif
    };

#if CORRELATION_LENGTHS
    const std::vector<std::string> string_names = { "Sigma Mag", "Sigma Mag2", "Sigma Mag4", "Sigma G(qmin)" ,
                                                    "Tau Mag",   "Tau Mag2",   "Tau Mag4", "Tau G(qmin)",
                                                    "Order Parameter", "Order Parameter2", "Order Parameter4", "Order Parameter G(qmin)",
                                                    "Nematicity", "Nematicity2", "Nematicity4", "Nematicity G(qmin)",
                                                    "Counts per Bin", "NUM OBS" };

    const std::vector<std::string> nonlinear_obs_strings = { "Sigma Susceptibility", "Sigma Binder Cumulant", "Sigma Correlation Length over L",
                                                             "Tau Susceptibility",   "Tau Binder Cumulant", "Tau Correlation Length over L",
                                                             "Nematicity Susceptibility",   "Nematicity Binder Cumulant", "Nematicity Correlation Length over L",
                                                             "Susceptibility", "Binder Cumulant", "Order Parameter Correlation Length over L"
    };
#else
    const std::vector<std::string> string_names = { "Sigma Mag", "Sigma Mag2", "Sigma Mag4", 
                                                    "Tau Mag",   "Tau Mag2",   "Tau Mag4",
                                                    "Order Parameter", "Order Parameter2", "Order Parameter4",
                                                    "Nematicity", "Nematicity2", "Nematicity4",
                                                    "Counts per Bin", "NUM OBS" };

    const std::vector<std::string> nonlinear_obs_strings = { "Sigma Susceptibility", "Sigma Binder Cumulant",
                                                             "Tau Susceptibility",   "Tau Binder Cumulant",
                                                             "Nematicity Susceptibility",   "Nematicity Binder Cumulant",
                                                             "Susceptibility", "Binder Cumulant" };
#endif
}

constexpr size_t convert(const Obs::enum_names obs_val)
{
    return static_cast<size_t>(obs_val);
}

constexpr size_t convert(const Obs::nonlinear_obs_enum obs_val)
{
    return static_cast<size_t>(obs_val);
}

template<typename data_t>
struct Ashkin_Teller2d_Obs
{
    const size_t num_bins; 

    data_t * obs_array = nullptr;

#if CORRELATION_LENGTHS
    Fourier_Correlator<data_t> correlator;
#endif

#if CORRELATION_LENGTHS
    Ashkin_Teller2d_Obs(const size_t nbins) : num_bins(nbins), correlator( Ashkin_Teller2d_Parameters::L )
#else
    Ashkin_Teller2d_Obs(const size_t nbins) : num_bins(nbins)
#endif
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
#if CORRELATION_LENGTHS
    void update_qmin_correlator(const data_t value, const Obs::enum_names ob, const size_t bin, const size_t counts ) const;
#endif

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

#if CORRELATION_LENGTHS
// This function is necessary so that the correlators will
// only be updated periodically whereas all of the other 
// observables can be computed each time.
template<typename data_t>
void Ashkin_Teller2d_Obs<data_t>::update_qmin_correlator(const data_t value, 
                                                         const Obs::enum_names ob, 
                                                         const size_t bin, const size_t counts ) const
{
    data_t current_avg = get_observable(ob, bin);
    current_avg = ( value + counts * current_avg ) / ( counts + 1 );
    set_observable( current_avg, ob, bin );
}
#endif

// Calculate the thermally-averaged nonlinear observables
// given a thermodynamics object and the thermally-averaged
// linear observables.
template<typename data_t, class thermo_t>
void calculate_nonlinear_observables( const size_t num_temps, const size_t system_size, const thermo_t * const thermo, data_t *& nonlinear_obs )
{
    // Wipe the nonlinear observables completely
    delete [] nonlinear_obs;
    // Create a new array (add an extra index for the temperature)
    const size_t num_nonlinear_obs = static_cast<size_t>(Obs::nonlinear_obs_enum::NUM_OBS);
    nonlinear_obs = new data_t [ num_nonlinear_obs * num_temps ];

#if CORRELATION_LENGTHS
    const size_t Lsize = static_cast<size_t>( sqrt(system_size) );
#endif

    for ( size_t Tidx = 0; Tidx != num_temps; ++Tidx )
    {
        const data_t temperature = static_cast<data_t>( thermo -> temperatures[Tidx] );
        
        // Calculate the sigma susceptibility
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::sigma_susc) ] = calculate_susceptibility<data_t>( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::sigma_mag2) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::sigma_mag) ), temperature, system_size ) / static_cast<data_t>(system_size);

        // Calculate the sigma Binder cumulant
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::sigma_binder) ] = calculate_Binder_cumulant( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::sigma_mag4) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::sigma_mag2) ), system_size );
            
        // Calculate the tau susceptibility
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::tau_susc) ] = calculate_susceptibility<data_t>( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::tau_mag2) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::tau_mag) ), temperature, system_size ) / static_cast<data_t>(system_size);

        // Calculate the tau Binder cumulant
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::tau_binder) ] = calculate_Binder_cumulant( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::tau_mag4) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::tau_mag2) ), system_size );
 
        // Calculate the nematicity susceptibility
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::nem_susc) ] = calculate_susceptibility<data_t>( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::nem_mag2) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::nem_mag) ), temperature, system_size ) / static_cast<data_t>(system_size);

        // Calculate the nematicity Binder cumulant
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::nem_binder_cumulant) ] = calculate_Binder_cumulant( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::nem_mag4) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::nem_mag4) ), system_size );
    
        // Calculate the order parameter susceptibility
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::susc) ] = calculate_susceptibility<data_t>( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::order_param2) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::order_param) ), temperature, system_size ) / static_cast<data_t>(system_size);

        // Calculate the two-component Binder cumulant of the order parameter
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::binder_cumulant) ] = calculate_two_component_Binder_cumulant( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::order_param4) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::order_param2) ), system_size );

#if CORRELATION_LENGTHS
        // Calculate the correlation lengths for all the variables here
        // TODO: This isn't super smart memory-access wise, but it looks better.

        // Sigma Correlation Length
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::sigma_corr_length) ] = calculate_correlation_length( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::sigma_mag2) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::sigma_corr_qmin) ), Lsize ) / Lsize;

        // Tau Correlation Length
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::tau_corr_length) ] = calculate_correlation_length( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::tau_mag2) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::tau_corr_qmin) ), Lsize ) / Lsize;

        // Nematicity Correlation Length
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::nem_corr_length) ] = calculate_correlation_length( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::nem_mag2) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::nem_corr_qmin) ), Lsize ) / Lsize;

        // Order Parameter Correlation Length
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::corr_length) ] = calculate_correlation_length( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::order_mag2) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::order_corr_qmin) ), Lsize ) / Lsize;
#endif

    }
}

#endif
