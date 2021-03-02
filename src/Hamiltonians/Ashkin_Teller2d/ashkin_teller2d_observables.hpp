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

#if AT_DENSITIES
// Include the 2d histograms for 
// sigma and tau
#include "Density_Plots/ashkin_teller_densities_parameters.cxx"
#include "Density_Plots/ashkin_teller_densities.hpp"
#endif

static constexpr float DATA_INITIALIZER = 0.;

// Set up an observables enum class.
// These variables are "non-energetic," meaning
// they are not directly obtained by the Wang 
// Landau sampling itself. They instead are 
// average quantities in each bin.
namespace Obs
{
    // The "Baxter variable" here is product sigma * tau.
    // The term comes from the J1-J2 antiferromagnetic
    // Heisenberg model.
    // Phi represents the radial variable in the sigma-
    // tau plane.
    enum class enum_names
    {
#if CORRELATION_LENGTHS
        sigma_mag, sigma_mag2, sigma_mag4, sigma_corr_qmin,
        tau_mag,   tau_mag2,   tau_mag4, tau_corr_qmin,
        phi, phi2, phi4, phi_corr_qmin,
        baxter_mag,   baxter_mag2,   baxter_mag4, baxter_corr_qmin,
        counts_per_bin, NUM_OBS
#else
        sigma_mag, sigma_mag2, sigma_mag4, 
        tau_mag,   tau_mag2,   tau_mag4, 
        phi, phi2, phi4,
        baxter_mag,   baxter_mag2,   baxter_mag4,
        counts_per_bin, NUM_OBS
#endif
    };

    enum class nonlinear_obs_enum
    {
#if CORRELATION_LENGTHS
        sigma_susc, sigma_binder, sigma_corr_length,
        tau_susc, tau_binder, tau_corr_length,
        baxter_susc, baxter_binder_cumulant, baxter_corr_length,
        phi_susc, phi_binder_cumulant, phi_corr_length, NUM_OBS
#else
        sigma_susc, sigma_binder, 
        tau_susc, tau_binder,
        baxter_susc, baxter_binder_cumulant,
        phi_susc, phi_binder_cumulant, NUM_OBS
#endif
    };

#if CORRELATION_LENGTHS
    const std::vector<std::string> string_names = { "Sigma Mag", "Sigma Mag2", "Sigma Mag4", "Sigma G(qmin)" ,
                                                    "Tau Mag",   "Tau Mag2",   "Tau Mag4", "Tau G(qmin)",
                                                    "Phi", "Phi2", "Phi4", "Phi G(qmin)",
                                                    "Baxter", "Baxter2", "Baxter4", "Baxter G(qmin)",
                                                    "Counts per Bin", "NUM OBS" };

    const std::vector<std::string> nonlinear_obs_strings = { "Sigma Susceptibility", "Sigma Binder Cumulant", "Sigma Correlation Length over L",
                                                             "Tau Susceptibility",   "Tau Binder Cumulant", "Tau Correlation Length over L",
                                                             "Baxter Susceptibility",   "Baxter Binder Cumulant", "Baxter Correlation Length over L",
                                                             "Phi Susceptibility", "Phi Binder Cumulant", "Phi Correlation Length over L"
    };
#else
    const std::vector<std::string> string_names = { "Sigma Mag", "Sigma Mag2", "Sigma Mag4", 
                                                    "Tau Mag",   "Tau Mag2",   "Tau Mag4",
                                                    "Phi", "Phi2", "Phi4",
                                                    "Baxter", "Baxter2", "Baxter4",
                                                    "Counts per Bin", "NUM OBS" };

    const std::vector<std::string> nonlinear_obs_strings = { "Sigma Susceptibility", "Sigma Binder Cumulant",
                                                             "Tau Susceptibility",   "Tau Binder Cumulant",
                                                             "Baxter Susceptibility",   "Baxter Binder Cumulant",
                                                             "Phi Susceptibility", "Phi Binder Cumulant" };
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

#if AT_DENSITIES
    density_int * density_histograms = nullptr;
    density_float * density_float_data = nullptr;
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

#if AT_DENSITIES
        density_histograms = new density_int [ nbins * AT_Density_Parameters::total_bins ]();
        density_float_data = new density_float [ nbins * AT_Density_Parameters::total_bins ]();
#endif
    }

    ~Ashkin_Teller2d_Obs()
    { 
        delete [] obs_array; 
#if AT_DENSITIES
        delete [] density_histograms;
        delete [] density_float_data;
#endif
    }

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

#if AT_DENSITIES
    // Export the density floats from the observables
    void export_density_plots( std::vector<std::vector<density_float> > & export_vectors )
    {
        // Resize the export vectors
        export_vectors.resize( num_bins );
        for ( auto &v : export_vectors )
            v.resize( AT_Density_Parameters::total_bins );
        
        // Now export the density plots
        for ( size_t idx = 0; idx != num_bins; ++idx )
        {
            export_vectors[ idx ] = std::vector<density_float> ( density_float_data + energy_bin_density_pointer(idx), density_float_data + energy_bin_density_pointer( idx + 1 ) );
        }
    }
#endif
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
 
        // Calculate the Baxter susceptibility
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::baxter_susc) ] = calculate_susceptibility<data_t>( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::baxter_mag2) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::baxter_mag) ), temperature, system_size ) / static_cast<data_t>(system_size);

        // Calculate the Baxter Binder cumulant
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::baxter_binder_cumulant) ] = calculate_Binder_cumulant( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::baxter_mag4) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::baxter_mag2) ), system_size );
    
        // Calculate the phi susceptibility
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::phi_susc) ] = calculate_susceptibility<data_t>( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::phi2) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::phi) ), temperature, system_size ) / static_cast<data_t>(system_size);

        // Calculate the two-component Binder cumulant of phi
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::phi_binder_cumulant) ] = calculate_two_component_Binder_cumulant( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::phi4) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::phi2) ), system_size );

#if CORRELATION_LENGTHS
        // Calculate the correlation lengths for all the variables here
        // TODO: This isn't super smart memory-access wise, but it looks better.

        // Sigma Correlation Length
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::sigma_corr_length) ] = calculate_correlation_length( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::sigma_mag2) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::sigma_corr_qmin) ), Lsize ) / Lsize;

        // Tau Correlation Length
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::tau_corr_length) ] = calculate_correlation_length( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::tau_mag2) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::tau_corr_qmin) ), Lsize ) / Lsize;

        // Baxter Correlation Length
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::baxter_corr_length) ] = calculate_correlation_length( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::baxter_mag2) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::baxter_corr_qmin) ), Lsize ) / Lsize;

        // phi Correlation Length
        nonlinear_obs[ Tidx * num_nonlinear_obs + convert(Obs::nonlinear_obs_enum::phi_corr_length) ] = calculate_correlation_length( thermo -> get_system_obs( Tidx, convert(Obs::enum_names::phi2) ), thermo -> get_system_obs( Tidx, convert(Obs::enum_names::phi_corr_qmin) ), Lsize ) / Lsize;
#endif

    }
}

#endif
