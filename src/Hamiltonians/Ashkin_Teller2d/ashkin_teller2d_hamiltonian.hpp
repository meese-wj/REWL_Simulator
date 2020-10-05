#ifndef ASKIN_TELLER2D_HAMILTONIAN
#define ASKIN_TELLER2D_HAMILTONIAN

// Include the parameters 
#include "ashkin_teller2d_parameters.cxx"
#include "ashkin_teller2d_parameter_string.hpp"

// Include the observables enum class
#include "ashkin_teller2d_observables.hpp"

// Include the grid setup from the 
// cmake include directories.
#include <grid_setup.hpp>

enum spin_type
{
    sigma, tau, NUM_SPIN_TYPES
};

// TODO: Upgrade the energies to allow for 
// double calculations. 
template<typename data_t>
struct State
{
    short which_to_update = spin_type::sigma;      
    float energy = 0;
    data_t sigma_magnetization = 0;
    data_t tau_magnetization = 0;
    data_t nematicity = 0.;
    data_t DoF [spin_type::NUM_SPIN_TYPES] = { 1., 1. }; // Store the local degree of freedom
};

template<typename data_t>
void print(const State<data_t> & stat)
{
    printf("\nState:");
    printf("\n\twhich to update     = %s", stat.which_to_update == spin_type::sigma ? "sigma" : "tau");
    printf("\n\tenergy              = %e", stat.energy);
    printf("\n\tsigma magnetization = %e", stat.sigma_magnetization);
    printf("\n\ttau magnetization   = %e", stat.tau_magnetization);
    printf("\n\tnematicity          = %e", stat.nematicity);
    printf("\n\tDoF                 = %e\n", stat.DoF);
}

template<typename data_t>
struct Ashkin_Teller2d
{
    State<data_t> current_state;

    data_t * spin_array = nullptr;
    size_t * neighbor_array = nullptr;

    // Add some Hamiltonian dependent functions
    float sigma_local_field(const size_t idx) const;
    float tau_local_field(const size_t idx) const;
    float sigma_local_energy(const size_t idx, const data_t spin_value) const;
    float tau_local_energy(const size_t idx, const data_t spin_value) const;
    void recalculate_state();
    void change_state(const size_t idx, State<data_t> & temp_state);
    void set_state(const size_t idx, const State<data_t> & _state);
    void update_observables(const size_t bin, Ashkin_Teller2d_Obs<data_t> * obs_ptr) const;

    // Finally add the constructor and destructor.
    Ashkin_Teller2d()
    {
        spin_array = new data_t [ static_cast<size_t>(spin_type::NUM_SPIN_TYPES) * Ashkin_Teller2d_Parameters::N ];
        // TODO: change this to be randomized?
        for ( size_t idx = 0; idx != static_cast<size_t>(spin_type::NUM_SPIN_TYPES) * Ashkin_Teller2d_Parameters::N; ++idx )
            spin_array[idx] = 1.;

        // TODO: Generalize this to different types 
        // of grids.
        // Allocate the neighbor array
        define_2d_square_periodic_neighbors(Ashkin_Teller2d_Parameters::L,
                                            Ashkin_Teller2d_Parameters::L,
                                            Ashkin_Teller2d_Parameters::num_neighbors_i,
                                            neighbor_array);

        recalculate_state();
    }

    // Get sigma or tau at a site
    data_t * spin_at_site(const size_t site, const spin_type type) const {   return &( spin_array[ static_cast<size_t>(spin_type::NUM_SPIN_TYPES) * site + static_cast<size_t>(type)] ); }

    ~Ashkin_Teller2d()
    {
        delete [] spin_array;
        delete [] neighbor_array;
    }

    void print_lattice() const;
};

// The local field is defined such that the local energy at site
// i is 
// sigma_i * sum_{j neighbors} ( -J * sigma_j - K * sigma_j * tau_i * tau_j )
template<typename data_t>
float Ashkin_Teller2d<data_t>::sigma_local_field(const size_t idx) const
{
    float field = 0.;
    for ( size_t nidx = 0; nidx != Ashkin_Teller2d_Parameters::num_neighbors_i; ++nidx )
    {
        const size_t neighbor = neighbor_array[ idx * Ashkin_Teller2d_Parameters::num_neighbors_i + nidx ];

        field += Ashkin_Teller2d_Parameters::J * static_cast<float>( *spin_at_site( neighbor, spin_type::sigma ) );

        field += Ashkin_Teller2d_Parameters::K * static_cast<float>(   (*spin_at_site(idx,      spin_type::tau))
                                                                     * (*spin_at_site(neighbor, spin_type::sigma)) 
                                                                     * (*spin_at_site(neighbor, spin_type::tau))   );
    }
    return field;
}

// The local field is defined such that the local energy at site
// i is 
// tau_i * sum_{j neighbors} ( -J * tau_j - K * tau_j * sigma_i * sigma_j )
template<typename data_t>
float Ashkin_Teller2d<data_t>::tau_local_field(const size_t idx) const
{
    float field = 0.;
    for ( size_t nidx = 0; nidx != Ashkin_Teller2d_Parameters::num_neighbors_i; ++nidx )
    {
        const size_t neighbor = neighbor_array[ idx * Ashkin_Teller2d_Parameters::num_neighbors_i + nidx ];

        field += Ashkin_Teller2d_Parameters::J * static_cast<float>( *spin_at_site( neighbor, spin_type::tau ) );

        field += Ashkin_Teller2d_Parameters::K * static_cast<float>(   (*spin_at_site(idx,      spin_type::sigma))
                                                                     * (*spin_at_site(neighbor, spin_type::tau)) 
                                                                     * (*spin_at_site(neighbor, spin_type::sigma)) );
    }
    return field;
}

template<typename data_t>
float Ashkin_Teller2d<data_t>::sigma_local_energy(const size_t idx, const data_t spin_value) const
{
    return -1. * static_cast<float>( spin_value ) * sigma_local_field(idx);
}

template<typename data_t>
float Ashkin_Teller2d<data_t>::tau_local_energy(const size_t idx, const data_t spin_value) const
{
    return -1. * static_cast<float>( spin_value ) * tau_local_field(idx);
}


template<typename data_t>
void Ashkin_Teller2d<data_t>::recalculate_state()
{
    float temp_energy      = 0.;
    data_t temp_sigma_mag  = 0.;
    data_t temp_tau_mag    = 0.;
    data_t temp_nematicity = 0.;

    for ( size_t idx = 0; idx != Ashkin_Teller2d_Parameters::N; ++idx )
    {
       temp_sigma_mag += *spin_at_site(idx, spin_type::sigma);
       temp_tau_mag   += *spin_at_site(idx, spin_type::tau);
       temp_nematicity += (*spin_at_site(idx, spin_type::sigma)) * (*spin_at_site(idx, spin_type::tau));
       temp_energy += 0.5 * ( sigma_local_energy(idx, *spin_at_site(idx, spin_type::sigma)) + tau_local_energy(idx, *spin_at_site(idx, spin_type::tau)) );
       // These next terms are necessary to avoid the double counting
       // that occurs in the K term.
       for ( size_t nidx = 0; nidx != Ashkin_Teller2d_Parameters::num_neighbors_i; ++nidx )
       {
           const size_t neighbor = neighbor_array[ idx * Ashkin_Teller2d_Parameters::num_neighbors_i + nidx ];
           temp_energy += 0.5 * Ashkin_Teller2d_Parameters::K * static_cast<float>(   (*spin_at_site(idx, spin_type::sigma)) * (*spin_at_site(idx, spin_type::tau))
                                                                                    * (*spin_at_site(neighbor, spin_type::sigma))
                                                                                    * (*spin_at_site(neighbor, spin_type::tau))                                     );
       }
    }

    current_state.energy              = temp_energy;
    current_state.sigma_magnetization = temp_sigma_mag;
    current_state.tau_magnetization   = temp_tau_mag;
    current_state.nematicity          = temp_nematicity;
}

// Change the state by changing the value of the spin
// at a particular site idx. Return what the new state
// would be if the change is accepted (the state is
// passed by reference to avoid extra copies).
template<typename data_t>
void Ashkin_Teller2d<data_t>::change_state(const size_t idx, State<data_t> & temp_state)
{
    switch(current_state.which_to_update)
    {
        case spin_type::sigma:
        {
            temp_state.which_to_update       = spin_type::sigma;     // Tell the temp_state that this is a sigma update
            temp_state.DoF[spin_type::sigma] = -1. * (*spin_at_site(idx, spin_type::sigma));    
            temp_state.DoF[spin_type::tau]   = *spin_at_site(idx, spin_type::tau);    

            temp_state.sigma_magnetization = current_state.sigma_magnetization + temp_state.DoF[spin_type::sigma] - (*spin_at_site(idx, spin_type::sigma));
            temp_state.tau_magnetization   = current_state.tau_magnetization;
            temp_state.nematicity          = current_state.nematicity + temp_state.DoF[spin_type::tau] * ( temp_state.DoF[spin_type::sigma] - (*spin_at_site(idx, spin_type::sigma)) );
            temp_state.energy              = current_state.energy + sigma_local_energy( idx, temp_state.DoF[spin_type::sigma] ) - sigma_local_energy( idx, *spin_at_site(idx, spin_type::sigma) );

            break;
        }
        case spin_type::tau:
        {
            temp_state.which_to_update       = spin_type::tau;       // Tell the temp_state that this is a tau update
            temp_state.DoF[spin_type::tau]   = -1. * (*spin_at_site(idx, spin_type::tau));    
            temp_state.DoF[spin_type::sigma] = *spin_at_site(idx, spin_type::sigma);
            
            temp_state.sigma_magnetization = current_state.sigma_magnetization;
            temp_state.tau_magnetization   = current_state.tau_magnetization + temp_state.DoF[spin_type::tau] - (*spin_at_site(idx, spin_type::tau));
            temp_state.nematicity          = current_state.nematicity + temp_state.DoF[spin_type::sigma] * ( temp_state.DoF[spin_type::tau] - (*spin_at_site(idx, spin_type::tau)) );
            temp_state.energy              = current_state.energy + tau_local_energy( idx, temp_state.DoF[spin_type::tau] ) - tau_local_energy( idx, *spin_at_site(idx, spin_type::tau) );

            break;
        }
    }

    // After a single sweep, switch the spin type being updated
    if ( idx == Ashkin_Teller2d_Parameters::N - 1 )
    {
        switch(current_state.which_to_update)
        {
            case spin_type::sigma: 
            {
                current_state.which_to_update = spin_type::tau;
                break;
            }
            case spin_type::tau: 
            {
                current_state.which_to_update = spin_type::sigma;
                break;
            }
        }
    }
}

// Set the current state
template<typename data_t>
void Ashkin_Teller2d<data_t>::set_state(const size_t idx, const State<data_t> & _state)
{
    current_state.energy                 = _state.energy;

    current_state.sigma_magnetization    = _state.sigma_magnetization;
    current_state.tau_magnetization      = _state.tau_magnetization;
    current_state.nematicity             = _state.nematicity;

    *spin_at_site(idx, spin_type::sigma) = _state.DoF[spin_type::sigma];
    *spin_at_site(idx, spin_type::tau)   = _state.DoF[spin_type::tau];
}

// Update the non-energetic observables
template<typename data_t>
void Ashkin_Teller2d<data_t>::update_observables(const size_t bin, Ashkin_Teller2d_Obs<data_t> * obs_ptr) const
{
    const data_t sigma_val = abs( current_state.sigma_magnetization );
    const data_t tau_val   = abs( current_state.tau_magnetization );
    const data_t nem_val   = abs( current_state.nematicity );
    
    obs_ptr -> update_observable_average(sigma_val, Obs::enum_names::sigma_mag, bin);
    obs_ptr -> update_observable_average(sigma_val * sigma_val, Obs::enum_names::sigma_mag2, bin);
    obs_ptr -> update_observable_average(sigma_val * sigma_val * sigma_val * sigma_val, Obs::enum_names::sigma_mag4, bin);

    obs_ptr -> update_observable_average(tau_val, Obs::enum_names::tau_mag, bin);
    obs_ptr -> update_observable_average(tau_val * tau_val, Obs::enum_names::tau_mag2, bin);
    obs_ptr -> update_observable_average(tau_val * tau_val * tau_val * tau_val, Obs::enum_names::tau_mag4, bin);

    obs_ptr -> update_observable_average(nem_val, Obs::enum_names::nem_mag, bin);
    obs_ptr -> update_observable_average(nem_val * nem_val, Obs::enum_names::nem_mag2, bin);
    obs_ptr -> update_observable_average(nem_val * nem_val * nem_val * nem_val, Obs::enum_names::nem_mag4, bin);
    
    obs_ptr -> increment_counts_per_bin(bin);
}

template<typename data_t>
void Ashkin_Teller2d<data_t>::print_lattice() const
{
    printf("\n\nSigmas +/-");
    for ( size_t idx = 0; idx != Ashkin_Teller2d_Parameters::L; ++idx )
    {
        printf("\n%ld\t    ", idx * Ashkin_Teller2d_Parameters::L);
        for ( size_t jdx = 0; jdx != Ashkin_Teller2d_Parameters::L; ++jdx )
        {
            printf("%c ", *spin_at_site( idx * Ashkin_Teller2d_Parameters::L + jdx, spin_type::sigma ) == 1. ? '+' : '-');   
        }
    }
    printf("\nTaus */o");
    for ( size_t idx = 0; idx != Ashkin_Teller2d_Parameters::L; ++idx )
    {
        printf("\n%ld\t    ", idx * Ashkin_Teller2d_Parameters::L);
        for ( size_t jdx = 0; jdx != Ashkin_Teller2d_Parameters::L; ++jdx )
        {
            printf("%c ", *spin_at_site( idx * Ashkin_Teller2d_Parameters::L + jdx, spin_type::tau ) == 1. ? '*' : 'o');   
        }
    }
}

#endif


