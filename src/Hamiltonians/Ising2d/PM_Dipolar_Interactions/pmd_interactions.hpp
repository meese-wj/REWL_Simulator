#ifndef PMD_INTERACTIONS
#define PMD_INTERACTIONS

#include <cstdint>
#include <iostream>
#include <string>
#include <mpi.h>
#include "../../grid/pbc_2d_vectors.hpp"

// Include helper functions for 
// the interaction itself.

// Return the value of the interaction without
// the coupling.
//
// interaction = cos(4 theta) / | site_1 - site_2 |^2
//
// where cos(theta) = (site_1 - site_2).x / | site_1 - site_2 |.
// Note that cos(4x) = 1 - 8 * cos(x)^2 + 8 * cos(x)^4.
template<typename energy_t>
energy_t interaction_value( const pbc_2d_vector<energy_t> & site_1, const pbc_2d_vector<energy_t> site_2, 
                            const std::uint32_t Lx, const std::uint32_t Ly )
{
    pbc_2d_vector<energy_t> difference = pbc_subtract<energy_t>(Lx, Ly, site_1, site_2);
    if (site_1.x == site_2.x && site_1.y == site_2.y)
        return 0.;

    int my_id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    energy_t xprojection = difference.x;
    energy_t cos_theta = xprojection / magnitude<energy_t, energy_t>( difference );
    energy_t cos_theta_sq = cos_theta * cos_theta;

    if ( site_1.x == 0. && site_1.y == 0. && my_id == 0 )
    {
        std::cout << "\n";
        std::cout << "\nsite_1      = (" << site_1.x << ", " << site_1.y << ")";
        std::cout << "\nsite_2      = (" << site_2.x << ", " << site_2.y << ")";
        std::cout << "\ndiff        = (" << difference.x << ", " << difference.y << ")";
        std::cout << "\ncos(theta)  = " << cos_theta;
        std::cout << "\ncos(4theta)  = " << 1. - 8. * cos_theta_sq + 8. * cos_theta_sq * cos_theta_sq;
        std::cout << "\n";
    }

    return (1. - 8. * cos_theta_sq + 8. * cos_theta_sq * cos_theta_sq) / square_magnitude<energy_t, energy_t>( difference );
}

// Encapsulate the massive phonon-mediated
// interactions inside of a struct/functor 
// thing because it will probably require 
// more work.
template<typename energy_t, typename spin_t>
struct PMDN_Interactions
{
    const std::uint32_t total_sites;
    energy_t * interaction_matrix = nullptr;
    PMDN_Interactions( const energy_t phonon_coupling, const std::uint32_t Lx, const std::uint32_t Ly ) : total_sites(Lx * Ly)
    {
        build_interaction_matrix( phonon_coupling, Lx, Ly);
    }

    void build_interaction_matrix( const energy_t phonon_coupling, const std::uint32_t Lx, const std::uint32_t Ly );

    energy_t interaction_strength( const std::uint32_t site, const std::uint32_t site_2 ) const
    {
        return interaction_matrix[ site * total_sites + site_2 ];
    }
    energy_t calculate_energy_per_spin( const std::uint32_t site, spin_t spin_val, const spin_t * const all_spins ) const;
    energy_t calculate_total_PMDN_energy( const spin_t * const all_spins ) const;

    ~PMDN_Interactions( )
    {
        delete [] interaction_matrix;
    }
};

template<typename energy_t, typename spin_t>
void PMDN_Interactions<energy_t, spin_t>::build_interaction_matrix(const energy_t phonon_coupling, const std::uint32_t Lx, const std::uint32_t Ly)
{
    std::uint32_t system_size = Lx * Ly;
    std::string printer = "\n";
    interaction_matrix = new energy_t [ system_size * system_size ]();
    int my_id = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, &my_id );
    for ( std::uint32_t x1 = 0; x1 != Lx; ++x1 )
        for ( std::uint32_t y1 = 0; y1 != Ly; ++y1 )
        {
            std::uint32_t index_1 = x1 * Lx + y1;
            pbc_2d_vector<energy_t> site_1( x1, y1 );
            
            printer += std::to_string( index_1 ) + ":  ";
            for ( std::uint32_t x2 = 0; x2 != Lx; ++x2 )
                for ( std::uint32_t y2 = 0; y2 != Ly; ++y2 )
                {    
                    std::uint32_t index_2 = x2 * Lx + y2;
                    pbc_2d_vector<energy_t> site_2( x2, y2 );
                    if (my_id == 0 && (x1 == 0 && y1 == 0))
                        std::cout << "\n(" << x2 << ", " << y2 << ") -> " << index_2 << "\n"; 
                    interaction_matrix[ index_1 * total_sites + index_2 ] = 0.;
                    if ( !( x1 == x2 && y1 == y2 ) )
                        interaction_matrix[ index_1 * total_sites + index_2 ] = phonon_coupling * interaction_value( site_1, site_2, Lx, Ly );
                    printer += std::to_string( interaction_matrix[ index_1 * total_sites + index_2 ] / phonon_coupling ) + " ";
                }
            printer += "\n";
        }

    printer += "\n";
    if ( my_id == 0 )
        std::cout << printer;
    return;
}

template<typename energy_t, typename spin_t>
energy_t PMDN_Interactions<energy_t, spin_t>::calculate_energy_per_spin( const std::uint32_t site, const spin_t spin_val, const spin_t * const all_spins ) const
{
    energy_t summand = 0.;
    // First sum along the interaction matrix
    for (std::uint32_t idx = 0; idx != total_sites; ++idx)
    {
        // The interaction strength of (site, site) is defined as zero
        // to not accidentally introduce code branches into this behemoth...
        summand += interaction_strength( site, idx ) * static_cast<energy_t>( all_spins[idx] );
    }
    return static_cast<energy_t>(spin_val) * summand;    
}

template<typename energy_t, typename spin_t>
energy_t PMDN_Interactions<energy_t, spin_t>::calculate_total_PMDN_energy( const spin_t * const all_spins ) const
{
    energy_t total_energy = 0.;
    for ( std::uint32_t site = 0; site != total_sites; ++site )
    {
        total_energy += calculate_energy_per_spin( site, all_spins[site], all_spins );
    }
    return 0.5 * total_energy;
}

// Use a proxy interaction scheme to calculate 
// the ground state contribution from these 
// interactions.
template<typename energy_t, typename spin_t>
energy_t PMNI_ground_state_contribution( const energy_t phonon_coupling, const std::uint32_t Lx, const std::uint32_t Ly, const energy_t moment_value=1. )
{
    // Assuming everything is ferromagnetic
    PMDN_Interactions<energy_t, spin_t> pmd_ints( phonon_coupling, Lx, Ly );
    spin_t * ferromagnetic_ground_state_spins = new spin_t [ Lx * Ly ];
    for (std::uint32_t site = 0; site != Lx * Ly; ++site )
        ferromagnetic_ground_state_spins[site] = moment_value;

    energy_t total_energy = pmd_ints.calculate_total_PMDN_energy( ferromagnetic_ground_state_spins );
    std::cout << "Total energy = " << total_energy << "\n";
    delete [] ferromagnetic_ground_state_spins;
    return total_energy;
}

#endif // PMD_INTERACTIONS
