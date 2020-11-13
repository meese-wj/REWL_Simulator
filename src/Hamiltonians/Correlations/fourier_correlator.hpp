#ifndef FOURIER_CORRELATOR
#define FOURIER_CORRELATOR
/* Use this file to set up the correlation
 * length measurements. To evaluate the 
 * correlation length, we employ the proxy 
 * given in 
 * https://doi.org/10.1103/PhysRevB.91.224201 */

#include <cmath>

constexpr size_t counts_per_transform = 1000;

// We assume that Lx = Ly for now.
template<typename obs_t>
struct Transform_Kernel
{
    const size_t L;
    const size_t N = L * L;
    const obs_t qmin = 2. * acos(-1.) / static_cast<obs_t>(L);
    
    obs_t * real_part = nullptr;
    obs_t * imag_part = nullptr;

    Transform_Kernel( const size_t _L ) : L(_L)
    {
        // For each qmin, there should be
        // a real part and imaginary part the
        // size of the lattice. But since 
        // qmin = (2pi/L, 0) or (0, 2pi/L),
        // the function only has L unique
        // values.
        real_part = new obs_t [L];
        imag_part = new obs_t [L];

        for ( size_t idx = 0; idx != L; ++idx )
        {
            real_part[ idx ] = cos( qmin * static_cast<obs_t>(idx) );
            imag_part[ idx ] = sin( qmin * static_cast<obs_t>(idx) );
            //real_part[ idx ] = cos( qmin * static_cast<obs_t>(idx) ) / sqrt( static_cast<obs_t>(N) );
            //imag_part[ idx ] = sin( qmin * static_cast<obs_t>(idx) ) / sqrt( static_cast<obs_t>(N) );
        }
    }
    
    ~Transform_Kernel()
    {
        delete [] real_part;
        delete [] imag_part;
    }
};

template<typename obs_t>
struct Fourier_Correlator
{
    const Transform_Kernel<obs_t> kernel;

    Fourier_Correlator( const size_t _L ) : kernel(_L)
    {}

    ~Fourier_Correlator()
    {}

    obs_t compute_correlator( const obs_t * const field_array ) const;
};

template<typename obs_t>
obs_t Fourier_Correlator<obs_t>::compute_correlator( const obs_t * const field_array ) const
{
    const size_t L = kernel.L;
    obs_t corr_real_x = 0.;
    obs_t corr_imag_x = 0.;
    obs_t corr_real_y = 0.;
    obs_t corr_imag_y = 0.;
    obs_t temp_corr_x = 0.;
    obs_t temp_kernel_real = 0.;
    obs_t temp_kernel_imag = 0.;
    
    // Calculate the Correlator by summing over
    // idx and jdx (x and y axes). The correlators
    // only vary along the x and y axes and so 
    // the correlators are evaluated in a way 
    // meant to speed things up. Hopefully it works
    // since it's so expensive.
    for ( size_t idx = 0; idx != L; ++idx )
    {
        temp_corr_x = 0.;
        for ( size_t jdx = 0; jdx != L; ++jdx )
        {
            temp_corr_x += field_array[ idx * L + jdx ];
            corr_real_y += field_array[ idx * L + jdx ] * kernel.real_part[ jdx ];
            corr_imag_y += field_array[ idx * L + jdx ] * kernel.imag_part[ jdx ];

            if ( idx == jdx )
            {
                temp_kernel_real = kernel.real_part[ jdx ];
                temp_kernel_imag = kernel.imag_part[ jdx ];
            }
        }
        corr_real_x += temp_corr_x * temp_kernel_real;
        corr_imag_x += temp_corr_x * temp_kernel_imag;
    }
   
    // Store the complex field transform square-magnitude in the real part.
    // This is the equivalent of the Fourier Correlator for 
    // qmin = (2pi/L, 0) and qmin = (0, 2pi/L).
    corr_real_x = corr_real_x * corr_real_x + corr_imag_x * corr_imag_x;
    corr_real_y = corr_real_y * corr_real_y + corr_imag_y * corr_imag_y;

    // Finally average the two directions for lattices with
    // tetragonal symmetry.
    return 0.5 * ( corr_real_x + corr_real_y );
}

#endif
