#ifndef THERMODYNAMICS
/* Create an object to calculate the 
 * canonical thermodynamics from a density
 * of states using trapezoidal rule. */

#include <cmath>

constexpr size_t NUM_ENERGY_OBS = 5; 

template<class Obs_enum_t>
constexpr size_t convert( const Obs_enum_t ob )
{
    return static_cast<size_t>(ob);
}

template<typename energy_t, typename logdos_t, typename obs_t,
         class Obs_enum_t, class Obs_container_t>
struct Thermodynamics
{
   const energy_t energy_min;
   const energy_t energy_max;
   const energy_t energy_bin_size;

   const energy_t Tmin;
   const energy_t Tmax; 
   const size_t num_T;
   const energy_t dT = (Tmax - Tmin) / static_cast<energy_t>(num_T);

   obs_t * canonical_observables = nullptr;

   Thermodynamics(const energy_t _emin, const energy_t _emax, const energy_t _ebsize,
                  const energy_t _Tmin, const energy_t _Tmax, const size_t _nT) : energy_min(_emin), energy_max(_emax), energy_bin_size(_ebsize),
                                                                                  Tmin(_Tmin), Tmax(_Tmax), num_T(_nT)
    {
        canonical_observables = new obs_t [ num_T * ( NUM_ENERGY_OBS + convert<Obs_enum_t>(Obs_enum_t::NUM_OBS) ) ];
    }

    ~Thermodynamics(){ if (canonical_observables != nullptr) delete [] canonical_observables; }

    obs_t energy_observable_value( const Energy_Obs e_ob, const energy_t Tvalue ) const; 
    
    obs_t observable_value( const Obs_enum_t ob, const energy_t Tvalue ) const; 

    void calculate_thermodynamics( const logdos_t * const logdos_array,
                                   const Obs_container_t * const observables ) const;


   
      
};

#endif
