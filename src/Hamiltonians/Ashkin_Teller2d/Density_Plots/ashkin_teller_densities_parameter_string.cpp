#ifndef AT_DENSITIES_PARAMETER_STRINGS
#define AT_DENSITIES_PARAMETER_STRINGS

#include "ashkin_teller_densities_parameters.cxx"

AT_Density_Parameters::Strings::Strings()
{
    header  = "# Ashkin-Teller 2d Density Plot\n";
    header += "#\n";
    header += "# Density Plot Parameters:\n";
    header += "#\n";
    header += "#     axis bins    = " + axis_bins + "\n";
    header += "#     total bins   = " + total_bins + "\n";
    header += "#     axis maximum = " + sigma_tau_max + "\n";
    header += "#     axis minimum = " + sigma_tau_min + "\n";
    header += "#     binwidth     = " + density_binwidth + "\n";
    header += "#";
}

#endif
