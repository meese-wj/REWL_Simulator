#ifndef SIMULATED_ANNEALING_PARAMETER_STRINGS
#define SIMULATED_ANNEALING_PARAMETER_STRINGS

#include <string>
#include "simulated_annealing_parameters.cxx"

struct SA_Strings
{
    const std::string initial_temperature  = std::to_string(SA_Parameters::initial_temperature);
    const std::string final_temperature    = std::to_string(SA_Parameters::final_temperature);
    const std::string num_iterations       = std::to_string(SA_Parameters::num_iterations);
    const std::string block_size           = std::to_string(SA_Parameters::block_size);
    const std::string energy_var_tolerance = std::to_string(SA_Parameters::energy_stdev_tolerance);
    const std::string frozen_criterion      = std::to_string(SA_Parameters::frozen_criterion);

    std::string header = "";

    SA_Strings();

    ~SA_Strings(){}

};

SA_Strings::SA_Strings()
{
    header  = "# Simulated Annealing Parameters:\n";
    header += "#\n";
    header += "#     Initial Temperature       = " + initial_temperature + "\n";
    header += "#     Final Temperature         = " + final_temperature + "\n";
    header += "#     Number of Iterations      = " + num_iterations + "\n";
    header += "#     Block Size                = " + block_size + "\n";
    header += "#     Energy Variance Tolerance = " + energy_var_tolerance + "\n";
    header += "#     Frozen Criterion          = " + frozen_criterion + "\n";
    header += "#";
}

#endif
