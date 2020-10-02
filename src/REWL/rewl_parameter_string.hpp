#ifndef REWL_PARAMETER_STRING
#define REWL_PARAMETER_STRING
/* Create the parameter string for
 * the REWL parameters. */
#include <string>
#include "rewl_parameters.hpp"

struct REWL_Parameter_String
{
    const std::string num_walkers         = std::to_string(REWL_Parameters::num_walkers);
    const std::string replicas_per_window = std::to_string(REWL_Parameters::replicas_per_window);
    const std::string window_overlap      = std::to_string(REWL_Parameters::window_overlap);

    const std::string sweeps_per_check    = std::to_string(REWL_Parameters::sweeps_per_check);
    const std::string final_increment     = std::to_string(REWL_Parameters::final_increment);
    const std::string flatness_criterion  = std::to_string(REWL_Parameters::flatness_criterion);

    std::string file_header;

    REWL_Parameter_String();

    ~REWL_Parameter_String(){}
};

REWL_Parameter_String::REWL_Parameter_String()
{
    file_header  = "# REWL Parameters";
    file_header += "\n#    num_walkers = " + num_walkers; 
    file_header += "\n#    replicas_per_window = " + replicas_per_window; 
    file_header += "\n#    window_overlap = " + window_overlap; 
    file_header += "\n#    sweeps_per_check = " + sweeps_per_check; 
    file_header += "\n#    final_increment = " + final_increment; 
    file_header += "\n#    flatness_criterion = " + flatness_criterion; 
}

#endif
