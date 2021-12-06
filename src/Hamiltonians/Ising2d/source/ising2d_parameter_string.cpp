#include "../ising2d_parameter_string.hpp"
#include <iostream>

#if JOB_ARRAYS
Ising2d_Parameter_String::Ising2d_Parameter_String( const std::string & job_id_string )
#else
Ising2d_Parameter_String::Ising2d_Parameter_String()
#endif
{
#if JOB_ARRAYS
    job_id = job_id_string;
#endif
    update_file_header();
}

void Ising2d_Parameter_String::update_file_header()
{
    file_name_base = "";
#if JOB_ARRAYS
    file_name_base += "JOBID-" + job_id + "_"; 
#endif
    file_name_base += "REWL_L-" + L + "_J-" + J + "_h-" + h;
#if PHONON_MEDIATED_NEMATIC_INTERACTIONS
    file_name_base += "_g-" + PMNI_Coupling;
#endif // PHONON_MEDIATED_NEMATIC_INTERACTIONS

    file_name_base += ".dat";

    std::cout << "\nCurrent model name: " << model_name << "\n";
    std::cout << "\nCurrent model name: " << model_name << "\n";
 
    file_header = "# REWL " + model_name + " on a Periodic Square Lattice\n#";
    file_header += "\n# Hamiltonian Parameters";
    file_header += "\n#    L = " + L;
    file_header += "\n#    N = " + N;
    file_header += "\n#    num_DoF = " + num_DoF;
    file_header += "\n#    J = " + J;
    file_header += "\n#    h = " + h;
#if PHONON_MEDIATED_NEMATIC_INTERACTIONS
    file_header += "\n#    PMNI = " + PMNI_Coupling;
#endif // PHONON_MEDIATED_NEMATIC_INTERACTIONS
    file_header += "\n#    num_neighbors = " + num_neighbors;
    file_header += "\n#";
    file_header += "\n# Wang Landau Parameters";
    file_header += "\n#    energy_min = " + energy_min;
    file_header += "\n#    energy_max = " + energy_max;
    file_header += "\n#    bin_size = " + energy_bin_size;
    file_header += "\n#    num_bins = " + num_bins;
}
