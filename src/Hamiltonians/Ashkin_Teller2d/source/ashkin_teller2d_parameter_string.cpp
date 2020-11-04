#include "../ashkin_teller2d_parameter_string.hpp"

#if JOB_ARRAYS
Ashkin_Teller2d_Parameter_String::Ashkin_Teller2d_Parameter_String( const std::string & job_id_string )
#else
Ashkin_Teller2d_Parameter_String::Ashkin_Teller2d_Parameter_String()
#endif
{
#if JOB_ARRAYS
    file_name_base = "JOBID-" + job_id_string + "_REWL_L-" + L + "_J-" + J + "_K-" + K + ".dat";
#else
    file_name_base = "REWL_L-" + L + "_J-" + J + "_K-" + K + ".dat";
#endif

    file_header = "# REWL " + model_name + " on a Periodic Square Lattice\n#";
    file_header += "\n# Hamiltonian Parameters";
    file_header += "\n#    L = " + L;
    file_header += "\n#    N = " + N;
    file_header += "\n#    num_DoF = " + num_DoF;
    file_header += "\n#    J = " + J;
    file_header += "\n#    K = " + K;
    file_header += "\n#    num_neighbors = " + num_neighbors;
    file_header += "\n#";
    file_header += "\n# Wang Landau Parameters";
    file_header += "\n#    energy_min = " + energy_min;
    file_header += "\n#    energy_max = " + energy_max;
    file_header += "\n#    bin_size = " + energy_bin_size;
    file_header += "\n#    num_bins = " + num_bins;
}
