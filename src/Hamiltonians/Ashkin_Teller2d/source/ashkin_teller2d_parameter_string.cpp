#include "../ashkin_teller2d_parameter_string.hpp"

Ashkin_Teller2d_Parameter_String::Ashkin_Teller2d_Parameter_String()
{
    file_name_base = "REWL_L-" + L + "_J-" + J + "_K-" + K + ".dat";

    file_header = "# REWL " + model_name + " on a Periodic Square Lattice\n#";
    file_header += "\n# Hamiltonian Parameters";
    file_header += "\n#    L = " + L;
    file_header += "\n#    N = " + N;
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