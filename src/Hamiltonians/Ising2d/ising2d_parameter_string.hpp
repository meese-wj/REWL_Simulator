#ifndef ISING2D_PARAMETER_STRING
#define ISING2D_PARAMETER_STRING

#include "ising2d_parameters.cxx"
#include <string>

struct Ising2d_Parameter_String
{
    const std::string L = std::to_string(Ising2d_Parameters::L);
    const std::string N = std::to_string(Ising2d_Parameters::N);
    const std::string num_DoF = std::to_string(Ising2d_Parameters::num_DoF);
    const std::string J = std::to_string(Ising2d_Parameters::J);
    const std::string h = std::to_string(Ising2d_Parameters::h);
    const std::string num_neighbors = std::to_string(Ising2d_Parameters::num_neighbors_i);

    std::string energy_min = std::to_string(Ising2d_Parameters::energy_min);
    std::string energy_max = std::to_string(Ising2d_Parameters::energy_max);
    const std::string energy_bin_size = std::to_string(Ising2d_Parameters::energy_bin_size);
    std::string num_bins = std::to_string(Ising2d_Parameters::num_bins);

#if RFIM
    const std::string model_name = "Ising2d_RFIM";
#else
    const std::string model_name = "Ising2d";
#endif
    const std::string size_string = "L-" + L;

    std::string file_name_base;
    std::string file_header;

#if JOB_ARRAYS
    std::string job_id;
    Ising2d_Parameter_String( const std::string & job_id_string );
#else
    Ising2d_Parameter_String();
#endif

    ~Ising2d_Parameter_String(){}

    void update_file_header();
};

#endif
