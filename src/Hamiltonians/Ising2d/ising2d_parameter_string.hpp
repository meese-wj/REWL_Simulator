#ifndef ISING2D_PARAMETER_STRING
#define ISING2D_PARAMETER_STRING

#include "ising2d_parameters.cxx"
#include <string>

static std::string change_Ising_model_name( std::string model_name )
{
    std::string output_name = model_name;
#if PHONON_MEDIATED_NEMATIC_INTERACTIONS
    output_name += "_PMNI";
#endif // PHONON_MEDIATED_NEMATIC_INTERACTIONS
#if RFIM
    output_name += "_RFIM";
#endif // RFIM
    return output_name;
}

struct Ising2d_Parameter_String
{
    const std::string L = std::to_string(Ising2d_Parameters::L);
    const std::string N = std::to_string(Ising2d_Parameters::N);
    const std::string num_DoF = std::to_string(Ising2d_Parameters::num_DoF);
    const std::string J = std::to_string(Ising2d_Parameters::J);
    const std::string h = std::to_string(Ising2d_Parameters::h);
#if PHONON_MEDIATED_NEMATIC_INTERACTIONS
    const std::string PMNI_Coupling = std::to_string(Ising2d_Parameters::PMNI_Coupling);
#endif // PHONON_MEDIATED_NEMATIC_INTERACTIONS
    const std::string num_neighbors = std::to_string(Ising2d_Parameters::num_neighbors_i);

    std::string energy_min = std::to_string(Ising2d_Parameters::energy_min);
    std::string energy_max = std::to_string(Ising2d_Parameters::energy_max);
    std::string energy_bin_size = std::to_string(Ising2d_Parameters::energy_bin_size);
    std::string num_bins = std::to_string(Ising2d_Parameters::num_bins);

    const std::string model_name = change_Ising_model_name( "Ising2d" );
    const std::string size_string = "L-" + L;

    std::string file_name_base;
    std::string file_header;

#if JOB_ARRAYS
    std::string job_id;
    Ising2d_Parameter_String( const std::string & job_id_string );
#else
    Ising2d_Parameter_String();
#endif // JOB_ARRAYS

    ~Ising2d_Parameter_String(){}

    void update_file_header();
};

#endif
