#ifndef ASHKIN_TELLER2D_PARAMETER_STRING 
#define ASHKIN_TELLER2D_PARAMETER_STRING 

#include "askin-teller2d_parameters.cxx"
#include <string>

struct Askin_Teller2d_Parameter_String
{
    const std::string L = std::to_string(Askin_Teller2d_Parameters::L);
    const std::string N = std::to_string(Askin_Teller2d_Parameters::N);
    const std::string J = std::to_string(Askin_Teller2d_Parameters::J);
    const std::string h = std::to_string(Askin_Teller2d_Parameters::h);
    const std::string num_neighbors = std::to_string(Askin_Teller2d_Parameters::num_neighbors_i);

    const std::string energy_min = std::to_string(Askin_Teller2d_Parameters::energy_min);
    const std::string energy_max = std::to_string(Askin_Teller2d_Parameters::energy_max);
    const std::string energy_bin_size = std::to_string(Askin_Teller2d_Parameters::energy_bin_size);
    const std::string num_bins = std::to_string(Askin_Teller2d_Parameters::num_bins);

    const std::string model_name = "Askin_Teller2d";
    const std::string size_string = "L-" + L;

    std::string file_name_base;
    std::string file_header;

    Askin_Teller2d_Parameter_String();

    ~Askin_Teller2d_Parameter_String(){}
};

#endif
