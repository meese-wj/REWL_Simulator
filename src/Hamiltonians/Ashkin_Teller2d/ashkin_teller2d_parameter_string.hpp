#ifndef ASHKIN_TELLER2D_PARAMETER_STRING 
#define ASHKIN_TELLER2D_PARAMETER_STRING 

#include "ashkin_teller2d_parameters.cxx"
#include <string>

struct Ashkin_Teller2d_Parameter_String
{
    const std::string L = std::to_string(Ashkin_Teller2d_Parameters::L);
    const std::string N = std::to_string(Ashkin_Teller2d_Parameters::N);
    const std::string J = std::to_string(Ashkin_Teller2d_Parameters::J);
    const std::string K = std::to_string(Ashkin_Teller2d_Parameters::K);
    const std::string num_neighbors = std::to_string(Ashkin_Teller2d_Parameters::num_neighbors_i);

    const std::string energy_min = std::to_string(Ashkin_Teller2d_Parameters::energy_min);
    const std::string energy_max = std::to_string(Ashkin_Teller2d_Parameters::energy_max);
    const std::string energy_bin_size = std::to_string(Ashkin_Teller2d_Parameters::energy_bin_size);
    const std::string num_bins = std::to_string(Ashkin_Teller2d_Parameters::num_bins);

    const std::string model_name = "Ashkin_Teller2d";
    const std::string size_string = "L-" + L;

    std::string file_name_base;
    std::string file_header;

    Ashkin_Teller2d_Parameter_String();

    ~Ashkin_Teller2d_Parameter_String(){}
};

#endif
