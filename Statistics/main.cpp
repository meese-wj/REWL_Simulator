/*  This is the main code for 
 *  the statistics routine that will 
 *  average the microcanonical observables
 *  from an ensemble of independent jobs */
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

struct IO_Strings
{
    enum CLA
    {
        executable, input, output
    };

    const fs::path input_path;
    const fs::path output_path;

    IO_Strings( const char * argv [] ) : input_path( argv[ CLA::input ] ), output_path( argv[ CLA::output ] )
    {}
};


int main( const int argc, const char * argv [] )
{

    return 0;
}

