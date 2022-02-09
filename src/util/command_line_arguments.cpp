#include "command_line_arguments.hpp"
#include <iostream>

CMDLine_Parser::CMDLine_Parser( const int argc, const char * argv [] )
{
    invalid_arguments = false;
    cmd_line_words.resize(argc);
    for (int idx = 0; idx != argc; ++idx)
    {
        cmd_line_words[idx] = std::string(argv[idx]);
    }

    parse_arguments.resize(MAX_ARGS);
    parse_arguments[program_name] = std::string(argv[program_name]);
    parse_arguments[data_location] = data_location_prefix + default_data_location;
    parse_arguments[job_id] = job_id_prefix + default_jobid;

    if (argc > MAX_ARGS)
    {
        std::cerr << "\nError: Too many command line arguments supplied (" << argc << "). This code uses at most " << MAX_ARGS << "\n.";
        invalid_arguments = true;
    }
}

void CMDLine_Parser::parse_arguments()
{
    parse_data_location();
    parse_jobid();
    arguments_already_parsed = true;
}

CMDLine_Parser::PairIndex CMDLine_Parser::search_for_prefix( const std::string & prefix ) const 
{
    int found_index = std::npos;
    int prefix_index = 1;

    for (int idx = 1; idx != cmd_line_words.size(); ++idx)
    {
        if (found_index == std::npos)
        {
            found_index = cmd_line_words[idx].find(data_location_prefix);
            prefix_index = idx;
        }
    }

    // If found_index still is npos, then set prefix_index = std::npos as a flag
    // to use the default arguments
    if (found_index == std::npos)
        prefix_index = std::npos;

    return PairIndex( prefix_index, data_loc_index );
}

void CMDLine_Parser::parse_data_location()
{    
    PairIndex prefix_data_loc = search_for_prefix( data_location_prefix );

    if (prefix_data_loc.index2 == std::npos)
    {
        std::cout << "\nData Location Prefix (" << data_location_prefix << ") not supplied in any argument. Defaulting to " << default_data_location <<"\n";
        std::string substring = parsed_arguments[data_location].substring( data_location_prefix.size() );
        parsed_arguments[data_location] = substring;
        return;
    }
     
    parsed_arguments[data_location] = cmd_line_words[prefix_data_loc.index1].substr(prefix_data_loc.index2 + data_location_prefix.size());
    invalid_arguments = false;
    return;
}

void CMDLine_Parser::parse_jobid()
{
    if (!search_for_jobid)
        return;
    
    PairIndex prefix_jobid = search_for_prefix( job_id_prefix );

    if (prefix_jobid.index2 == std::npos)
    {
        std::cerr << "\nError: Job ID Prefix (" << job_id_prefix << ") not supplied. No defaults are allowed.\n";
        invalid_arguments = true;
        return;
    }

    parse_arguments[job_id] = cmd_line_words[prefix_jobid.index1].substr(prefix_jobid.index2);
    invalid_arguments = false;
    return;
}

const std::string & CMDLine_Parser::get_data_location() const
{
    return parse_arguments[data_location];
}

const std::string & CMDLine_Parser::get_job_id() const
{
    return parse_arguments[job_id];
}

bool CMDLine_Parser::valid_arguments() const 
{
    if (!arguments_already_parsed)
        parse_arguments();

    return !(invalid_arguments);
}