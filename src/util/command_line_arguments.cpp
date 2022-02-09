#include "command_line_arguments.hpp"
#include <iostream>

CMDLine_Parser::CMDLine_Parser( const int argc, char * argv [] )
{
    invalid_arguments = false;
    cmd_line_words.resize(argc);
    for (int idx = 0; idx != argc; ++idx)
    {
        cmd_line_words[idx] = std::string(argv[idx]);
    }

    parsed_arguments.resize(MAX_ARGS);
    parsed_arguments[program_name] = std::string(argv[program_name]);
    parsed_arguments[data_location] = data_location_prefix + default_data_location;
    parsed_arguments[job_id] = job_id_prefix + default_jobid;

    if (argc > MAX_ARGS)
    {
        std::cerr << "\nError: Too many command line arguments supplied (" << argc << "). This code uses at most " << MAX_ARGS << "\n.";
        invalid_arguments = true;
    }
}

void CMDLine_Parser::parse_through_arguments()
{
    parse_data_location();
    parse_jobid();
    arguments_already_parsed = true;
}

CMDLine_Parser::CMDLine_Parser::PairIndex CMDLine_Parser::search_for_prefix( const std::string & prefix ) const 
{
    std::string::size_type found_index = std::string::npos;
    std::string::size_type prefix_index = std::string::npos;

    for (std::string::size_type idx = 1; idx != cmd_line_words.size(); ++idx)
    {
        if (found_index == std::string::npos)
        {
            found_index = cmd_line_words[idx].find(data_location_prefix);
            prefix_index = idx;
        }
    }

    // If found_index still is npos, then set prefix_index = std::string::npos as a flag
    // to use the default arguments

    return CMDLine_Parser::PairIndex{ prefix_index, found_index };
}

void CMDLine_Parser::parse_data_location()
{    
    CMDLine_Parser::PairIndex prefix_data_loc = search_for_prefix( data_location_prefix );

    if (prefix_data_loc.index2 == std::string::npos)
    {
        std::cout << "\nData Location Prefix (" << data_location_prefix << ") not supplied in any argument. Defaulting to " << default_data_location <<"\n";
        std::string substring = parsed_arguments[data_location].substr( data_location_prefix.size() );
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
    
    CMDLine_Parser::PairIndex prefix_jobid = search_for_prefix( job_id_prefix );

    if (prefix_jobid.index2 == std::string::npos)
    {
        std::cerr << "\nError: Job ID Prefix (" << job_id_prefix << ") not supplied. No defaults are allowed.\n";
        invalid_arguments = true;
        return;
    }

    parsed_arguments[job_id] = cmd_line_words[prefix_jobid.index1].substr(prefix_jobid.index2);
    invalid_arguments = false;
    return;
}

const std::string & CMDLine_Parser::get_data_location() const
{
    return parsed_arguments[data_location];
}

const std::string & CMDLine_Parser::get_job_id() const
{
    return parsed_arguments[job_id];
}

bool CMDLine_Parser::arguments_invalid() 
{
    if (!arguments_already_parsed)
        parse_through_arguments();

    return invalid_arguments;
}