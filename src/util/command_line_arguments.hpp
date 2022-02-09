#ifndef _COMMAND_LINE_ARGUMENTS_H_
#define _COMMAND_LINE_ARGUMENTS_H_
// This header defines a class that operates on the
// command line arguments passed to main().
//
// This should be pretty lightweight considering
// right now (Feb 8, 2022), I only care about JOB_ARRAYS
// and any non-default data directories.
//
// The arguments should be passed as the following
//
//   --data_location=path/to/data/location
//
//   --job_id=jobid

#include <vector>
#include <string>
#include <filesystem>  // Requires C++17

#if JOB_ARRAYS
    static const bool SEARCH_FOR_JOBID = true;
#else  /* JOB_ARRAYS */
    static const bool SEARCH_FOR_JOBID = false;
#endif /* JOB_ARRAYS */

struct CMDLine_Parser
{
    // State boolean to keep track of whether the
    // arguments are appropriate or not.
    bool invalid_arguments = false;
    bool arguments_already_parsed = false;
    const bool search_for_jobid = SEARCH_FOR_JOBID;

    const std::string default_data_location = std::filesystem::current_path().parent_path().string();  // This will by default place the data directory node outside of the build directory
    const std::string default_jobid = "None";
    
    const std::string data_location_prefix = "--data_location=";
    const std::string job_id_prefix = "--job_id=";
    
    std::vector<std::string> cmd_line_words;
    std::vector<std::string> parsed_arguments;

    enum CMDLine_Args
    {
        program_name=0, data_location, job_id, MAX_ARGS
    };

    CMDLine_Parser( const int argc, const char * argv [] );


    void parse_through_arguments( );
    void parse_data_location();
    void parse_jobid();

    const std::string & get_data_location() const;
    const std::string & get_job_id() const;

    bool valid_arguments();
    
    virtual ~CMDLine_Parser() {}

    struct PairIndex
    {
        std::string::size_type index1;
        std::string::size_type index2;
    };
    
    PairIndex search_for_prefix( const std::string & prefix ) const;
};

#endif /* _COMMAND_LINE_ARGUMENTS_H_ */