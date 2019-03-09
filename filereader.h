#pragma once
#include "parameters.h"

#include <string>
#include <fstream>
#include <stdexcept>
#include <cstdio> // std::sscanf

#define DEBUG_SET_PARAMETERS
#undef DEBUG_SET_PARAMETERS

class FileReader
{
public:
    FileReader(std::string filename, Parameters * parameters);

    void loadContent( std::vector<std::string> & content );
    void parseContent( std::vector<std::string> & content );
    void parseLine( std::string line, bool & is_assignment, std::string & variable, std::string & value, size_t lineNumber );

    double string_to_double( std::string const & value, size_t line );
    int string_to_int( std::string const & value, size_t line );

private:
    std::string filename;
    Parameters * parameters;
};
