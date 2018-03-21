#ifndef FILEREADER_H
#define FILEREADER_H

#include "parameters.h"

#include <string>
#include <fstream>
#include <stdexcept>
#include <cstdio> // std::sscanf

class FileReader
{
public:
    FileReader(std::string filename, Parameters * parameters);

    void loadContent( std::vector<std::string> & content );
    void parseContent( std::vector<std::string> & content );
    void parseLine( std::string line, bool & is_assignment, std::string & variable, std::string & value, int lineNumber );
    double string_to_double( std::string & value, size_t line );

private:
    std::string filename;
    Parameters * parameters;
};

#endif // FILEREADER_H
