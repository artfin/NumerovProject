#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <map>

class MatrixReader
{
public:
    explicit MatrixReader( std::string const & filename );
    void loadFormatFile( );
    void parseContent( );
    void parseLine( std::string & line, bool & is_empty, std::string & variable, std::string & value, size_t lineNumber );

    std::pair<int, int> parseTwoIntegers( const std::string & s, size_t line );

    inline bool isTwoIntegers( const std::string &s );
    inline bool isInteger(const std::string & s);
    int string_to_int(std::string const & value, size_t line);

    void fillMatrix( Eigen::MatrixXd & m, double d );

    void resetFile( std::string const & filename );

    int get_factor() const { return factor; }
    int get_d_power() const { return d_power; }
    std::vector<std::pair<int, double>> const& get_mtxFmt() const { return mtxFmt; }

private:
    std::string filename;
    std::vector<std::string> content;

    int factor; // множитель перед матрицей
    int d_power = 0; // степень d перед матрицей
    std::vector<std::pair<int, double>> mtxFmt;

    std::map<std::pair<int, int>, double> extreme;
};
