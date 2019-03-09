#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <vector>
#include <stdexcept>

class MatrixReader
{
public:
    explicit MatrixReader( std::string const & filename );
    void loadFormatFile( std::string const & filename, std::vector<std::string> & content );
    void parseFormatFile( std::vector<std::string> & content );
    void parseLine( std::string & line, bool & is_empty, std::string & variable, std::string & value, int lineNumber );

    inline bool isInteger(const std::string & s);
    int string_to_int(std::string const & value, size_t line);

    void fillMatrix( Eigen::MatrixXd & m, double d );

    void resetFile( std::string const & filename );

    int get_factor() const { return factor; }
    int get_d_power() const { return d_power; }
    std::vector<std::pair<int, double>> const& get_mtxFmt() const { return mtxFmt; }

private:
    int factor; // множитель перед матрицей
    int d_power = 0; // степень d перед матрицей
    std::vector<std::pair<int, double>> mtxFmt; // вектор пар для заполнения матрицы
};
