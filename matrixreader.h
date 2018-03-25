#ifndef MATRIXREADER_H
#define MATRIXREADER_H

#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <vector>
#include <stdexcept>

class MatrixReader
{
public:
    MatrixReader( std::string const & filename );
    void loadFormatFile( std::string const & filename, std::vector<std::string> & content );
    void parseFormatFile( std::vector<std::string> & content );
    void parseLine( std::string & line, bool & is_empty, std::string & variable, std::string & value, int lineNumber );

    inline bool isInteger(const std::string & s);
    int string_to_int(std::string const & value, size_t line);

    void showMatrixFormat();

    void fillMatrix( Eigen::MatrixXd & m, double d );

    void resetFile( std::string const & filename );

private:
    int factor; // множитель перед матрицей
    int d_power = 0; // степень d перед матрицей
    std::vector<std::pair<int, double>> mtxFmt; // вектор пар для заполнения матрицы
};

#endif // MATRIXREADER_H
