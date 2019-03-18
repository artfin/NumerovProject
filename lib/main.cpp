//
// Created by artfin on 17.03.19.
//

#include <iostream>
#include <string>
#include <iomanip>
#include "matrixreader.h"

#include <Eigen/Dense>

int main()
{
    std::string filename = "../10order_central_difference/a.mtx";

    MatrixReader matrixReader( filename );

    const int N = 7;
    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(N, N);

    const double h = 1.0;
    matrixReader.fillMatrix(m, h);

    std::cout << std::fixed << std::setprecision(9);
    std::cout << "matrix: " << std::endl << m << std::endl;

    return 0;
}
