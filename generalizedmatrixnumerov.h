#ifndef GENERALIZEDMATRIXNUMEROV_H
#define GENERALIZEDMATRIXNUMEROV_H

#include "parameters.h"
#include <Eigen/Dense>
#include <cmath>
#include <string>
#include "matrixreader.h"

class GeneralizedMatrixNumerov
{
public:
    GeneralizedMatrixNumerov( Parameters * parameters, std::string const & dir );

    void allocateMatrices( );
    void fillMatrices( );
    Eigen::VectorXd diagonalizeHamiltonian( );

private:
    Parameters * parameters;
    std::string dir;

    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
    Eigen::MatrixXd V;
    Eigen::MatrixXd H;
};

#endif // GENERALIZEDMATRIXNUMEROV_H
