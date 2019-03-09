#pragma once

#include "parameters.h"
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

class MatrixNumerov
{    
public:
    MatrixNumerov( Parameters * parameters );

    void allocateMatrices();
    void fillBasicNumerovMatrices();
    Eigen::VectorXd diagonalizeHamiltonian();

    void calculatePartitionFunction( Eigen::VectorXd & eigenvalues );

private:
    Parameters * parameters;

    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
    Eigen::MatrixXd V;

    Eigen::MatrixXd H;
};

