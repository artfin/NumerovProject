#ifndef GENERALIZEDMATRIXNUMEROV_H
#define GENERALIZEDMATRIXNUMEROV_H

#include "parameters.h"
#include <Eigen/Dense>
#include <cmath>

class GeneralizedMatrixNumerov
{
public:
    GeneralizedMatrixNumerov( Parameters * parameters );

    void allocateMatrices();
    void fillMatrices();

private:
    Parameters * parameters;

    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
    Eigen::MatrixXd V;
    Eigen::MatrixXd H;
};

#endif // GENERALIZEDMATRIXNUMEROV_H
