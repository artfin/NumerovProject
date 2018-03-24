#include "generalizedmatrixnumerov.h"

GeneralizedMatrixNumerov::GeneralizedMatrixNumerov( Parameters * parameters )
    : parameters(parameters)
{

}

void GeneralizedMatrixNumerov::allocateMatrices()
{
    int size = parameters->get_N();

    A.resize(size, size);
    B.resize(size, size);
    V.resize(size, size);
    H.resize(size, size);

    for ( int i = 0; i < size; i++ )
    {
        for ( int j = 0; j < size; j++ )
        {
            A(i, j) = 0.0;
            B(i, j) = 0.0;
            V(i, j) = 0.0;
            H(i, j) = 0.0;
        }
    }
}

void GeneralizedMatrixNumerov::fillMatrices()
{
    int size = parameters->get_N();

    Eigen::MatrixXd sys;
    sys.resize(size, size);

    for ( int i = 0; i < size; i++ )
    {
        for ( int j = 0; j < size; j++ )
        {
            sys(i, j) = std::pow(i * parameters->get_d(), 2 * (j + 1));
        }
    }
}

