#include "matrixnumerov.h"

MatrixNumerov::MatrixNumerov(Parameters * parameters)
    : parameters(parameters)
{
}

void MatrixNumerov::allocateMatrices()
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

void MatrixNumerov::fillBasicNumerovMatrices()
{   
    int size = A.rows();

    double value = 1.0 / parameters->get_d() / parameters->get_d();
    for ( int i = 0; i < size; i++ )
    {
        if ( i - 1 >= 0 )
            A(i, i - 1) = value;
        if ( i + 1 < size )
            A(i, i + 1) = value;

        A(i, i) = -2 * value;
    }

    value = 1.0 / 12.0;
    for ( int i = 0; i < size; i++ )
    {
        if ( i - 1 >= 0 )
            B(i, i - 1) = value;
        if ( i + 1 < size )
            B(i, i + 1) = value;

        B(i, i) = 10.0 / 12.0;
    }

    double x0 = - parameters->get_d() * (parameters->get_N() - 1) / 2.0;
    for ( int i = 0; i < size; i++ )
    {
        double x = x0 + parameters->get_d() * i;
        V(i, i) = parameters->potential(x);
    }

    H = - B.inverse() * A / 2.0 + V;
}

Eigen::VectorXd MatrixNumerov::diagonalizeHamiltonian()
{
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > eigensolver( H );

    if ( eigensolver.info() != Eigen::Success )
        abort();

    std::cout << "Number \t Eigenvalue" << std::endl;
    std::cout << "-----------------------------------" << std::endl;
    for ( int i = 0; i < 50; i++ )
    {
        std::cout << std::fixed << std::setprecision(10);
        std::cout << i << " " << eigensolver.eigenvalues()[i] << std::endl;
    }

    Eigen::VectorXd eigenvalues;
    eigenvalues.resize( H.rows() );
    for ( int i = 0; i < H.rows(); i++ )
        eigenvalues(i) = eigensolver.eigenvalues()[i];

    return eigenvalues;
}

void MatrixNumerov::calculatePartitionFunction(Eigen::VectorXd & eigenvalues)
{

}



