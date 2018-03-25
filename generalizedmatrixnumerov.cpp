#include "generalizedmatrixnumerov.h"

GeneralizedMatrixNumerov::GeneralizedMatrixNumerov( Parameters * parameters, std::string const & dir )
    : parameters(parameters), dir(dir)
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
   MatrixReader matrixReader(dir + "/a.mtx");
   matrixReader.fillMatrix(A, parameters->get_d());
   //std::cout << "A: " << std::endl << A << std::endl;

   matrixReader.resetFile(dir + "/b.mtx");
   matrixReader.fillMatrix(B, parameters->get_d());
   //std::cout << "B: " << std::endl << B << std::endl;

   double x0 = - parameters->get_d() * (parameters->get_N() - 1) / 2.0;
   for ( int i = 0; i < V.rows(); i++ )
   {
       double x = x0 + parameters->get_d() * i;
       V(i, i) = parameters->potential(x);
   }

   H = - B.inverse() * A / (2.0 * parameters->get_mass()) + V;
}

Eigen::VectorXd GeneralizedMatrixNumerov::diagonalizeHamiltonian()
{
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > eigensolver( H );

    if ( eigensolver.info() != Eigen::Success )
        abort();

    std::cout << "Number \t Eigenvalue" << std::endl;
    std::cout << "-----------------------------------" << std::endl;
    std::cout << std::fixed << std::setprecision(10);
    for ( int i = 0; i < 10; i++ )
        std::cout << i << " " << eigensolver.eigenvalues()[i] << std::endl;

    Eigen::VectorXd eigenvalues;
    eigenvalues.resize( H.rows() );
    for ( int i = 0; i < H.rows(); i++ )
        eigenvalues(i) = eigensolver.eigenvalues()[i];

    return eigenvalues;
}

