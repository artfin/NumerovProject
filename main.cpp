#include "parameters.h"
#include "filereader.h"
#include "matrixnumerov.h"

#include <iostream>
#include <cmath>
#include <ctime>

// возьмем потенциал f(x) = |x|
double potential( double x )
{
    return 0.5 * x * x; // std::abs(x);
}

int main()
{
    Parameters parameters;
    FileReader fileReader( "/home/artfin/Desktop/CPPcode/NumerovProject/parameters.inp", &parameters );

    parameters.show();

    parameters.setPotential( potential );
    parameters.findTurningPoints();
    parameters.setGridParameters();

    MatrixNumerov matrixNumerov( &parameters );

    matrixNumerov.allocateMatrices();
    matrixNumerov.fillBasicNumerovMatrices();

    std::clock_t start = std::clock();
    Eigen::VectorXd eigenvalues = matrixNumerov.diagonalizeHamiltonian();
    std::cout << "Diagonalization done in: " << (double) (std::clock() - start) / CLOCKS_PER_SEC << " s" << std::endl;

    matrixNumerov.calculatePartitionFunction( eigenvalues );

    return 0;
}
