#include "parameters.h"
#include "filereader.h"
#include "matrixnumerov.h"
#include "generalizedmatrixnumerov.h"
#include "matrixreader.h"

#include <iostream>
#include <cmath>
#include <ctime>

#include <string>
#include <limits.h>
#include <unistd.h>

double potential( double x )
{
    return 0.5 * x * x;
}

std::string getApplicationBinPath()
{
  char result[ PATH_MAX ];
  ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
  return std::string( result, (count > 0) ? count : 0 );
}

int main()
{
    Parameters parameters;
    FileReader fileReader( "../parameters.inp", &parameters );

    parameters.show();

    parameters.setPotential( potential );
    parameters.findTurningPoints();
    parameters.setGridParameters();

    /*
    MatrixNumerov matrixNumerov( &parameters );

    matrixNumerov.allocateMatrices();
    matrixNumerov.fillBasicNumerovMatrices();

    std::clock_t start = std::clock();
    Eigen::VectorXd eigenvalues = matrixNumerov.diagonalizeHamiltonian();
    std::cout << "Diagonalization done in: " << (double) (std::clock() - start) / CLOCKS_PER_SEC << " s" << std::endl;

    matrixNumerov.calculatePartitionFunction( eigenvalues );
    */

    std::string dir = "../14order";
    std::cout << std::endl << "--- Using " << dir << " Numerov Method --- " << std::endl << std::endl;

    GeneralizedMatrixNumerov generalizedMatrixNumerov( &parameters, dir );
    generalizedMatrixNumerov.allocateMatrices();
    generalizedMatrixNumerov.fillMatrices();

    std::clock_t start = std::clock();
    Eigen::VectorXd eigenvalues = generalizedMatrixNumerov.diagonalizeHamiltonian();
    std::cout << "Diagonalization done in: " << (double) (std::clock() - start) / CLOCKS_PER_SEC << " s" << std::endl;

    return 0;
}
