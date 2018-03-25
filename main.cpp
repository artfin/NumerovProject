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

double testPotential( double x )
{
    return 0.5 * x * x;
}

namespace morse
{
    const double De = 0.1744;
    const double omega = 1.02764;
    const double re = 1.40201;
    const double mu = 918.6446;

    double potential( const double r )
    {
        return De * pow(1 - std::exp(-omega * pow(r - re, 2)), 2);
    }
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

    /*
    std::ofstream outFile("h2_morsePotential.txt");
    for ( double x0 = 0.1; x0 < 15.0; x0 += 0.1 )
        outFile << x0 << " " << morse::potential(x0) << std::endl;
    outFile.close();
    */

    parameters.set_mass( morse::mu );
    parameters.setPotential( morse::potential );
    parameters.findTurningPoints();
    parameters.setGridParameters();

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
