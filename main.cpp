#include <iostream>
#include <cmath>
#include <chrono>
#include <string>

#include "parameters.h"
#include "filereader.h"
#include "matrixnumerov.h"
#include "generalizedmatrixnumerov.h"
#include "eigenvalue.h"

#define DEBUG_SHOW_PARAMETERS
#undef DEBUG_SHOW_PARAMETERS

namespace harmonic
{
    const double mu = 1.0;
    const double omega = 1.0;

    double potential( const double x )
    {
        return 0.5 * mu * std::pow(omega, 2) * x * x;
    }
}

namespace morse
{
    const double De = 0.1744;
    const double omega = 1.02764;
    const double re = 1.40201;
    const double mu = 918.6446;

    double potential( const double r )
    {
        return De * std::pow(1 - std::exp(-omega * (r-re)*(r-re)), 2);
    }
}

std::string getDirByOrder( const int ORDER )
{
    return "../" + std::to_string(ORDER) + "order";
}

Eigen::VectorXd computeEigenvalues_arnoldi( GeneralizedMatrixNumerov & generalizedMatrixNumerov, const int n )
{
    generalizedMatrixNumerov.allocateMatrices();
    generalizedMatrixNumerov.fillMatrices();

    Eigen::VectorXd eigvals = Eigen::VectorXd::Zero( n );
    generalizedMatrixNumerov.diagonalize_arnoldi( n, eigvals );

    return eigvals;
}

std::vector<Eigenvalue> computeEigenvalues( GeneralizedMatrixNumerov & generalizedMatrixNumerov, const int n )
// performs MatrixNumerov computation of eigenvalues;
// returns first n eigenvalues
{
    generalizedMatrixNumerov.allocateMatrices();
    generalizedMatrixNumerov.fillMatrices();

    int size_ = generalizedMatrixNumerov.get_N();
    Eigen::VectorXd eigvals = Eigen::VectorXd::Zero( size_ );
    Eigen::MatrixXd eigvecs = Eigen::MatrixXd::Zero( size_, size_ );
    generalizedMatrixNumerov.diagonalize( eigvals, eigvecs );

    std::vector<Eigenvalue> result;
    for ( int k = 0; k < n; ++k )
        result.emplace_back(k, eigvals(k));

    return result;
}

void dump_row( size_t i, std::vector<double> const& R )
{
    std::cout << std::setprecision(15);

    std::cout << "R[" << i << "] = ";
    for(size_t j = 0; j <= i; ++j)
        std::cout << R[j] << " ";

    std::cout << std::endl;
}

double FRE( const double Eh, const double Eh2, const int N )
// first Richardson extrapolation with asymptotic order N
{
    double pt = std::pow(2.0, N);
    return (pt * Eh2 - Eh) / (pt - 1);
}

double SRE( const double Eh, const double Eh2, const double Eh4, const int N )
// second Richardson extrapolation with asymptotic order N
{
   double pt = std::pow(2.0, N);
   double fpt = 5 * pt;
   return Eh4 + ( (fpt-1)*Eh4 - fpt*Eh2 + Eh) / ( 4*pt*pt - fpt + 1 );
}

int main()
{
    Parameters parameters;
    FileReader fileReader( "../parameters.inp", &parameters );

    parameters.set_mu( harmonic::mu );
    parameters.setPotential( harmonic::potential );
    //parameters.findTurningPoints();

    double leftTurningPoint = -10.0;
    double rightTurningPoint = 10.0;
    double h = 0.25;
    parameters.setGridParameters( leftTurningPoint, rightTurningPoint, h );

#ifdef DEBUG_SHOW_PARAMETERS
    int size_ = parameters.get_N();
    std::cout << std::setprecision(12);
    std::cout << " PARAMETERS USED INSIDE MATRIX NUMEROV ALGORITHM" << std::endl;
    std::cout << " -----------------------------------------------------" << std::endl;
    std::cout << " REDUCED MASS MU = " << parameters.get_mu() << std::endl;
//    std::cout << " UPPER BOUND OF EIGENVALUES TO BE FIND = " << parameters.get_maxEnergy() << std::endl;
//    std::cout << " MIN DE BROGLIE WAVELENGTH = " << parameters.get_lambda() << std::endl;
    std::cout << " GRID SIZE N = " << size_ << std::endl;
    std::cout << " GRID STEP d = " << parameters.get_d() << std::endl;
    std::cout << " -----------------------------------------------------" << std::endl << std::endl;
#endif

    const int ORDER = 4;
    std::string dir = getDirByOrder(ORDER);

    std::cout << " USING MATRIX NUMEROV METHOD OF THE ORDER = " << ORDER << std::endl;

    GeneralizedMatrixNumerov generalizedMatrixNumerov( &parameters, dir );

    int N = 10; // number of eigenvalues to compute
    int niter = 3; // number of iterations to reduce the step

    std::vector<std::vector<Eigenvalue>> eigenvalues;
    for ( int iter = 0; iter < niter; ++iter )
    {
        parameters.setGridParameters( leftTurningPoint, rightTurningPoint, h );
        std::vector<Eigenvalue> eigs = computeEigenvalues( generalizedMatrixNumerov, N );
        eigenvalues.push_back( eigs );
        h /= 2.0;
    }

    int nlevel = 9;
    std::vector<Eigenvalue> level;
    for ( auto const& eig : eigenvalues )
        level.push_back( eig[nlevel] );

    double exact = nlevel + 0.5;

    double Eh = level[0].get_value();
    double Eh2 = level[1].get_value();
    double Eh4 = level[2].get_value();
    std::cout << std::setprecision(15);
    std::cout << "Eh: " << Eh << std::endl;
    std::cout << "Eh2: " << Eh2 << std::endl;
    std::cout << "Eh4: " << Eh4 << std::endl;

    const int REORDER = ORDER - 2;
    double FREV = FRE( Eh, Eh2, REORDER );
    double SREV = SRE( Eh, Eh2, Eh4, REORDER );

    std::cout << "FRE: " << FREV << std::endl;
    std::cout << "SRE: " << SREV << std::endl;

    double D1 = std::abs(exact - FREV);
    double D2 = std::abs(exact - SREV);
    std::cout << "Diff with FREV: " << D1 << std::endl;
    std::cout << "Diff with SREV: " << D2 << std::endl;
    std::cout << "Diff with smallest h level: " << std::abs(exact - level[2].get_value()) << std::endl;

    // AAdHP correction
    double mukn = 12.0 * (1 - std::cos(nlevel*h)) / (h*h*(5.0 + std::cos(nlevel*h)));
    double corrected = Eh + nlevel*nlevel - mukn;
    double mukn4 = 12.0 * (1 - std::cos(nlevel*h/4.0)) / (h/4.0*h/4.0*(5.0 + std::cos(nlevel*h/4.0)));
    double corrected_h4 = Eh4 + nlevel*nlevel - mukn4;

    std::cout << "Eh corrected: " << corrected << std::endl;
    std::cout << "Diff corrected: " << std::abs(exact - corrected) << std::endl;
    std::cout << "Eh4 corrected: " << corrected_h4 << std::endl;
    std::cout << "Diff corrected h4: " << std::abs(exact - corrected_h4) << std::endl;

    return 0;
}



//auto start = std::chrono::high_resolution_clock::now();
//auto end = std::chrono::high_resolution_clock::now();
//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( end - start ).count() / 1000.0;

