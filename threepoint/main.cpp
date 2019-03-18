#include <iostream>
#include <cmath>
#include <chrono>
#include <string>

#include "../lib/parameters.h"
#include "../lib/filereader.h"
#include "../lib/generalizedmatrixnumerov.h"
#include "../lib/eigenvalue.h"

#define DEBUG_SHOW_PARAMETERS
//#undef DEBUG_SHOW_PARAMETERS

namespace harmonic
{
    const double mu = 1.0;
    const double omega = 1.0;

    double potential( const double x )
    {
        return 0.5 * mu * std::pow(omega, 2) * x * x;
    }
}

namespace exponential
{
    const double mu = 0.5;
    double potential( const double x )
    {
        return std::exp(x);
    }
}

#define EPSILON 1.0E-15
bool AreSame(double a, double b)
{
    return std::fabs(a - b) < EPSILON;
}

namespace nopotential
{
    const double mu = 0.5;
    double potential( const double x )
    {
        if ( x == 0 ) return std::numeric_limits<double>::max();
        else if ( AreSame(x, M_PI) ) return std::numeric_limits<double>::max();
        return 0.0;
    }
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

double FRE( const double Eh, const double Eh2, const int N )
// first Richardson extrapolation with asymptotic order N
{
    double pt = std::pow(2.0, N);
    return (pt * Eh2 - Eh) / (pt - 1);
}

std::vector<double> FREV( std::vector<Eigenvalue> const& Eh, std::vector<Eigenvalue> const& Eh2, const int N )
// performing first Richardson extrapolation with asymptotic order N with vector of energies
{
    size_t size_ = Eh.size();
    std::vector<double> extrapolated;
    extrapolated.reserve( size_ );

    for ( size_t k = 0; k < size_; ++k )
        extrapolated.push_back( FRE(Eh[k].get_value(), Eh2[k].get_value(), N) );

    return extrapolated;
}

double SRE( const double Eh, const double Eh2, const double Eh4, const int N )
// second Richardson extrapolation with asymptotic order N
{
    double pt = std::pow(2.0, N);
    double fpt = 5 * pt;
    return Eh4 + ( (fpt-1)*Eh4 - fpt*Eh2 + Eh) / ( 4*pt*pt - fpt + 1 );
}

std::vector<double> SREV( std::vector<Eigenvalue> const& Eh, std::vector<Eigenvalue> const& Eh2, std::vector<Eigenvalue> const& Eh4, const int N )
// performing second Richardson extrapolation with asymptotic order N with vector of energies
{
    size_t size_ = Eh.size();
    std::vector<double> extrapolated;
    extrapolated.reserve( size_ );

    for ( size_t k = 0; k < size_; ++k )
        extrapolated.push_back( SRE(Eh[k].get_value(), Eh2[k].get_value(), Eh4[k].get_value(), N) );

    return extrapolated;
}

int main()
{
    Parameters parameters;

    // Test parameters for Harmonic oscillator
    /*
    parameters.set_mu( harmonic::mu );
    parameters.setPotential( harmonic::potential );
    double leftTurningPoint = -10.0;
    double rightTurningPoint = 10.0;
    double h = 0.125;
    parameters.setGridParameters( leftTurningPoint, rightTurningPoint, h );
    */

    // Test parameters for Exponential potential
    parameters.set_mu( nopotential::mu );
    parameters.setPotential( nopotential::potential );
    double leftTurningPoint = 0.0;
    double rightTurningPoint = M_PI;
    double L = rightTurningPoint - leftTurningPoint;
    int n = 7;
    double h = L / n;
    std::cout << "h: " << h << std::endl;
    parameters.setGridParameters( leftTurningPoint, rightTurningPoint, n );

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

    //std::string dir = "../../2order_central_difference/";
    //const int ORDER = 12;
    //std::string dir = "../../" + std::to_string(ORDER) + "order_central_difference/";
    //std::string dir = "../../12order_central_difference/";
    //std::string dir = "../../pade_3_3/";
    //std::string dir = "../../4order_central_difference/";
    std::string dir = "../../8order_corrected/";

    const int REORDER = 2;
    std::cout << " USING MATRIX 2ND ORDER CENTRAL DIFFERENCE METHOD" << std::endl;

    GeneralizedMatrixNumerov generalizedMatrixNumerov( &parameters, dir );

    int N = 5; // number of eigenvalues to compute
    //std::vector<Eigenvalue> eigs = computeEigenvalues( generalizedMatrixNumerov, N );
    std::vector<Eigenvalue> eigs = computeEigenvalues( generalizedMatrixNumerov, N );

    Eigen::MatrixXd const& A = generalizedMatrixNumerov.get_A();
    std::cout << "A: " << std::endl << A << std::endl;

    std::vector<double> exact = { 1.0, 4.0, 9.0, 16.0, 25.0 };
    //std::vector<double> exact = { 11.5424, 41.1867, 90.5404, 159.6296, 248.4569 };

    //std::ofstream ofs("../12order_diff_k.txt");
    for ( int k = 0; k < n - 1; ++k ) {
        std::cout << (k+1)*(k+1) << " " << eigs[k].get_value() << std::endl;
        //std::cout << k+1 << " " << std::abs(eigs[k].get_value() - (k+1)*(k+1)) << std::endl;
    }

    //std::ofstream ofs("../" + std::to_string(ORDER) + "order_central_diff_h.txt");

    std::ofstream ofs("../8order_corrected_diff_h.txt");
    ofs << std::fixed << std::setprecision(15);

    for ( int k = 10; k < 400; k += 10 ) {
        parameters.setGridParameters( leftTurningPoint, rightTurningPoint, k );
        std::vector<Eigenvalue> eigs = computeEigenvalues( generalizedMatrixNumerov, N );
        double diff = std::abs(eigs[0].get_value() - exact[0]);
        double h = L / k;

        std::cout << h << " " << diff << std::endl;
        ofs << h << " " << diff << std::endl;
    }

    /*
    parameters.setGridParameters( leftTurningPoint, rightTurningPoint, 2*n );
    std::vector<Eigenvalue> eigs2 = computeEigenvalues( generalizedMatrixNumerov, N );

    parameters.setGridParameters( leftTurningPoint, rightTurningPoint, 4*n );
    std::vector<Eigenvalue> eigs3 = computeEigenvalues( generalizedMatrixNumerov, N );

    std::vector<double> fre = FREV( eigs2, eigs3, REORDER );
    std::vector<double> sre = SREV( eigs, eigs2, eigs3, REORDER );

    std::cout << std::fixed << std::setprecision(15);
    */
    /*
    for ( size_t k = 0; k < eigs.size(); ++k ) {
        std::cout << std::abs(eigs[k].get_value() -exact[k]) << " " << std::abs(eigs2[k].get_value()-exact[k]) << " " <<
                    std::abs(fre[k] - exact[k]) << std::endl;
    }
    */
    /*
    std::cout << std::scientific << std::setprecision(3);
    for ( size_t k = 0; k < eigs.size(); ++k ) {
        std::cout << std::abs(eigs[k].get_value() - exact[k]) << " "
                  << std::abs(eigs2[k].get_value() - exact[k]) << " "
                  << std::abs(eigs3[k].get_value() - exact[k]) << " "
                  << std::abs(fre[k] - exact[k]) << " "
                  << std::abs(sre[k] - exact[k]) << std::endl;
    }
    */

    return 0;
}

