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

std::string getDirByOrder( const int ORDER )
{
    return "../" + std::to_string(ORDER) + "order";
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
    parameters.set_mu( harmonic::mu );
    parameters.setPotential( harmonic::potential );
    double leftTurningPoint = -10.0;
    double rightTurningPoint = 10.0;
    double L = rightTurningPoint - leftTurningPoint;
    int n = 40;
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

    const int ORDER = 4;
    std::string dir = getDirByOrder(ORDER);
    //std::string dir = "../" + std::to_string(ORDER) + "order_central_difference/";
    //std::string dir = "../" + std::to_string(ORDER) + "order_corrected/";
    //std::string dir = "../pade_2_2/";

    //std::cout << " USING MATRIX NUMEROV METHOD OF THE ORDER = " << ORDER << std::endl;

    GeneralizedMatrixNumerov generalizedMatrixNumerov( &parameters, dir );

    /*
    int N = 10; // number of eigenvalues to compute
    int niter = 3; // number of iterations to reduce the step

    std::vector<std::vector<Eigenvalue>> eigenvalues;
    for ( int iter = 0; iter < niter; ++iter, n *= 2 )
    {
        parameters.setGridParameters( leftTurningPoint, rightTurningPoint, n );
        std::vector<Eigenvalue> eigs = computeEigenvalues( generalizedMatrixNumerov, N );
        eigenvalues.push_back( eigs );
    }

    const int REORDER = ORDER - 2;

    std::vector<double> first_extrapolation = FREV( eigenvalues[1], eigenvalues[2], REORDER );
    std::vector<double> second_extrapolation = SREV( eigenvalues[0], eigenvalues[1], eigenvalues[2], REORDER );

    std::cout << std::setprecision(15);
    std::cout << "h/4 \t 1st extrapolation \t 2nd extrapolation" << std::endl;
    for ( int k = 0; k < N; ++k ) {
        double exact = k + 0.5;
        std::cout << std::abs(exact - eigenvalues.back()[k].get_value()) << " \t " <<
                     std::abs(exact - first_extrapolation[k]) << " \t " <<
                     std::abs(exact - second_extrapolation[k]) << std::endl;
    }
    */

    //std::string filename = "./pade_2_2_diff_k.txt";
    //std::string filename = "./" + std::to_string(ORDER) + "order_numerov_diff_h.txt";
    //std::string filename = "./" + std::to_string(ORDER) + "order_central_diff_h.txt";
    //std::string filename = "./" + std::to_string(ORDER) + "order_corrected_diff_h.txt";
    /*
    std::ofstream ofs(filename);
    ofs << std::fixed << std::setprecision(15);

    for ( int k = 40; k < 400; k += 10 ) {
        parameters.setGridParameters( leftTurningPoint, rightTurningPoint, k );
        std::vector<Eigenvalue> eigs = computeEigenvalues( generalizedMatrixNumerov, 1 );
        double diff = std::abs(eigs[0].get_value() - 0.5);
        double h = L / k;

        ofs << h << " " << diff << std::endl;
    }
     */

    //std::string filename = "./plot_k/" + std::to_string(ORDER) + "order_diff_k3.txt";
    //std::string filename = "./pade_1_1_diff_k.txt";
    std::string filename = "./" + std::to_string(ORDER) + "order_central_difference_diff_k.txt";

    std::ofstream ofs(filename);
    ofs << std::fixed << std::setprecision(15);

    parameters.setGridParameters( leftTurningPoint, rightTurningPoint, n );
    std::vector<Eigenvalue> eigs = computeEigenvalues( generalizedMatrixNumerov, n - 1 );

    /*
    Eigen::MatrixXd const& A = generalizedMatrixNumerov.get_A();
    std::cout << "A: " << std::endl << A << std::endl;

    Eigen::MatrixXd const& B = generalizedMatrixNumerov.get_B();
    std::cout << "B: " << std::endl << B << std::endl;
    */

    /*
    for ( size_t k = 0; k < eigs.size(); ++k ) {
        ofs << k+1 << " " << std::abs(eigs[k].get_value() - (k+0.5)) << std::endl;
    }
     */


    return 0;
}

/*
// AAdHP correction
double mukn = 12.0 * (1 - std::cos(nlevel*h)) / (h*h*(5.0 + std::cos(nlevel*h)));
double corrected = Eh + nlevel*nlevel - mukn;
double mukn4 = 12.0 * (1 - std::cos(nlevel*h/4.0)) / (h/4.0*h/4.0*(5.0 + std::cos(nlevel*h/4.0)));
double corrected_h4 = Eh4 + nlevel*nlevel - mukn4;
*/

//auto start = std::chrono::high_resolution_clock::now();
//auto end = std::chrono::high_resolution_clock::now();
//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( end - start ).count() / 1000.0;

