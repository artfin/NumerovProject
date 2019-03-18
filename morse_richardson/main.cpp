#include <iostream>
#include <cmath>
#include <chrono>
#include <string>

#include "../lib/parameters.h"
#include "../lib/filereader.h"
#include "../lib/generalizedmatrixnumerov.h"
#include "../lib/eigenvalue.h"

#include "morse.h"

#define DEBUG_SHOW_PARAMETERS
//#undef DEBUG_SHOW_PARAMETERS

std::string getDirByOrder( const int ORDER )
{
    return "../../" + std::to_string(ORDER) + "order";
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

double morse_exact_nth( const int n )
{
    return Morse::omega0 * (n + 0.5) -  \
            std::pow(Morse::omega0*(n+0.5),2.0) / 4.0 / Morse::De - Morse::De;
}

std::vector<double> morse_exact( )
{
    std::vector<double> levels;
    levels.push_back( morse_exact_nth(0) );
    double level = morse_exact_nth(1);

    for ( int n = 2; levels.back() < level; ++n ) {
        levels.push_back( level );
        level = morse_exact_nth(n);
    }

    return levels;
}

void write_potential( const int J )
{
    std::ofstream ofs("../..//potential_profiles/danila_morse_J=" + std::to_string(J) + ".txt");
    double x = 0.1;
    Morse morse;
    morse.setJ(J);
    for ( int i = 1; i < 200; i++ ) {
        ofs << x << " " << morse.potential(x) << std::endl;
        x += 0.1;
    }
    ofs.close();
}

int main()
{
    Parameters parameters;
    Morse morse;

    const int JMIN = 0;
    const int JMAX = 0;

    morse.setJ( JMIN );
    //std::cout << " ROTATIONAL NUMBER J = " << J << std::endl;

    // Parameters for Danila morse potential
    parameters.set_mu( Morse::mu );
    parameters.setPotential( Morse::potential );
    double leftTurningPoint = 4.0;
    double rightTurningPoint = 40.0;
    double h0 = 0.25;
    double h = h0;
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

    //const int ORDER = 12;
    //std::string dir = getDirByOrder(ORDER);
    std::string dir = "../../pade_3_3/";

    morse.setJ(5);
    //std::cout << " USING MATRIX NUMEROV METHOD OF THE ORDER = " << ORDER << std::endl;

    GeneralizedMatrixNumerov generalizedMatrixNumerov( &parameters, dir );

    int N = 20; // some magical number
    int niter = 3; // number of iterations to reduce the step
    //const int REORDER = ORDER - 2; // order of Richardson extrapolation

    const int REORDER = 12;
    h = h0;

    std::vector<std::vector<Eigenvalue>> eigenvalues;
    for ( int iter = 0; iter < niter; ++iter ) {
        int n = static_cast<int>((rightTurningPoint - leftTurningPoint) / h);
        parameters.setGridParameters( leftTurningPoint, rightTurningPoint, n );
        std::vector<Eigenvalue> eigs = computeEigenvalues( generalizedMatrixNumerov, N );

        eigenvalues.push_back( eigs );

        h /= 2.0;
    }

    std::vector<double> exact = morse_exact();

    std::vector<double> first_extrapolation = FREV( eigenvalues[1], eigenvalues[2], REORDER );
    std::vector<double> second_extrapolation = SREV( eigenvalues[0], eigenvalues[1], eigenvalues[2], REORDER );

    /*
    for( size_t k = 0; k < eigenvalues[0].size(); ++k ) {
        std::cout << std::abs(exact[k] - eigenvalues[0][k].get_value()) << " "
                  << std::abs(exact[k] - eigenvalues[1][k].get_value()) << " "
                  << std::abs(exact[k] - eigenvalues[2][k].get_value()) << " "
                  << std::abs(exact[k] - second_extrapolation[k]) << std::endl;
    }
     */
    /*
    for( size_t k = 0; k < eigenvalues[0].size(); ++k ) {
        std::cout << std::fixed << std::setprecision(8)
                  << exact[k] << " "
                  << second_extrapolation[k] << " "
                  << std::scientific << std::setprecision(3)
                  << std::abs(exact[k] - eigenvalues[2][k].get_value()) << " "
                  << std::min( std::abs(exact[k] - first_extrapolation[k]),
                               std::abs(exact[k] - second_extrapolation[k]) ) << std::endl;
    }
     */

    for ( size_t k = 0; k < eigenvalues[0].size(); ++k ) {
        std::cout << std::fixed << std::setprecision(8)
                  << second_extrapolation[k] << " "
                  << std::scientific << std::setprecision(3)
                  << std::abs(eigenvalues[2][k].get_value() -second_extrapolation[k]) << std::endl;
    }

    /*
    std::ofstream ofs("../levels.txt");
    ofs << "# ORDER OF GMN: " << ORDER << std::endl;
    ofs << "# ORDER OF RICHARDSON EXTRAPOLATION: " << REORDER << std::endl;
    ofs << "# H VALUES: " << h0 << "; " << h0/2.0 << "; " << h0/4.0 << std::endl;
    ofs << "# LEFT RIGID BOUNDARY: " << leftTurningPoint << std::endl;
    ofs << "# RIGHT RIGID BOUNDARY: " << rightTurningPoint << std::endl;
    ofs << std::fixed << std::setprecision(6);

    std::vector<std::vector<Eigenvalue>> eigenvalues;
    std::vector<double> second_extrapolation;
    for ( int J = JMIN; J <= JMAX; ++J ) {
        std::cout << "ROTATIONAL NUMBER J = " << J << std::endl;

        ofs << "J = " << J << std::endl;
        morse.setJ( J );

        // resetting step size
        h = h0;

        for ( int iter = 0; iter < niter; ++iter )
        {
            parameters.setGridParameters( leftTurningPoint, rightTurningPoint, h );
            std::vector<Eigenvalue> eigs = computeEigenvalues( generalizedMatrixNumerov, N );

            eigenvalues.push_back( eigs );
            h /= 2.0;
        }

        second_extrapolation = SREV( eigenvalues[0], eigenvalues[1], eigenvalues[2], REORDER );
        for ( double & level : second_extrapolation ) {
            if ( level < 0.0 ) {
                ofs << level << std::endl;
            }
        }
        ofs << std::endl;

        second_extrapolation.clear();
        eigenvalues.clear();
    }
    */
    /*
    std::cout << std::fixed << std::setprecision(4);
    for ( int n = 0; n < N; ++n ) {
        double eig = eigenvalues.back()[n].get_value();
        if ( eig < 0.0 ) {
            std::cout << eig + Morse::De << "  " << second_extrapolation[n] + Morse::De << std::endl;
        }
    }
     */

    /*
    ofs << "Exact \t 2nd extrapolation \t Delta; h/4 \t Delta; 2nd extrapolation" << std::endl;
    for ( int n = 0; n < N; ++n ) {
        ofs << eigs_exact[n] << "   " <<
               second_extrapolation[n] << "  " <<
               std::scientific <<
               std::abs(eigs_exact[n] - eigenvalues.back()[n].get_value()) << "  " <<
               std::abs(eigs_exact[n] - second_extrapolation[n]) << std::endl;
    }
     */

    return 0;
}
