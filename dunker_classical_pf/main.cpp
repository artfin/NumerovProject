#include <iostream>
#include <chrono>
#include <fstream>

#include <hep/mc-mpi.hpp>
#include <boost/math/special_functions.hpp>

#include "../constants.h"
#include "dunker.h"

double integrand_(hep::mc_point<double> const& x, const double Temperature)
{
    double R_new = x.point()[0];
    double R = std::tan(M_PI / 2 * R_new);

    if ( R < 5.0 )
        return 0.0;

    double U_CM = Dunker::potential( R );
    double U = U_CM / constants::HTOCM;

    //if ( U_CM < -0.030794 ) {
    if ( U < 0 ) {
        double U_kT = U * constants::HTOJ / (constants::BOLTZCONST * Temperature);
        double exp_ = std::exp(-U_kT);

        double gamma = boost::math::gamma_p(1.5, -U_kT);
        double jacR = M_PI / 2 * (1 + std::pow(R, 2));

        return R * R * gamma * exp_ * jacR;
    }

    return 0.0;
}


int main( int argc, char * argv[] )
{
    MPI_Init( &argc, &argv );

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    int size;
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    auto start = std::chrono::high_resolution_clock::now();

    double LTEMP = 5000000;
    double HTEMP = 5000000;
    double STEP = 10;

    std::vector<double> temps;
    std::vector<double> pfs;

    for ( double TEMP = LTEMP; TEMP <= HTEMP; TEMP += STEP )
    {
        // creating integrand for current temperature
        auto integrand = std::bind(integrand_, std::placeholders::_1, TEMP);

        // set the verbose vegas callback function
        hep::mpi_vegas_callback<double>( hep::mpi_vegas_verbose_callback<double> );

        int Ncycles = 5;
        int Npoints = static_cast<int>(1e6);

        auto results = hep::mpi_vegas(
                MPI_COMM_WORLD,
                hep::make_integrand<double>(integrand, 1),
                std::vector<std::size_t>(Ncycles, Npoints)
        );

        auto result = hep::accumulate<hep::weighted_with_variance>( results.begin() + 1, results.end() );
        double chi_square_dof = hep::chi_square_dof<hep::weighted_with_variance>( results.begin() + 1, results.end() );

        if ( rank == 0 )
        {
            std::cout << "TEMPERATURE RESULTS: " << TEMP << std::endl;
            std::cout << ">> Cumulative result: " << std::endl;
            std::cout << ">> N = " << result.calls() << "; I = " << result.value() << " +- " << result.error() << std::endl;
            std::cout << ">> chi^2/dof = " << chi_square_dof << std::endl;
        }

        double multiplier = 4.0 * M_PI * std::pow(2.0 * M_PI * Dunker::MU * constants::AMU * constants::BOLTZCONST * TEMP / constants::PLANCKCONST /
                                                  constants::PLANCKCONST, 1.5) * std::pow(constants::ALU, 3.0);

        double ans = result.value() * multiplier;
        temps.push_back( TEMP );
        pfs.push_back( ans );
    }

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "TEMPERATURE: " << temps.back( ) << "; PF: " << pfs.back() << std::endl;
    /*
    std::ofstream ofs( "../classical_pf_to_upper_level.txt" );
    ofs << std::fixed << std::setprecision(12);

    for ( size_t k = 0; k < temps.size(); ++k )
        ofs << temps[k] << " " << pfs[k] << std::endl;
    */

    return 0;
}
