#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>

#include "../constants.h"

struct Level
{
public:
    Level( double value, int j ) : value(value), j(j) { }

    double value;
    int j;
};


void parse_file( std::string const & filename, std::vector<Level> & levels )
{
    std::ifstream ifs( filename );
    std::string line;

    double value;
    int current_j = 0;

    while( getline(ifs, line) )
    {
        if ( line.empty() )
            continue;

        if ( line.find('#') != std::string::npos )
            continue;

        if ( line.find("J =") != std::string::npos )
        {
            size_t len = line.size();
            size_t pos = line.find( "J =" ) + 3;
            std::string j_str = line.substr( pos, len - pos );
            current_j = std::stoi( j_str );
        }
        else
        {
            std::stringstream ss( line );
            ss >> value;
            Level lvl( value, current_j );
            levels.push_back( lvl );
        }
    }
}

double CALCPARTSUM( const double Temperature, const std::vector<Level> & lvls )
// Calculates ro-vibrational partition sum at a given Temperature
{
    double partsum = 0.0;

    for ( auto const& lvl : lvls )
    {
        // wavenumber*PLANCKCONST*SPLIGHTCM = energy in Joules
        double exp_ = std::exp( - lvl.value * constants::PLANCKCONST * constants::SPLIGHTCM / constants::BOLTZCONST / Temperature );
        partsum += (2.0*lvl.j+1) * exp_;
    }

    return partsum;
}

int main()
{
    std::string filename = "../levels_ord12.txt";
    std::vector<Level> levels;

    parse_file( filename, levels );
    std::cout << "Levels parsed: " << levels.size() << std::endl;

    int total_levels = 0;
    for ( auto const& level : levels )
        total_levels += (2 * level.j + 1);
    std::cout << "Total levels: " << total_levels << std::endl;

    Level max_level = levels[0];
    for ( auto const& level : levels ) {
        if ( level.value > max_level.value ) {
            max_level = level;
        }
    }
    std::cout << "Max level: " << max_level.value << "; j: " << max_level.j <<
                 "; degeneracy: " << (2.0*max_level.j + 1) << std::endl;

    double LTEMP = 500000.0;
    double HTEMP = 500000.0;
    double STEP = 10.0;

    //std::ofstream ofs( "../quantum_pf.txt" );
    //ofs << std::fixed << std::setprecision(10);

    std::cout << std::fixed << std::setprecision(10);
    for ( double TEMP = LTEMP; TEMP <= HTEMP; TEMP += STEP )
    {
        double partsum = CALCPARTSUM( TEMP, levels );
        std::cout << TEMP << "   " << partsum << std::endl;
    }

    return 0;
}
