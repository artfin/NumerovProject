//
// Created by artfin on 10.03.19.
//

#pragma once

#include <cmath>
#include "../constants.h"

// Danila test morse potential
class Morse
{
public:
    Morse() = default;
    void setJ( const int J ) { this->J = J; }
    int getJ() const { return J; }
    static double potential( const double r ) {
        return De*(1.0-std::exp(-a*(r-re)))*(1.0-std::exp(-a*(r-re))) + J*(J+1.0)/(2.0*mu*r*r) - De;
    }

    static int J;
    static const double De;
    static const double a;
    static const double re;
    static const double mu;
    static const double MU;
    static const double nu0;
    static const double omega0;
};