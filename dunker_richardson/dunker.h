//
// Created by artfin on 10.03.19.
//

#pragma once

#include <cmath>
#include <gsl/gsl_sf_legendre.h>

#include "../constants.h"

class Dunker
{
public:
    void setJ( const int J ) { this->J = J; }
    int getJ() const { return J; }

    static double potential( const double r )
    // returns cm-1
    {
        double Rang = r * constants::BOHRTOANG;

        double mult1 = -epsilon * alpha / (alpha - 6.0) * std::pow(Rm / Rang, 6.0);
        double mult2 = epsilon * 6.0 / (alpha - 6.0) * std::exp(alpha * (1.0 - Rang / Rm));

        // mu is in fact MU(a.m.u)/HTOCM, hence centrifugal term is in cm-1
        return (mult1 * P0A + mult2 * P0R) + J*(J+1.0)/(2.0 * mu * r * r);
    }

    static int J;

    static const double mu;
    static const double MU;

    static const double epsilon;
    static const double alpha;
    static const double Rm;

    static const double P0A;
    static const double P0R;
};
