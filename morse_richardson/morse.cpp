//
// Created by artfin on 10.03.19.
//

#include "morse.h"

int Morse::J = 0;
const double Morse::De = 200.0; // cm-1
const double Morse::a = 0.5;
const double Morse::re = 6.5; // bohr
const double Morse::mu = 0.17; // in what units??
const double Morse::nu0 = a / (2.0 * M_PI) * std::sqrt( 2.0 * De / mu );
const double Morse::omega0 = 2.0 * M_PI * nu0;
