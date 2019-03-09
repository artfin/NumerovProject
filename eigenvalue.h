//
// Created by artfin on 08.03.19.
//
#pragma once

class Eigenvalue
{
public:
    Eigenvalue( int n, double value ) : n(n), value(value)
    {
    }

    int get_n() const { return n; }
    double get_value() const { return value; }

private:
    int n;
    double value;
};


