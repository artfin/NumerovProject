#pragma once

#include <iostream>
#include <vector>
#include <functional>
#include <random>
#include <iomanip>

#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#define DEBUG_FIND_TURNING_POINTS
//#undef DEBUG_FIND_TURNING_POINTS

class Parameters
{
public:
    explicit Parameters() = default;

    void setPotential(std::function<double(const double)> potential);

    double brent( std::function<double(double)> f, double xb1, double xb2, double eps );
    void findTurningPoints();
    void setGridParameters();
    void setGridParameters( double leftTurningPoint, double rightTurningPoint, double h );

    void set_d( double d ) { this->d = d; }
    void set_mu( double mu ) { this->mu = mu; }
    void set_maxEnergy(double maxEnergy) { this->maxEnergy = maxEnergy; }
    void set_epsilon(double epsilon) { this->epsilon = epsilon; }

    void set_leftPointleftBound( double value ) { leftPointLeftBound = value; }
    void set_leftPointrightBound( double value ) { leftPointRightBound = value; }
    void set_rightPointleftBound( double value ) { rightPointLeftBound = value; }
    void set_rightPointrightBound( double value ) { rightPointRightBound = value; }

    double get_leftTurningPoint() const { return leftTurningPoint; }
    double get_rightTurningPoint() const { return rightTurningPoint; }

    double get_lambda() const { return lambda; }
    double get_d() const { return d; }
    int get_N() const { return N; }
    double get_mu() const { return mu; }
    double get_maxEnergy() const { return maxEnergy; }

    double call_potential( const double x ) { return potential(x); }

    bool FIXED_GRID = false;
    bool ENERGY_BASED_GRID = false;

private:
    double maxEnergy; // максимальная энергия
    std::function<double(const double)> potential;

    double mu; // масса тела

    // параметры поиска поворотных точек
    double epsilon;

    double leftPointLeftBound;
    double leftPointRightBound;
    double rightPointLeftBound;
    double rightPointRightBound;

    // левая и правая разворотные точки
    double rightTurningPoint;
    double leftTurningPoint;

    // параметры сетки
    double lambda;
    double d; // расстояние между точками сетки
    int N; // количество точек
};
