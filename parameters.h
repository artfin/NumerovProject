#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <iostream>
#include <vector>
#include <functional>
#include <random>
#include <iomanip>

class Parameters
{
public:
    explicit Parameters();
    Parameters( double maxEnergy );

    void show();
    void setPotential(std::function<double(double)> potential);
    void findTurningPoints();
    void setGridParameters();

    void set_d( double d ) { this->d = d; }
    void set_mass( double mass ) { this->mass = mass; }
    void set_maxEnergy(double maxEnergy) { this->maxEnergy = maxEnergy; }
    void set_epsilon(double epsilon) { this->epsilon = epsilon; }
    void set_upperBound(double upperBound ) { this->upperBound = upperBound; }
    void set_lowerBound(double lowerBound ) { this->lowerBound = lowerBound; }

    double getLowerBound( ) const { return lowerBound; }
    double getUpperBound( ) const { return upperBound; }

    double get_leftTurningPoint() const { return leftTurningPoint; }
    double get_rightTurningPoint() const { return rightTurningPoint; }

    double get_d() const { return d; }
    int get_N() const { return N; }
    double get_mass() const { return mass; }

    std::function<double(double)> potential;

private:
    double maxEnergy; // максимальная энергия

    double mass; // масса тела

    // параметры поиска поворотных точек
    double epsilon;
    double lowerBound;
    double upperBound;

    // левая и правая разворотные точки
    double rightTurningPoint;
    double leftTurningPoint;

    // параметры сетки
    double d; // расстояние между точками сетки
    int N; // количество точек
};

#endif // PARAMETERS_H
