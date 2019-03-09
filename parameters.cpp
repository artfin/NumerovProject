#include "parameters.h"

void Parameters::setPotential( std::function<double(const double)> potential )
{
    this->potential = potential;
}

double Parameters::brent( std::function<double(double)> f, double xb1, double xb2, const double eps )
{
    int status;
    int iter = 0, max_iter = 100;

    const gsl_root_fsolver_type *T;
    T = gsl_root_fsolver_brent;

    gsl_root_fsolver *s;
    s = gsl_root_fsolver_alloc(T);

    gsl_function F =
            {
                    [](double d, void *vf) -> double {
                        auto &func = *static_cast<std::function<double(double)> *>(vf);
                        return func(d);
                    },
                    &f
            };

    gsl_root_fsolver_set(s, &F, xb1, xb2);

    double x_lo, x_hi, r = 0;

    do {
        iter++;

        gsl_root_fsolver_iterate(s);
        r = gsl_root_fsolver_root(s);

        x_lo = gsl_root_fsolver_x_lower(s);
        x_hi = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(x_lo, x_hi, eps, eps);

        /*
        if ( status == GSL_SUCCESS ) {
            printf( "Converged:\n");
        }

        printf( "%5d [%.12f, %.12f] %.12f %.12f\n",
                iter, x_lo, x_hi, r, x_hi - x_lo );
        */

    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);

    return r;
}

void Parameters::findTurningPoints()
{
    auto f = [=](double x) { return potential(x) - maxEnergy; };

#ifdef DEBUG_FIND_TURNING_POINTS
    std::cout << std::setprecision(12);
    std::cout << std::endl;
    std::cout << " CALLING BRENT TO FIND TURNING POINTS " << std::endl;
    std::cout << " INTERVAL FOR FINDING LEFT TURNING POINT: (" << leftPointLeftBound << ", " << rightPointLeftBound << ")" << std::endl;
#endif

    leftTurningPoint = brent(f, leftPointLeftBound, leftPointRightBound, epsilon);

#ifdef DEBUG_FIND_TURNING_POINTS
    std::cout << " LEFT TURNING POINT FOUND: " << leftTurningPoint << std::endl;
#endif

#ifdef DEBUG_FIND_TURNING_POINTS
    std::cout << " INTERVAL FOR FINDING RIGHT TURNING POINT: (" << leftPointRightBound << ", " << rightPointRightBound << ")" << std::endl;
#endif

    rightTurningPoint = brent(f, rightPointLeftBound, rightPointRightBound, epsilon);

#ifdef DEBUG_FIND_TURNING_POINTS
    std::cout << " RIGHT TURNING POINT FOUND: " << rightTurningPoint << std::endl;
    std::cout << std::endl;
#endif
}

void Parameters::setGridParameters()
{
    ENERGY_BASED_GRID = true;

    lambda = 2.0 * M_PI / std::sqrt(2.0 * mu * maxEnergy); // minimum De'Broglie length
    d =  lambda / (2.0 * M_PI); // hbar * hbar / sqrt(2 * mu * E_max)
    N = static_cast<int>( (rightTurningPoint - leftTurningPoint) / d + 8.0 * M_PI );

    // recalculating grid step after we set the number of point in grid
    d = (rightTurningPoint - leftTurningPoint + 4.0 * lambda) / (N - 1);
}

void Parameters::setGridParameters( const double leftTurningPoint, const double rightTurningPoint, const double d )
{
    FIXED_GRID = true;
    this->leftTurningPoint = leftTurningPoint;
    this->rightTurningPoint = rightTurningPoint;
    this->d = d;
    N = static_cast<int>( (rightTurningPoint - leftTurningPoint) / d  );
}
