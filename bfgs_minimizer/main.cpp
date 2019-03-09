#include <iostream>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include <functional>
#include <iomanip>
#include <cmath>

double brent( std::function<double(double)> f, double xb1, double xb2, const double eps )
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

        if ( status == GSL_SUCCESS ) {
            printf( "Converged:\n");
        }

        printf( "%5d [%.12f, %.12f] %.12f %.12f\n",
        iter, x_lo, x_hi, r, x_hi - x_lo );
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);

    return r;
}

int main()
{
    double b = 2.0;
    auto f = [=](double x) { return x*x - b; };

    double lb = 1.0;
    double rb = 5.0;
    double eps = 1.0e-12;

    double r = brent(f, lb, rb, eps);
    std::cout << std::setprecision(12);
    std::cout << "(main) r: " << r << std::endl;
    std::cout << "(main) b: " << b << std::endl;
    std::cout << "(main) sqrt(b): " << std::sqrt(b) << std::endl;

    return 0;
}
