#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <cassert>

void polint( double xa[], double ya[], int n, double x, double *y, double *dy )
// array xa[1..n], y[1..n], and given a value x, this routine returns a value y, and and an error estimate
// dy. The returned value y = P(x).
{
    int i, m, ns = 1;
    double den, dif, dift, ho, hp, w;

    dif = std::abs( x - xa[1] );
    std::vector<double> c, d;
    c.resize(n); d.resize(n);

    for ( i = 1; i <= n; i++ )
    {
        if ( (dift = std::abs(x - xa[i])) < dif ) {
            ns = i;
            dif = dift;
        }
    }

    *y = ya[ns--];

    for ( m = 1; m < n; m++ ) {
        for ( i = 1; i <= n-m; i++ ) {
            ho = xa[i] - x;
            hp = xa[i+m] - x;
            w = c[i+1] - d[i];

            // two xa's are (to within roundoff) identical
            if ( (den = ho-hp) == 0.0 ) {
                std::cerr << "Error in routine polint." << std::endl;
                exit( 1 );
            }

            den = w /den;
            d[i] = hp * den;
            c[i] = ho * den;
        }

        *y += (*dy = (2*ns < (n - m) ? c[ns + 1] : d[ns--]));
    }
}

double neville( std::vector<double> const& x, std::vector<double> const& y, const double X )
// polynomial interpolation using Neville algorithm
{
    assert( x.size() == y.size() );

    size_t size_ = x.size();
    double Q[size_][size_];

    for ( int i = 0; i < size_; i++ )
        Q[i][0] = y[i];

    for ( int i = 1; i < size_; i++ ) {
        for ( int j = 1; j <= i; j++ ) {
            Q[i][j] =(((X - x[i - j]) * Q[i][j-1]) - (X - x[i]) * Q[i - 1][j - 1]) / (x[i] - x[i - j]);
        }
    }

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Resultant Q matrix: " << std::endl;
    for ( int i = 0; i < size_; i++ ) {
        for ( int j = 0; j <= i; j++ ) {
            std::cout << Q[i][j] << " ";
        }
        std::cout << std::endl;
    }

    return Q[size_ - 1][size_ - 1];
}

int main()
{
    std::cout << "Neville algorithm." << std::endl;

    int n = 4;
    double Q[n][n];

    double X = 8.4; // we want to approximate f(X)
    std::vector<double> x{ 8.1, 8.3, 8.6, 8.7 };
    std::vector<double> y{ 16.94410, 17.56492, 18.50515, 18.82091 };


    double interp = neville( x, y, X );
    std::cout << "Interpolated value: " << interp << std::endl;

    return 0;
}