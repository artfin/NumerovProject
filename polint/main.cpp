#include <iostream>
#include <cmath>
#include <vector>

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

int main()
{
    std::cout << "Neville algorithm." << std::endl;

    int i, j;
    int n = 4;
    double X; // X will be used to approximate f(X)
    std::vector<double> x;
    return 0;
}