#include <iostream>
#include <functional>
#include <cmath>
#include <vector>
#include <iomanip>

double f( double x )
{
    return std::log(x);
}
void dump_row(size_t i, double *R)
{
   std::cout << "R[" << i << "] = ";
   std::cout << std::setprecision(15);
   for(size_t j = 0; j <= i; ++j)
   {
       std::cout << R[j] << " ";
   }
   printf("\n");
}

double romberg( std::function<double(double)> f, double a, double b, size_t max_steps, double acc )
{
    std::vector<double> R1, R2;
    R1.resize(max_steps); R2.resize(max_steps);

    double *Rp = &R1[0], *Rc = &R2[0]; //Rp is previous row, Rc is current row
    double h = (b-a); //step size
    Rp[0] = (f(a) + f(b))*h*.5; //first trapezoidal step

    dump_row(0, Rp);

    for(size_t i = 1; i < max_steps; ++i)
    {
        h /= 2.;
        double c = 0;
        size_t ep = 1 << (i-1); //2^(n-1)
        for(size_t j = 1; j <= ep; ++j){
            c += f(a+(2*j-1)*h);
        }
        Rc[0] = h*c + .5*Rp[0]; //R(i,0)

        for(size_t j = 1; j <= i; ++j){
            double n_k = pow(4, j);
            Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1); //compute R(i,j)
        }

        //Dump ith column of R, R[i,i] is the best estimate so far
        dump_row(i, Rc);

        if(i > 1 && fabs(Rp[i-1]-Rc[i]) < acc)
        {
            return Rc[i-1];
        }

        //swap Rn and Rc as we only need the last row
        double *rt = Rp;
        Rp = Rc;
        Rc = rt;
    }

    return Rp[max_steps-1]; //return our best guess
}

int main()
{
    double a = 1.0;
    double b = 2.0;
    int n = 4;
    double acc = 1.0e-12;

    romberg(f, a, b, n, acc);

    std::cout << std::setprecision(15);
    std::cout << "Exact: " << 2.0 * std::log(2.0) - 1.0 << std::endl;

    return 0;
}

/*
function r=romberg(f,a,b,n)
h=(b-a)./(2.ˆ(0:n-1));
r(1,1)=(b-a)*(f(a)+f(b))/2;
for j=2:n
subtotal = 0;
for i=1:2ˆ(j-2)
subtotal = subtotal + f(a+(2*i-1)*h(j));
end
r(j,1) = r(j-1,1)/2+h(j)*subtotal;
for k=2:j
r(j,k)=(4ˆ(k-1)*r(j,k-1)-r(j-1,k-1))/(4ˆ(k-1)-1);
end
end
*/