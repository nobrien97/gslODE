#include <ode.h>
#include <math.h>

int NARODE (double t, const double z[], double dzdt[], void* params)
{
    // Params is length 7
    double* p = (double*)params;

    // Extract parameters for clarity
    double tstart = p[0];
    double tend = p[1];
    double b = p[2];

    // Get n to do the power calculations now
    double n = p[5];

    double KXZn = pow(p[3], n); // Raise to power of n
    double KZn = pow(p[4], n);
    double a = p[6];
    int X = t >= tstart && t < tend; 

    double Xn = pow(X, n);
    double Zn = pow(z[0], n);

    dzdt[0] = b * (Xn / (KXZn + Xn)) * (KZn / (KZn + Zn)) - a * z[0];
    return GSL_SUCCESS;
}

int solve(gsl_odeiv2_driver* d, int maxIter, double y[], size_t nTraits)
{
    double t = 0;

    for (int i = 0; i < maxIter; ++i)
    {
        double ti = i * maxIter / 100.0;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);

        if (status != GSL_SUCCESS)
        {
            fprintf(stderr, "error, return value = %d\n", status);
            break;
        }
        // print solution to std out
        for (size_t j = 0; j < nTraits-1; ++j)
        {
            printf("%.2lf, %.3lf,", ti, y[i]);
        }
        printf("%.2lf, %.3lf\n", ti, y[nTraits-1]);
    }
    return GSL_SUCCESS;
}