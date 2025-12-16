#include <ode.h>
#include <math.h>
#include <string.h>

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
    int X = t > tstart && t <= tend; 

    double Xn = pow(X, n);
    double Zn = pow(z[0], n);

    dzdt[0] = b * (Xn / (KXZn + Xn)) * (KZn / (KZn + Zn)) - a * z[0];
    return GSL_SUCCESS;
}

int solve(gsl_odeiv2_driver* d, int max_time, double y[], size_t nTraits, int par_id, double measure_interval)
{
    double t = 0;
    int status = 0;
    double* p = (double*)d->sys->params;
    double tstart = p[0];
    double tend = p[1];
    int X = 0;

    int max_iter = max_time / measure_interval;
    for (int i = 0; i < max_iter; ++i)
    {
        double ti = i * measure_interval;
        status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        
        if (status != GSL_SUCCESS)
        {
            fprintf(stderr, "Unable to apply fixed step to stepper function %s - error code %d\n", d->s->type->name, status);
            break;
        }
        
        X = t >= tstart && t < tend;
        // print solution to std out
        for (size_t j = 0; j < nTraits-1; ++j)
        {
            printf("%d, %.2lf, %d, %.3lf,", par_id, ti, X, y[j]);
        }
        printf("%d, %.2lf, %d, %.3lf\n", par_id, ti, X, y[nTraits-1]);
    }
    return GSL_SUCCESS;
}


const gsl_odeiv2_step_type* get_stepper_from_input(char* input_string)
{
    // Check the input is fine
    if (input_string == NULL)
    {
        return NULL;
    }

    if (strcmp("rk4", input_string) == 0)
    {
        return gsl_odeiv2_step_rk4;
    }

    if (strcmp("rkf45", input_string) == 0)
    {
        return gsl_odeiv2_step_rkf45;
    }

    if (strcmp("msadams", input_string) == 0)
    {
        return gsl_odeiv2_step_msadams;
    }

    if (strcmp("rk4-fs", input_string) == 0)
    {
        return gsl_odeiv2_step_rk4;
    }

    return NULL;
}