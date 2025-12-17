#include <ode.h>
#include <math.h>
#include <string.h>

int ODE_NAR(double t, const double z[], double dzdt[], void* params)
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

// From GSL Example 
int ODE_VanDerPol(double t, const double y[], double f[], void* params)
{
    (void)(t); /* avoid unused parameter warning */
    double mu = *(double *)params;
    f[0] = y[1];
    f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
    return GSL_SUCCESS;
}

int ODE_Lorenz(double t, const double y[], double f[], void* params)
{
    // Parameters
    (void)(t);
    double* p = (double*)params;
    double sigma = p[0];
    double rho = p[1];
    double beta = p[2];

    // Equations
    f[0] = sigma * (y[1] - y[0]);
    f[1] = y[0] * (rho - y[2]) - y[1];
    f[2] = y[0] * y[1] - beta * y[2];
}

int ODE_Robertson(double t, const double y[], double f[], void* params)
{
    // Parameters
    (void)(t);
    (void)(params);

    // Equations
    double powY = pow(10,4) * f[1];
    double powYSqr = pow(10,7) * f[1] * f[1];
    f[0] = -0.04 * f[0] + powY * y[2];
    f[1] = 0.04 * f[1] - powY * y[2] - 3 * powYSqr;
    f[2] = 3 * powYSqr;
}

/// @brief Solves an ODE.
/// @param d The ODE driver.
/// @param max_time Maximum time to simulate.
/// @param y Array of initial states.
/// @param n_y Length of y.
/// @param par_id Index of this solution in the matrix of input parameters.
/// @param measure_interval How often to measure during the simulation.
/// @return GSL_SUCCESS if successful, a GSL error code if not.
int solve(gsl_odeiv2_driver* d, int max_time, double y[], size_t n_y, int par_id, double measure_interval)
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
        printf("%d, %.2lf, %d, ", par_id, ti, X);
        for (size_t j = 0; j < n_y-1; ++j)
        {
            printf("%.3lf,", y[j]);
        }
        printf("%.3lf\n", y[n_y-1]);
    }
    return GSL_SUCCESS;
}

/// @brief Chooses a stepper function from command line input option.
/// @param input_string A user input string to match against.
/// @return A gsl_odeiv2_step_type* stepper for the odeiv2 solver.
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

/// @brief Sets a function pointer to the ODE for use by GSL.
/// @param input_string A user input string to match against.
/// @param start_conditions A pointer to an array of starting conditions. 
/// @return 0 if success, -1 if fail.
const int get_ODE_from_input(char* input_string, struct ode* ODE)
{
    // Check the input is fine
    if (input_string == NULL)
    {
        fprintf(stderr, "Invalid ODE input - see -h for available models.\n");
        return 1;
    }

    if (strcmp("NAR", input_string) == 0)
    {
        ODE->ODE_fn_ptr = &ODE_NAR;
        ODE->n_pars = 7;
        ODE->n_y = 1;
        ODE->pars = (double*)malloc(sizeof(double) * ODE->n_pars);
        ODE->y = (double*)malloc(sizeof(double) * ODE->n_y);

        ODE->y[0] = 0.0;
        return 0;
    }   

    if (strcmp("VanDerPol", input_string) == 0)
    {
        ODE->ODE_fn_ptr = &ODE_VanDerPol;
        ODE->n_pars = 1;
        ODE->n_y = 2;
        ODE->pars = (double*)malloc(sizeof(double) * ODE->n_pars);
        ODE->y = (double*)malloc(sizeof(double) * ODE->n_y);

        ODE->y[0] = 1.0;
        ODE->y[1] = 0.0;
        return 0;
    }

    if (strcmp("Lorenz", input_string) == 0)
    {
        ODE->ODE_fn_ptr = &ODE_Lorenz;
        ODE->n_pars = 3;
        ODE->n_y = 3;
        ODE->pars = (double*)malloc(sizeof(double) * ODE->n_pars);
        ODE->y = (double*)malloc(sizeof(double) * ODE->n_y);

        ODE->y[0] = 1.0;
        ODE->y[1] = 1.0;
        ODE->y[2] = 1.0;
        return 0;
    }

    if (strcmp("Robertson", input_string) == 0)
    {
        ODE->ODE_fn_ptr = &ODE_Robertson;
        ODE->n_pars = 0;
        ODE->n_y = 3;
        ODE->pars = (double*)malloc(sizeof(double) * ODE->n_pars);
        ODE->y = (double*)malloc(sizeof(double) * ODE->n_y);

        ODE->y[0] = 1.0;
        ODE->y[1] = 0.0;
        ODE->y[2] = 0.0;
        return 0;
    }

    // Failed to match - return error
    return 1;
}

int free_ode_system(struct ode* ODE)
{
    free(ODE->pars);
    free(ODE->y);
    free(ODE);
    ODE = NULL;

    return 0;
}

int update_ode_pars(struct ode* ODE, double* new_pars, size_t n_new_pars)
{
    if (ODE->n_pars != n_new_pars)
    {
        fprintf(stderr, "New pars length does not match ODE system par length!\n");
        return -1;
    }

    for (int i = 0; i < ODE->n_pars; ++i)
    {
        ODE->pars[i] = new_pars[i];
    }

    return 0;

}

int update_ode_start_conditions(struct ode* ODE, double* new_y, size_t n_new_y)
{
    if (ODE->n_y != n_new_y)
    {
        fprintf(stderr, "New start conditions length does not match ODE system length!\n");
        return -1;
    }

    for (int i = 0; i < ODE->n_pars; ++i)
    {
        ODE->y[i] = new_y[i];
    }

    return 0;
}

int update_ode(struct ode* ODE, double* new_pars, size_t n_new_pars, double* new_y, size_t n_new_y)
{
    int pars_check = update_ode_pars(ODE, new_pars, n_new_pars);
    int y_check = update_ode_start_conditions(ODE, new_y, n_new_y);

    if (pars_check || y_check)
    {
        return 1;
    }

    return 0;
}