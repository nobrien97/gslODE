#include <ode.h>
#include <math.h>
#include <string.h>


int ODE_Poly(double t, const double y[], double dydt[], void* params)
{
    // Polynomial function
    double* p = (double*)params;

    // First "parameter" is the length of the array (inclusive of this)
    size_t n = (int)p[0];
    double x = y[0];    // GSL syntax
    double sum = 0.0;
    for (int i = 1; i < n; ++i) 
    {
        sum += p[i] * pow(x, n-1-i);
    }
    dydt[0] = sum;
    return GSL_SUCCESS;

}

int ODE_NAR(double t, const double y[], double dydt[], void* params)
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
    double Zn = pow(y[0], n);

    dydt[0] = b * (Xn / (KXZn + Xn)) * (KZn / (KZn + Zn)) - a * y[0];
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
    return GSL_SUCCESS;
}

int ODE_Robertson(double t, const double y[], double f[], void* params)
{
    // Parameters
    (void)(t);
    double* p = (double*)params;
    double c1 = p[0];
    double c2 = p[1];
    double c3 = p[2];

    double t1 = c1 * y[0];
    double t2 = c2 * y[1] * y[2];
    double t3 = c3 * y[1] * y[1];

    // Equations
    f[0] = -t1 + t2;
    f[1] = t1 - t3 - t2;
    f[2] = t3;
    return GSL_SUCCESS;
}

int jac_robertson(double t, const double y[], double* dfdy, double dfdt[], void* params)
{
    (void)(t); /* avoid unused parameter warning */
    double* p = (double*)params;
    double c1 = p[0];
    double c2 = p[1];
    double c3 = p[2];

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 3, 3);
    gsl_matrix *m = &dfdy_mat.matrix;

    double xdot_y = c2 * y[2];
    double xdot_z = c2 * y[1];
    double zdot_y = 2 * c3 * y[1];
   
    gsl_matrix_set (m, 0, 0, -c1);
    gsl_matrix_set (m, 0, 1, xdot_y);
    gsl_matrix_set (m, 0, 2, xdot_z);

    gsl_matrix_set (m, 1, 0, c1);
    gsl_matrix_set (m, 1, 1, -(xdot_y + zdot_y));
    gsl_matrix_set (m, 1, 2, -xdot_z);

    gsl_matrix_set (m, 2, 0, 0.0);
    gsl_matrix_set (m, 2, 1, zdot_y);
    gsl_matrix_set (m, 2, 2, 0.0);

    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    dfdt[2] = 0.0;
    return GSL_SUCCESS;
}

int jac_nar(double t, const double y[], double* dfdy, double dfdt[], void* params)
{
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 1, 1);
    gsl_matrix *m = &dfdy_mat.matrix;

    // Extract parameters for clarity
    double* p = (double*)params;
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
    double Zn = pow(y[0], n);
       
    gsl_matrix_set (m, 0, 0, -b * (Xn / (KXZn + Xn)) * (KZn * n * pow(y[0], n-1)) / ((KZn + Zn) * (KZn + Zn)) - a);

    dfdt[0] = 0.0;

    return GSL_SUCCESS;
}

int jac_vanderpol(double t, const double y[], double* dfdy, double dfdt[], void* params)
{
    (void)(t); /* avoid unused parameter warning */
    double mu = *(double *)params;
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, 1.0);
    gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
    gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
}

int jac_lorenz(double t, const double y[], double* dfdy, double dfdt[], void* params)
{
    (void)(t); /* avoid unused parameter warning */
    double* p = (double*)params;
    double sigma = p[0];
    double rho = p[1];
    double beta = p[2];

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 3, 3);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, -sigma);
    gsl_matrix_set (m, 0, 1, sigma);
    gsl_matrix_set (m, 0, 2, 0.0);
    
    gsl_matrix_set (m, 1, 1, rho - y[2]);
    gsl_matrix_set (m, 1, 1, -1.0);
    gsl_matrix_set (m, 1, 2, -y[0]);

    gsl_matrix_set (m, 2, 0, y[1]);
    gsl_matrix_set (m, 2, 1, y[0]);
    gsl_matrix_set (m, 2, 2, -beta);

    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    dfdt[2] = 0.0;
    return GSL_SUCCESS;
}



/// @brief Solves an ODE.
/// @param d The ODE driver.
/// @param max_time Maximum time to simulate.
/// @param y Array of initial states.
/// @param n_y Length of y.
/// @param par_id Index of this solution in the matrix of input parameters.
/// @param measure_interval How often to measure during the simulation.
/// @return GSL_SUCCESS if successful, a GSL error code if not.
int solve(gsl_odeiv2_driver* d, int max_time, struct ode* ODE, int par_id, double measure_interval, int benchmark)
{
    double t = 0;
    int status = 0;
    double* y = ODE->y;

    int max_iter = max_time / measure_interval;
    for (int i = 0; i <= max_iter; ++i)
    {
        double ti = i * measure_interval;
        status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS)
        {
            fprintf(stderr, "Stepper function %s failed on step %d - error %s\n", d->s->type->name, i, gsl_strerror(status));
            break;
        }
        
        // print solution to std out
        if (benchmark == 0)
        {
            print_ode(ODE, par_id, t);
        }
    }
    return GSL_SUCCESS;
}

int print_ode(struct ode* ODE, int par_id, double t)
{
    size_t n_y = ODE->n_y;
    double* y = ODE->y;
    
    // Handle ODE-specific printing
    if (ODE->ODE_fn_ptr == ODE_NAR)
    {
        // Print X as well
        double tstart = ODE->pars[0];
        double tend = ODE->pars[1];
        int X = t >= tstart && t < tend;
        printf("%d, %.5le, %d, ", par_id, t, X);
    }
    else
    {
        printf("%d, %.5le, ", par_id, t);
    }


    // Print results
    for (size_t j = 0; j < n_y-1; ++j)
    {
        printf("%.5le,", y[j]);
    }
    // Last result with newline
    printf("%.5le\n", y[n_y-1]);
    
    return 0;
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

    if (strcmp("msbdf", input_string) == 0)
    {
        return gsl_odeiv2_step_msbdf;
    }

    if (strcmp("bsimp", input_string) == 0)
    {
        return gsl_odeiv2_step_bsimp;
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
        ODE->jac = &jac_nar;
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
        ODE->jac = &jac_vanderpol;
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
        ODE->jac = &jac_lorenz;
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
        ODE->jac = &jac_robertson;
        ODE->n_pars = 3;
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

int method_requires_jacobian(const char* method_name)
{
    if (method_name == NULL)
    {
        fprintf(stderr, "Invalid method name - see -h for available models.\n");
        return 1;
    }

    if (strcmp(method_name, "msbdf") == 0)
    {
        return 1;
    }

    if (strcmp(method_name, "bsimp") == 0)
    {
        return 1;
    }

    return 0;

}