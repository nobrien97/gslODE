/*
Example ODE for a negative autoregulation gene regulatory network motif.
Modelling two nodes X and Z with X activating Z and Z inhibiting itself.
X is active between timepoints t_start and t_end, other parameters are
production rate b, activation coefficient KXZ, repression coefficient KZ,
Hill coefficient n, and degradation rate a. 
*/
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

struct ode
{
    int (*ODE_fn_ptr)(double, const double*, double*, void*);
    int (*jac)(double, const double*, double*, double*, void*);
    double* pars;
    size_t n_pars;
    double* y;
    size_t n_y;
};

// ODEs
int ODE_NAR(double t, const double z[], double dzdt[], void* params);
int ODE_VanDerPol(double t, const double z[], double dzdt[], void* params);
int ODE_Lorenz(double t, const double z[], double dzdt[], void* params);
int ODE_Robertson(double t, const double z[], double dzdt[], void* params);

// Jacobian
int jac_robertson(double t, const double y[], double* dfdy, double dfdt[], void* params);
int jac_nar(double t, const double y[], double* dfdy, double dfdt[], void* params);
int jac_vanderpol(double t, const double y[], double* dfdy, double dfdt[], void* params);

// Solver and helper functions
int solve(gsl_odeiv2_driver* d, int max_time, struct ode* ODE, int par_id, double measure_interval);

const gsl_odeiv2_step_type* get_stepper_from_input(char* input_string);

int get_ODE_from_input(char* input_string, struct ode* ODE);

int free_ode_system(struct ode* ODE);

int update_ode_pars(struct ode* ODE, double* new_pars, size_t n_new_pars);

int update_ode_start_conditions(struct ode* ODE, double* new_y, size_t n_new_y);

int update_ode(struct ode* ODE, double* new_pars, size_t n_new_pars, double* new_y, size_t n_new_y);

int print_ode(struct ode* ODE, int par_id, double t);