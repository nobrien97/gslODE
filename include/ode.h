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

int NARODE (double t, const double z[], double dzdt[], void* params);

int solve(gsl_odeiv2_driver* d, int max_time, double y[], size_t nTraits, int par_id, double measure_interval);

const gsl_odeiv2_step_type* get_stepper_from_input(char* input_string);