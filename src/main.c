#include <time.h>
#include <stdio.h>
#include "gsl/gsl_odeiv2.h"
#include "ode.h"




int main(int argc, char* argv[])
{
    printf("-------- Solving ODE with GSL --------\n");

    // Time solver
    clock_t start_t, end_t;
    double total_t;
    start_t = clock();

    // Setup default parameters
    double pars[7] = {1.0, 6.0, 1.0, 1.0, 1.0, 4.0, 1.0};

    gsl_odeiv2_system sys = {NARODE, NULL, 7, &pars};

    gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,
                                                1e-6, 1e-6, 0.0);

    // Setup output
    double y[1] = {0.0};
    // No arguments, proceed with a default ODE parameter set
    if (argc == 1) 
    {
        printf("No arguments, solving default system...\n");

        int success = solve(d, 100, y, 1);
        gsl_odeiv2_driver_free(d);

        if (success != 0)
        {
            fprintf(stderr, "Error in solution, return value %d", success);
            return success;
        }
        
        // End timer
        end_t = clock();

        printf("Time taken = %.3lf\n", (double)(end_t - start_t) / CLOCKS_PER_SEC);
        return 0;
    }

    // With arguments, read in CSV and evaluate line by line
    // TODO
    printf("Loading parameter combinations...\n");
    return 0;
}