#include <time.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include "gsl/gsl_odeiv2.h"
#include "ode.h"

#define BUFF_SIZE 1024

char** str_split(char* a_str, const char a_delim)
{
    char** result    = 0;
    size_t count     = 0;
    char* tmp        = a_str;
    char* last_comma = 0;
    char delim[2];
    delim[0] = a_delim;
    delim[1] = 0;

    /* Count how many elements will be extracted. */
    while (*tmp)
    {
        if (a_delim == *tmp)
        {
            count++;
            last_comma = tmp;
        }
        tmp++;
    }

    /* Add space for trailing token. */
    count += last_comma < (a_str + strlen(a_str) - 1);

    /* Add space for terminating null string so caller
       knows where the list of returned strings ends. */
    count++;

    result = malloc(sizeof(char*) * count);

    if (result)
    {
        size_t idx  = 0;
        char* token = strtok(a_str, delim);

        while (token)
        {
            assert(idx < count);
            *(result + idx++) = strdup(token);
            token = strtok(0, delim);
        }
        assert(idx == count - 1);
        *(result + idx) = 0;
    }

    return result;
}

void show_help()
{
    printf(
        "---------------------NAR GSL ODE solver-----------------------\n"
        "Options:\n"
        "-h:                                    Print this help screen.\n"
        "-i:                                Provide an input .csv file.\n"
        "Example:                                          -i input.csv\n"
        "-s:                                Provide a stepper function,\n"
        "                                     one of 'rk4', 'rkf45', or\n"
        "                                                    'msadams'.\n"
        "Example:                                                -s rk4\n"
        "-f:                                     Specify an interval to\n"
        "                                        measure the system at.\n"
        "Example:                                                -f 0.1\n"
        "-t                                      Maximum solution time.\n"
        "Example:                                                -t 100\n"
        "\n"
    );
    return;
}


int main(int argc, char* argv[])
{
    fprintf(stderr, "-------- Solving ODE with GSL --------\n");

    // Time solver
    clock_t start_t, end_t;
    double total_t;
    start_t = clock();

    // Setup default parameters
    double pars[7] = {1.0, 6.0, 1.0, 1.0, 1.0, 4.0, 1.0};
    size_t n_pars = sizeof(pars) / sizeof(pars[0]);
    double t = 0.0;
    gsl_odeiv2_system sys = {NARODE, NULL, 7, &pars};

    int error_code = 0;

    // Setup output
    double y[1] = {0.0};
    // No arguments, proceed with a default ODE parameter set
    if (argc == 1) 
    {
        fprintf(stderr, "No arguments, solving default system...\n");
        gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,
                                                    1e-6, 1e-6, 0.0);

        error_code = solve(d, 100, y, 1, 1, 0);
        gsl_odeiv2_driver_free(d);

        if (error_code != 0)
        {
            fprintf(stderr, "Error in solution, return value %d\n", error_code);
            return error_code;
        }
        
        // End timer
        end_t = clock();

        fprintf(stderr, "Time taken = %.3lf\n", (double)(end_t - start_t) / CLOCKS_PER_SEC);
        return error_code;
    }

    // With arguments, read in CSV and evaluate line by line
    // TODO
    fprintf(stderr, "Loading parameter combinations...\n");
    FILE* fp;
    char* line = NULL;
    size_t len = 0;
    long read;

    char* filepath = NULL;
    const gsl_odeiv2_step_type* stepper = gsl_odeiv2_step_rkf45;
    double measure_interval = 0.1;
    double time = 10.0;
    double stepsize = 0.1;

    int opt = 0;

    while ((opt = getopt(argc, argv, ":i:s:hf:t:")) != -1)
    {
        switch (opt)
        {
        case 'i':
            filepath = optarg;
            break;
        case 's':
            stepper = get_stepper_from_input(optarg);
            break;
        case 'h':
            show_help();
            return 0;
        case 't':
            time = atof(optarg);
            break;
        case 'f':
            measure_interval = atof(optarg);
            break;
        default:
            break;
        }
    }
    
    if (stepper == NULL)
    {
        fprintf(stderr, "Invalid stepper function provided. Refer to -h for list of available steppers.\n");
        return 1;
    }

    // Read data line by line
    char buff[BUFF_SIZE];
    char** tokens;
    fp = fopen(filepath, "r");

    // ID for parameter combination
    int par_id = 1;

    gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, stepper,
                            1e-6, 1e-8, 0);

    while(fgets(buff, sizeof(buff), fp))
    {
        // Reset pars values so we can check for invalid rows with too few entries
        for (int i = 0; i < n_pars; ++i)
        {
            pars[i] = -1.0;
        }

        // Split row by comma
        tokens = str_split(buff, ',');
        if (tokens)
        {
            // Fill parameter array
            for (int i = 0; *(tokens + i); i++)
            {
                if (i > n_pars - 1)
                {
                    fprintf(stderr, "Error, too many values in row\n");
                    return 1;
                }
                pars[i] = atof(*(tokens + i));
                free(*(tokens + i));
            }

            // Check parameter array is properly filled
            for (int i = 0; i < n_pars; ++i)
            {
                if (pars[i] < 0.0)
                {
                    fprintf(stderr, "Error, either too few supplied parameters in row or negative parameters given\n");
                    return 1;

                }
            }

            // Solve this iteration
            error_code = solve(d, time, y, 1, par_id++, measure_interval);

            // Reset the driver/state
            gsl_odeiv2_driver_reset(d);
            y[0] = 0.0;
            
            if (error_code != 0)
            {
                fprintf(stderr, "Error in solution, return value %d\n", error_code);
                return error_code;
            }
        }
    }
    
    // Free the ODE driver
    gsl_odeiv2_driver_free(d);

    // End timer
    end_t = clock();
    fprintf(stderr, "Time taken = %.3lf\n", (double)(end_t - start_t) / CLOCKS_PER_SEC);

    return error_code;
}