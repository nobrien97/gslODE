# gslODE
Tests for GSL ODE solver implementation

This repo implements four ODE systems
- Negative autoregulation (NAR), a common gene regulatory network motif
- The Lorenz system, a highly chaotic weather model
- The Robertson problem, an extremely stiff ODE system modelling chemical kinetics
- The Van Der Pol oscillator, an oscillating ODE

and tests four different stepper functions
- Explicit embedded Runge-Kutta-Fehlberg (4, 5) method (rkf45)
- Implicit Bulirsch-Stoer method of Bader and Deuflhard (bsimp)
- Variable-coefficient linear multistep Adams method in Nordsieck form (msadams)
- Variable-coefficient linear multistep backward differentiation formula (BDF) method in Nordsieck form (msbdf)

I ran performance tests on a AMD Ryzen 7 5800HS @3.2GHz with the test app compiled with -O3 via gcc. Each model was run from time 0 to 100, with measurements taken every 0.01 units.

I measured performance across 100 replicate runs at two absolute error targets:
![](./plt_benchmark_reps.png)

And across many error targets:
![](./plt_benchmark_err.png)

Stepper performance was dependent on the system, but in general was fast (apart from msadams, which was slow across the board, although fairly stable in performance across all models except the NAR - perhaps for larger systems it might compare better?)

We also measured solution similarities vs the R ```LSODA``` solver implemented in ```deSolve```:
![](./plt_sims_sep.png)

With each output being one of the solutions for the ODE.
The total difference across all equations is given below:

![](./plt_sims_com.png)


Solutions were similar across solvers, except for the Lorenz system which is extremely chaotic and known to be sensitive to solver choice (so not surprising).

Overall, the tests show the GSL library produces comparable results to deSolve in most cases, and with considerable performance advantages. It should be a fine inclusion into SLiM. 


## SLiM Implementation

The Eidos API is proposed as

```
object<DataFrame>$ solveODESystem(s system, numeric initialValues, numeric parValues, numeric maxTime, numeric measureInterval, [s stepper = "rkf45"], [f epsRel = 1e-10], [f epsAbs = 1e-10])
```

<p class="p3"><span class="s5">Solves a system of <b>first-order ordinary differential equations</b> from time 0 to maxTime, specified by the argument system with starting conditions initialValues and parameters parValues. Outputs are measured in intervals of measureInterval.</span> Returns a dataframe with columns id for the solution index, t, for the solution time and a column for each variable. id is meaningful only when initialValues and parValues are matrices (see below).

Available systems are "polynomial" for a polynomial equation of form <i>dx/dt</i> = a<sub>n</sub>x<sup>n</sup> + a<sub>n-1</sub>x<sup>n-1</sup> + ... a<sub>2</sub>x<sup>2</sup> + a<sub>1</sub>x + a<sub>0</sub>, "simplereg" for a simple gene regulation model <i>dX/dt</i> = β - αX, "+Hillreg" for a positive Hill function gene regulation model <i>dX/dt</i> = β X<sup>h</sup>/(K<sup>h</sup>+X<sup>h</sup>)-αX, "-Hillreg" for a negative Hill function gene regulation model <i>dX/dt</i> = β K<sup>h</sup>/(K<sup>h</sup>+X<sup>h</sup>)-αX, or a path to a dynamic library (.dll or .so file) containing the function for a user-defined system.</p>
<p>Available steppers are "rkf45" (the default), for a Runge-Kutta-Fehlberg (4, 5) stepper, "bsimp" for an Implicit Bulirsch-Stoer method, "msadams" for a variable-coefficient linear multistep Adams method in Nordsieck form, or "msbdf" for a variable-coefficient linear multistep backward differentiation method in Nordsieck form.</p>
<p>epsRel and epsAbs determine the relative and absolute acceptable error for each step, respectively. If the solver cannot attain this error, it will reduce the step size and try again.</p>
<p>initialValues and parValues may be supplied as vectors or as n*m and n*k matrices, where n is the number of different initial value/parameter values to evaluate, m is the number of solutions, and k is the number of parameters. When initialValues and parValues are matrices, the output id column refers to the row number in these matrices. When system, maxTime, measureInterval, stepper, epsRel, and/or epsAbs are length n, each value will be applied according to each solution. When these values are singleton, they are applied across all solutions.</p>
<p>Polynomial, simple regulation, and Hill Function systems have one solution output; when using these systems, initialValues should either be singleton, or length n. Polynomial systems have p parameters, where p is the order of the polynomial plus one (e.g. a quadratic function would have p=2+1=3). Hence, parValues should either be length p or a n * p matrix. Parameters should be input in decreasing order (a<sub>n</sub>, a<sub>n-1</sub>, ... a<sub>0</sub>). Simple regulation has two parameter values, α and β, and should be length 2 or a n * 2 matrix. Hill function systems have 4 parameters, α, β, K, and h, and should be length 4 or a n * 4 matrix. parValues should be ordered according to the order given here.</p>
<p>User-defined systems are constructed by pre-compiled libraries containing the ODE function with signature int func(double t, const double y[], double dydt[], void* params), integers n_pars_ode and n_vars_ode for the number of parameters and variables, and optionally a function for the Jacobian matrix, with signature int jac(double t, const double y[], double* dfdy, double dfdt[], void* params). The Jacobian is required to use bsimp and msbdf step functions. The library should be constructed with C-style linkage and link against gsl_errno.h, and gsl_matrix.h (gsl_matrix.h is only required if using the Jacobian). An example is given below.</p>

#include <gsl/gsl_errno.h><br>
#include <gsl/gsl_matrix.h><br>
extern "C"<br>
<br>
// Van Der Pol Oscillator<br>
<br>
int n_pars_ode = 1;<br>
int n_vars_ode = 2;<br>
int func(double t, const double y[], double f[], void\* params)<br>
{<br>
&nbsp;&nbsp;&nbsp;&nbsp;    double mu = \*(double \*)params;<br>
&nbsp;&nbsp;&nbsp;&nbsp;    f[0] = y[1];<br>
&nbsp;&nbsp;&nbsp;&nbsp;    f[1] = -y[0] - mu\*y[1]\*(y[0]\*y[0] - 1);<br>
&nbsp;&nbsp;&nbsp;&nbsp;    return GSL_SUCCESS;<br>
}<br>
<br>
// Jacobian<br>
int jac(double t, const double y[], double\* dfdy, double dfdt[], void\* params)<br>
{<br>
&nbsp;&nbsp;&nbsp;&nbsp;    double mu = \*(double \*)params;<br>
&nbsp;&nbsp;&nbsp;&nbsp;    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);<br>
&nbsp;&nbsp;&nbsp;&nbsp;    gsl_matrix \* m = &dfdy_mat.matrix;<br>
&nbsp;&nbsp;&nbsp;&nbsp;    gsl_matrix_set (m, 0, 0, 0.0);<br>
&nbsp;&nbsp;&nbsp;&nbsp;    gsl_matrix_set (m, 0, 1, 1.0);<br>
&nbsp;&nbsp;&nbsp;&nbsp;    gsl_matrix_set (m, 1, 0, -2.0\*mu\*y[0]\*y[1] - 1.0);<br>
&nbsp;&nbsp;&nbsp;&nbsp;    gsl_matrix_set (m, 1, 1, -mu\*(y[0]\*y[0] - 1.0));<br>
&nbsp;&nbsp;&nbsp;&nbsp;    dfdt[0] = 0.0;<br>
&nbsp;&nbsp;&nbsp;&nbsp;    dfdt[1] = 0.0;<br>
&nbsp;&nbsp;&nbsp;&nbsp;    return GSL_SUCCESS;<br>
}<br>




-- End docstring --
