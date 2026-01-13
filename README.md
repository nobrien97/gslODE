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

We measured performance across 100 replicate runs at two absolute error targets:
![](./plt_benchmark_reps.png)

We also measured the solution similarity vs the R ```LSODA``` solver implemented in ```deSolve```