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

I measured performance across 100 replicate runs at two absolute error targets:
![](./plt_benchmark_reps.png)

Stepper performance was dependent on the system, but in general was fast (apart from msadams, which was slow across the board, although very stable in performance across models - perhaps for larger systems it might compare better?)

We also measured solution similarities vs the R ```LSODA``` solver implemented in ```deSolve```:
![](./plt_sims_sep.png)
![](./plt_sims_com.png)

Solutions were similar across solvers, except for the Lorenz system which is extremely chaotic and known to be sensitive to solver choice (so not surprising).


Overall, the tests show the GSL library produces comparable results to deSolve in most cases, and with considerable performance advantages. It should be a fine inclusion into SLiM. 