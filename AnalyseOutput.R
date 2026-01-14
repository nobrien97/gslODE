library(deSolve)
library(tidyverse)
library(microbenchmark)
library(cowplot)
library(plotly)
library(paletteer)
library(ggh4x)
library(see)
library(latex2exp)

# Seed
#seed <- sample(0:.Machine$integer.max, 1)
# [1] 1503521039
seed <- 1503521039
set.seed(seed)

# Generate example parameters
pars <- matrix(c(rep(c(1, 6), each = 10000),
                 abs(rnorm(50000))), nrow = 10000)

write_csv(as.data.frame(pars), "./tests/test_nar_input.csv", col_names = F)

# Load example parameters
pars <- read_csv("./tests/test_nar_input.csv",
                 col_names = c("Xstart", "Xstop", "b", "KXZ", "KZ", "n", "a"))

# Run R solution for comparison
## R code modified from Jan Engelstaedter

# NAR system
ODE_NAR <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- (t >= Xstart && t <= Xstop)
    dZ <- b * X * (X^n)/(KXZ^n + X^n) * (KZ^n)/(KZ^n + Z^n) - a*Z
    return(list(dZ))
  })
}

# Other systems
## Van der Pol
ODE_VDP <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    du <- v
    dv <- mu * (1 - u * u) * v - u
    return(list(c(du, dv)))
  })
}

## Lorenz
ODE_Lorenz <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    dX <- sigma*(Y - X)
    dY <- X*(rho - Z) - Y
    dZ <- X*Y - beta*Z
    return(list(c(dX, dY, dZ)))
  })
}

## Robertson
ODE_Robertson <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    dX <- -0.04 * X + 1e4 * Y * Z
    dY <- 0.04 * X - 1e4 * Y * Z - 3e7 * Y * Y
    dZ <- 3e7 * Y * Y
    return(list(c(dX, dY, dZ)))
  })
}


# Benchmark
a_err <- 10^(-16:-6)
r_err <- 10^(-16:-6)
methods <- c("rkf45", "msbdf", "bsimp", "msadams")
models <- c("NAR", "Lorenz", "Robertson", "VanDerPol")
model_inputs <- c("./tests/test_nar_input.csv", "./tests/test_lorenz_input.csv",
                  "./tests/test_rob_input.csv", "./tests/test_vdp_input.csv")


err <- expand.grid(a_err, r_err, methods, models)

err <- err %>%
  mutate(method_index = as.integer(as.factor(Var4)),
         filename = model_inputs[method_index])

base_call <- "./build/GSLODE -f 0.01 -t 100 -b"

out_bench <- data.frame(i = numeric(nrow(err)),
                            a_err = numeric(nrow(err)),
                            r_err = numeric(nrow(err)),
                            t = numeric(nrow(err)))

pb <- progress::progress_bar$new(
  format = "[:bar] :current/:total (:percent eta: :eta)", total = nrow(err))
pb$tick(0)

for (i in seq_len(nrow(err))) {
  str_call <- paste(base_call, "-i", err[i,]$filename, "-o", err[i,]$Var4,
                    "-a", err[i,]$Var1, "-r", err[i,]$Var2, "-s", err[i,]$Var3)

  out <- system(str_call, intern = T, ignore.stderr = T)
  out_bench[i,] <- read_csv(paste0(out, "\n"), col_names = F, show_col_types = F)
  out_bench[i,1] <- i
  pb$tick(1)
}

out_bench <- out_bench %>%
  mutate(method = err[i,]$Var3,
         model = err[i,]$Var4)


ggplot(out_bench %>%
         rename(err_a = a_err,
                err_r = r_err) %>%
         pivot_longer(cols = c(err_a, err_r), names_prefix = "err_",
                      names_to = "err_type", values_to = "err_value"),
       aes(x = as.factor(err_value), y = t, colour = err_type)) +
  facet_nested("ODE System" + model ~ "Stepper" + method) +
  geom_boxplot() +
  scale_y_log10() +
  scale_x_discrete(labels = paste(16:6)) +
  labs(x = TeX("Error $(10^{-x})$"), y = "Time (s)", colour = "Error type") +
  scale_colour_manual(values = c("#09283CFF", "#CB1724FF"), labels = c("Absolute", "Relative")) +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "bottom") -> plt_bench_all

plt_bench_all
ggsave("plt_benchmark_err.png", device = png, bg = "white", width = 9*1.5, height = 5*1.5)





# For a high and a low amount of absolute error, run many replicates to compare repeatability
a_err <- c(10^-16, 10^-6)
err <- expand.grid(a_err, methods, models)

err <- err %>%
  mutate(method_index = as.integer(as.factor(Var3)),
         filename = model_inputs[method_index])

# Repeat each call for replicates
err <- err %>%
  slice(rep(1:n(), each = 100))

base_call <- "./build/GSLODE -f 0.01 -t 100 -b"

out_bench_reps <- data.frame(i = numeric(nrow(err)),
                        a_err = numeric(nrow(err)),
                        r_err = numeric(nrow(err)),
                        t = numeric(nrow(err)))

pb <- progress::progress_bar$new(
  format = "[:bar] :current/:total (:percent eta: :eta)", total = nrow(err))
pb$tick(0)

for (i in seq_len(nrow(err))) {
  str_call <- paste(base_call, "-i", err[i,]$filename, "-o", err[i,]$Var3,
                    "-a", err[i,]$Var1, "-s", err[i,]$Var2)

  out <- system(str_call, intern = T, ignore.stderr = T)
  out_bench_reps[i,] <- read_csv(paste0(out, "\n"), col_names = F, show_col_types = F)
  out_bench_reps[i,1] <- i
  pb$tick(1)
}

out_bench_reps <- out_bench_reps %>%
  mutate(method = err[i,]$Var2,
         model = err[i,]$Var3)

ggplot(out_bench_reps %>%
         rename(err_a = a_err,
                err_r = r_err) %>%
         pivot_longer(cols = c(err_a, err_r), names_prefix = "err_",
                      names_to = "err_type", values_to = "err_value"),
       aes(x = method, y = t, colour = model)) +
  facet_nested("Absolute error" + as.factor(err_value)~.) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x = "Stepping function", y = "Time (s)", colour = "Model") +
  scale_colour_paletteer_d("nationalparkcolors::Everglades") +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "bottom") -> plt_bench_reps
plt_bench_reps
ggsave("plt_benchmark_reps.png", device = png, bg = "white", width = 6, height = 8)

# Similarity
dt <- 0.01
tmax <- 100
times <- seq(0,tmax,by=dt)
err <- expand.grid(methods, models)

err <- err %>%
  mutate(model_index = as.integer(as.factor(Var2)),
         filename = model_inputs[model_index])

base_call <- "./build/GSLODE -f 0.01 -t 100"

out_comp_gsl <- data.frame(id = numeric(nrow(err) * length(times)),
                          time = numeric(nrow(err) * length(times)),
                          X = numeric(nrow(err) * length(times)),
                          Y = numeric(nrow(err) * length(times)),
                          Z = numeric(nrow(err) * length(times)))

pb <- progress::progress_bar$new(
  format = "[:bar] :current/:total (:percent eta: :eta)", total = nrow(err))
pb$tick(0)

for (i in seq_len(nrow(err))) {
  out_this_iter <- system(paste(base_call, "-s", err[i,]$Var1, "-o",
                                err[i,]$Var2, "-i", err[i,]$filename),
                          intern = T, ignore.stderr = T)
  out_comp_gsl[((i-1)*length(times) + 1):(length(times) * i),] <- read_csv(paste0(out_this_iter, collapse = "\n"),
                                                                          col_names = c("id", "time", "X", "Y", "Z"),
                                                                          show_col_types = F)
  out_comp_gsl[((i-1)*length(times) + 1):(length(times) * i),]$id <- i
  pb$tick(1)
}

out_comp_gsl <- out_comp_gsl %>%
  mutate(method = err[id,]$Var1,
         model = err[id,]$Var2)

# Run models in deSolve

# Run R solution for comparison

# Set parameter sets
params <- list(NAR = c(Xstart = 1, Xstop = 80, b = 1, KXZ = 1, KZ = 1, n = 1, a = 1),
               Lorenz = c(sigma = 10, rho = 28, beta = 2.666667),
               Robertson = NULL,
               VanDerPol = c(mu = 1))

times <- seq(0,tmax,by=dt)

getODEFromModelName <- function(x) {
  switch (x,
    "NAR" = { return(ODE_NAR) },
    "Lorenz" = { return(ODE_Lorenz) },
    "Robertson" = { return(ODE_Robertson) },
    "VanDerPol" = { return(ODE_VDP) }
  )
}

getIniStateFromModelName <- function(x) {
  switch (x,
    "NAR" = { return(c(Z = 0)) },
    "Lorenz" = { return(c(X = 1, Y = 1, Z = 1)) },
    "Robertson" = { return(c(X = 1, Y = 0, Z = 0)) },
    "VanDerPol" = { return(c(u = 1.0, v = 0.0)) }
  )
}

out_comp_desolve <- data.frame(id = integer(length(models) * length(times)),
                          time = double(length(models) * length(times)),
                          X = double(length(models) * length(times)),
                          Y = double(length(models) * length(times)),
                          Z = double(length(models) * length(times)))

pb <- progress::progress_bar$new(
  format = "[:bar] :current/:total (:percent eta: :eta)", total = nrow(err))
pb$tick(0)
for (i in seq_along(models)) {
  ODE <- getODEFromModelName(models[i])
  iniState <- getIniStateFromModelName(models[i])
  solution <- ode(iniState, times, ODE, params[[i]]) %>%
    as.data.frame() %>%
    as_tibble()
  range <- (((i-1) * nrow(solution)) + 1):(i * nrow(solution))
  out_comp_desolve[range,]$time <- solution$time
  out_comp_desolve[range,3:5] <- c(solution[,-1], vector("double", 3 - ncol(solution[,-1])))
  if (models[i] == "NAR") {
    out_comp_desolve[range,4] <- out_comp_desolve[range,3]
    out_comp_desolve[range,3] <- out_comp_desolve[range,2] >= 1 & out_comp_desolve[range,2] < 80
  }
  out_comp_desolve[range,]$id <- i
  pb$tick(1)
}

# Combine datasets
out_comp_desolve <- out_comp_desolve %>%
  mutate(model = models[id])

out_comp_desolve <- out_comp_desolve %>%
  select(-id) %>%
  mutate(method = "deSolve (LSODA)")

out_comp_gsl <- out_comp_gsl %>%
  select(-id)

out_comp_gsl$time <- factor(out_comp_gsl$time)

# Repeat each group
out_diff <- out_comp_desolve %>%
  group_by(model) %>%
  slice(rep(1:n(), length(methods)))

out_diff$method <- out_comp_gsl$method
out_diff$time <- factor(out_diff$time)
out_diff <- full_join(out_diff, out_comp_gsl, by = c("time", "model", "method"))

out_diff_sep <- out_diff %>%
  rename(X_desolve = X.x,
         Y_desolve = Y.x,
         Z_desolve = Z.x,
         X_gsl = X.y,
         Y_gsl = Y.y,
         Z_gsl = Z.y) %>%
  group_by(time, model, method) %>%
  mutate(X_diff = X_gsl - X_desolve,
         Y_diff = Y_gsl - Y_desolve,
         Z_diff = Z_gsl - Z_desolve)

out_diff_sep <- out_diff_sep %>%
  pivot_longer(c(X_diff, Y_diff, Z_diff), names_to = "output", values_to = "diff",
               names_pattern = "(.*)_diff")

# Remove rows which don't make sense from NAR and VDP
out_diff_sep <- out_diff_sep %>%
  filter(!(model == "NAR" & output == "Z"),
         !(model == "NAR" & output == "X"),
         !(model == "VanDerPol" & output == "Z"))

# Rename outputs
out_diff_sep[out_diff_sep$model == "NAR" & out_diff_sep$output == "Y",]$output <- "Z"
out_diff_sep[out_diff_sep$model == "VanDerPol" & out_diff_sep$output == "X",]$output <- "U"
out_diff_sep[out_diff_sep$model == "VanDerPol" & out_diff_sep$output == "Y",]$output <- "V"

# Boxplot
ggplot(out_diff_sep %>%
         mutate(output = factor(output, levels = c("X", "Y", "Z", "U", "V"))),
       aes(x = method, y = diff, colour = output)) +
  facet_nested("Model" + model ~ .,
               scales = "free_y") +
  geom_boxplot(whisker.linewidth = 0.1,
               outliers = F) +
  theme_bw() +
  scale_colour_paletteer_d("calecopal::superbloom2") +
  guides(colour = guide_legend(
    nrow = 2, reverse = F, byrow = T
  )) +
  labs(x = "Stepper",
       y = "Difference between deSolve vs GSL",
       colour = "Output") +
  theme(text = element_text(size = 12),
        legend.position = "bottom") -> plt_sim_sep
plt_sim_sep
ggsave("plt_sims_sep.png", device = png, bg = "white", width = 6, height = 7)


# Total difference across all traits
out_diff_com <- out_diff %>%
  rename(X_desolve = X.x,
         Y_desolve = Y.x,
         Z_desolve = Z.x,
         X_gsl = X.y,
         Y_gsl = Y.y,
         Z_gsl = Z.y) %>%
  group_by(time, model, method) %>%
  # Since some models have different numbers of outputs,
  # reset invalid comparisons so they don't contribute
  mutate(Z_gsl = if_else(model == "NAR" | model == "VanDerPol", 0, Z_gsl),
    diff = abs(X_gsl - X_desolve) + abs(Y_gsl - Y_desolve) + abs(Z_gsl - Z_desolve))

# Boxplot
ggplot(out_diff_com,
       aes(x = method, y = diff)) +
  facet_nested("Model" + model ~ .,
               scales = "free_y") +
  geom_boxplot(whisker.linewidth = 0.1,
               outliers = F) +
  theme_bw() +
  labs(x = "Stepper",
       y = "Total difference between deSolve vs GSL") +
  theme(text = element_text(size = 12),
        legend.position = "bottom") -> plt_sim_com
plt_sim_com
ggsave("plt_sims_com.png", device = png, bg = "white", width = 6, height = 7)

