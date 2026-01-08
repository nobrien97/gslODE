library(deSolve)
library(tidyverse)
library(microbenchmark)
library(cowplot)
library(plotly)

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
    X <- (t > Xstart && t <= Xstop)
    dZ <- b * X * (X^n)/(KXZ^n + X^n) * (KZ^n)/(KZ^n + Z^n) - a*Z
    return(list(dZ))
  })
}

# Setup parameters
dt <- 0.1
tmax <- 9.9

params <- split(pars, seq(nrow(pars)))
params <- lapply(params, function(x) {
  unlist(x)
})

iniState <- c(Z=0)
times <- seq(0,tmax,by=dt)

out_desolve <- data.frame(id = integer(length(params) * length(times)),
                          time = double(length(params) * length(times)),
                          X = double(length(params) * length(times)),
                          Z = double(length(params) * length(times)))

for (i in seq_along(params)) {
  solution <- ode(iniState, times, ODE_NAR, params[[i]]) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(id = i,
           X = ifelse(time > params[[i]]["Xstart"] & time <= params[[i]]["Xstop"], 1, 0)) %>%
    select(id, time, X, Z)
  range <- (((i-1) * nrow(solution)) + 1):(i * nrow(solution))
  out_desolve[range,] <- solution
}

# Plotting colours
colX   <- "#00AAFF"
pal <- c("#FF3300", "#5500FF")

plotZ <- ggplot(out_desolve %>% filter(id < 3 & id > 0)) +
  annotate("rect", xmin = params[[1]]["Xstart"], xmax = params[[1]]["Xstop"], ymin = 0, ymax = 1.05,
           alpha = .2, fill = colX) +
  geom_line(aes(time, Z, colour = factor(id)), linewidth = 1) +
  #scale_y_continuous(limits = c(0,1.05)) +
  scale_colour_manual(values = pal) +
  labs(x = "Time", y = "Z expression", colour = "Parameter combo") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")
plotZ

# Run other examples
out <- system("./build/GSLODE -i './tests/test_nar_input.csv' -f 0.1 -s rk4",
       intern = T, ignore.stderr = T)

out <- read_csv(paste0(out, collapse = "\n"),
                col_names = c("id", "time", "X", "Z"))

plotZ_gsl <- ggplot(out %>% filter(id < 3)) +
  annotate("rect", xmin = params[[1]]["Xstart"], xmax = params[[1]]["Xstop"], ymin = 0, ymax = 1.05,
           alpha = .2, fill = colX) +
  geom_line(aes(time, Z, colour = factor(id)), linewidth = 1) +
  #scale_y_continuous(limits = c(0,1.05)) +
  scale_colour_manual(values = pal) +
  labs(x = "Time", y = "Z expression", colour = "Parameter combo") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")
plotZ_gsl

# Correlation
out$id <- factor(out$id)
out$time <- factor(round(out$time, digits = 1))

out_desolve_bak <- out_desolve

out_desolve$id <- factor(out_desolve$id)
out_desolve$time <- factor(round(out_desolve$time, digits = 1))

out_diff <- full_join(out, out_desolve, by = c("id", "time"))

out_diff <- out_diff %>%
  rename(X_gsl = X.x,
         Z_gsl = Z.x,
         X_desolve = X.y,
         Z_desolve = Z.y) %>%
  group_by(id, time) %>%
  mutate(Z_diff = Z_gsl - Z_desolve)

# box plot
ggplot(out_diff,
       aes(x = time, y = Z_diff)) +
  geom_boxplot(whisker.linewidth = 0.1,
               outliers = F) +
  theme_bw() +
  labs(x = "Time",
       y = "Difference in Z expression\ndeSolve vs GSL (RK4)") +
  theme(text = element_text(size = 12))

# Doesn't appear to change over time, look at overall
ggplot(out_diff,
       aes(y = Z_diff)) +
  geom_boxplot(whisker.linewidth = 0.1,
               outliers = F) +
  theme_bw() +
  labs(y = "Difference in Z expression\ndeSolve vs GSL (RK4)") +
  theme(text = element_text(size = 12))

t.test(out_diff$Z_gsl, out_diff$Z_desolve)

# GSL is more likely to come up with solutions greater than deSolve, but they are
# on average very similar

# What if we try a different solving method
# Run other examples
out_msadams <- system("./build/GSLODE -i './tests/test_nar_input.csv' -f 0.1 -s msadams",
              intern = T, ignore.stderr = T)

out_msadams <- read_csv(paste0(out_msadams, collapse = "\n"),
                col_names = c("id", "time", "X", "Z"))

out_msadams_diff <- out_msadams
out_msadams_diff$id <- factor(out_msadams_diff$id)
out_msadams_diff$time <- factor(round(out_msadams_diff$time, digits = 1))

out_msadams_diff <- full_join(out_msadams_diff, out_desolve, by = c("id", "time"))

out_msadams_diff <- out_msadams_diff %>%
  rename(X_gsl = X.x,
         Z_gsl = Z.x,
         X_desolve = X.y,
         Z_desolve = Z.y) %>%
  group_by(id, time) %>%
  mutate(Z_diff = Z_gsl - Z_desolve)

# Boxplot
ggplot(out_msadams_diff,
       aes(x = time, y = Z_diff)) +
  geom_boxplot(whisker.linewidth = 0.1,
               outliers = F) +
  theme_bw() +
  labs(x = "Time",
       y = "Difference in Z expression\ndeSolve vs GSL (Adams)") +
  theme(text = element_text(size = 12))

# Doesn't appear to change over time, look at overall
ggplot(out_msadams_diff,
       aes(y = Z_diff)) +
  geom_boxplot(whisker.linewidth = 0.1,
               outliers = F) +
  theme_bw() +
  labs(y = "Difference in Z expression\ndeSolve vs GSL (Adams)") +
  theme(text = element_text(size = 12))

t.test(out_msadams_diff$Z_gsl, out_msadams_diff$Z_desolve)


# Benchmark

microbenchmark(system("./build/GSLODE -i './tests/test_nar_input.csv'",
                             intern = T, ignore.stderr = T)
               )
microbenchmark(
  for (i in seq_along(params)) {
  solution <- ode(iniState, times, ODE_NAR, params[[i]]) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(id = i,
           X = ifelse(time > params[[i]]["Xstart"] & time <= params[[i]]["Xstop"], 1, 0)) %>%
    select(id, time, X, Z)
  range <- (((i-1) * nrow(solution)) + 1):(i * nrow(solution))
  out_desolve[range,] <- solution
  }
)

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


# Plot VDP for interval 0-10
# Setup parameters
dt <- 0.1
tmax <- 9.9

pars <- data.frame(mu = c(1, 0))

params <- split(pars, seq(nrow(pars)))
params <- lapply(params, function(x) {
  unlist(x)
})

iniState <- c(u = 1.0, v = 0.0)
times <- seq(0,tmax,by=dt)

out_vdp_desolve <- data.frame(id = integer(length(params) * length(times)),
                          time = double(length(params) * length(times)),
                          u = double(length(params) * length(times)),
                          v = double(length(params) * length(times)))

for (i in seq_along(params)) {
  solution <- ode(iniState, times, ODE_VDP, params[[i]]) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(id = i) %>%
    select(id, time, u, v)
  range <- (((i-1) * nrow(solution)) + 1):(i * nrow(solution))
  out_vdp_desolve[range,] <- solution
}

# Plotting colours
colX   <- "#00AAFF"
pal <- c("#FF3300", "#5500FF")

plotu <-ggplot(out_vdp_desolve) +
  geom_line(aes(time, u, colour = factor(id)), linewidth = 1) +#,
            #colour = "grey") +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "u") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")

plotv <- ggplot(out_vdp_desolve) +
  geom_line(aes(time, v, colour = factor(id)), linewidth = 1)+#,
            #colour = "grey") +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "v") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")
plot_grid(plotu, plotv,
          nrow = 2, labels = "AUTO")

# Plot in GSL
out_vdp_gsl <- system("./build/GSLODE -i './tests/test_vdp_input.csv' -f 0.1 -s rk4 -o VanDerPol",
              intern = T, ignore.stderr = T)

out_vdp_gsl <- read_csv(paste0(out_vdp_gsl, collapse = "\n"),
                col_names = c("id", "time", "u", "v"))

plotu_gsl <-ggplot(out_vdp_gsl) +
  geom_line(aes(time, u, colour = factor(id)), linewidth = 1)+#,
            #colour = "grey") +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "u") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")

plotv_gsl <- ggplot(out_vdp_gsl) +
  geom_line(aes(time, v, colour = factor(id)), linewidth = 1)+#,
            #colour = "grey") +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "v") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")
plot_grid(plotu_gsl, plotv_gsl,
          nrow = 2, labels = "AUTO")

# Similarity
out_vdp_gsl$id <- factor(out_vdp_gsl$id)
out_vdp_gsl$time <- factor(round(out_vdp_gsl$time, digits = 1))
out_vdp_desolve$id <- factor(out_vdp_desolve$id)
out_vdp_desolve$time <- factor(round(out_vdp_desolve$time, digits = 1))

out_vdp_diff <- full_join(out_vdp_gsl, out_vdp_desolve, by = c("id", "time"))

out_vdp_diff <- out_vdp_diff %>%
  rename(u_gsl = u.x,
         v_gsl = v.x,
         u_desolve = u.y,
         v_desolve = v.y) %>%
  group_by(id, time) %>%
  mutate(u_diff = u_gsl - u_desolve,
         v_diff = v_gsl - v_desolve) %>%
  pivot_longer(cols = c(u_diff, v_diff), names_to = "solution", values_to = "diff",
               names_pattern = "(.*)_diff")

# Boxplot
ggplot(out_vdp_diff,
       aes(x = time, y = diff, colour = solution)) +
  geom_boxplot(whisker.linewidth = 0.1,
               outliers = T) +
  theme_bw() +
  labs(x = "Time",
       y = "Difference in u \ndeSolve vs GSL (RK4)",
       colour = "Variable") +
  theme(text = element_text(size = 12),
        legend.position = "bottom")

# Random over time
ggplot(out_vdp_diff,
       aes(y = diff, colour = solution)) +
  geom_boxplot(whisker.linewidth = 0.1,
               outliers = T) +
  theme_bw() +
  labs(y = "Difference in u \ndeSolve vs GSL (RK4)",
       colour = "Solution") +
  theme(text = element_text(size = 12),
        legend.position = "bottom")

# ZeroNo difference between them in these solutions
t.test(out_vdp_diff[out_vdp_diff$solution == "u",]$u_gsl,
       out_vdp_diff[out_vdp_diff$solution == "u",]$u_desolve)

t.test(out_vdp_diff[out_vdp_diff$solution == "v",]$v_gsl,
       out_vdp_diff[out_vdp_diff$solution == "v",]$v_desolve)


# Plot Robertson for interval t in [0, 40)
# Setup parameters
dt <- 0.1
tmax <- 40

params <- NULL

iniState <- c(X = 1, Y = 0, Z = 0)
times <- seq(0,tmax,by=dt)

out_desolve <- data.frame(id = integer(length(params) * length(times)),
                          time = double(length(params) * length(times)),
                          X = double(length(params) * length(times)),
                          Y = double(length(params) * length(times)),
                          Z = double(length(params) * length(times)))

out_rob_desolve <- ode(iniState, times, ODE_Robertson, params[[1]]) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(id = 1) %>%
    select(id, time, X, Y, Z)


plotX_rob <- ggplot(out_rob_desolve) +
  geom_line(aes(time, X), linewidth = 1) +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "X") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")

plotY_rob <- ggplot(out_rob_desolve) +
  geom_line(aes(time, Y), linewidth = 1) +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "Y") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")


plotZ_rob <- ggplot(out_rob_desolve) +
  geom_line(aes(time, Z), linewidth = 1) +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "Z") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")

plot_grid(plotX_rob, plotY_rob, plotZ_rob,
          nrow = 3, labels = "AUTO")

# GSL solution
out_rob_gsl <- system("./build/GSLODE -i './tests/test_rob_input.csv' -f 0.1 -t 40 -s msbdf -o Robertson",
                      intern = T, ignore.stderr = T)

out_rob_gsl <- read_csv(paste0(out_rob_gsl, collapse = "\n"),
                        col_names = c("id", "time", "X", "Y", "Z"))

plotX_rob_gsl <- ggplot(out_rob_gsl) +
  geom_line(aes(time, X), linewidth = 1) +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "X") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")

plotY_rob_gsl <- ggplot(out_rob_gsl) +
  geom_line(aes(time, Y), linewidth = 1) +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "Y") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")


plotZ_rob_gsl <- ggplot(out_rob_gsl) +
  geom_line(aes(time, Z), linewidth = 1) +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "Z") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")

plot_grid(plotX_rob_gsl, plotY_rob_gsl, plotZ_rob_gsl,
          nrow = 3, labels = "AUTO")




# Plot Lorenz system
dt <- 0.1
tmax <- 10

pars <- data.frame(sigma = 10,
                   rho = 28,
                   beta = 8/3)

params <- split(pars, seq(nrow(pars)))
params <- lapply(params, function(x) {
  unlist(x)
})

iniState <- c(X = 1.0, Y = 1.0, Z = 1.0)
times <- seq(0,tmax,by=dt)

out_lor_desolve <- data.frame(id = integer(length(params) * length(times)),
                              time = double(length(params) * length(times)),
                              X = double(length(params) * length(times)),
                              Y = double(length(params) * length(times)),
                              Z = double(length(params) * length(times)))

for (i in seq_along(params)) {
  solution <- ode(iniState, times, ODE_Lorenz, params[[i]]) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(id = i) %>%
    select(id, time, X, Y, Z)
  range <- (((i-1) * nrow(solution)) + 1):(i * nrow(solution))
  out_lor_desolve[range,] <- solution
}

# 3d plot
plot_ly(out_lor_desolve, x = ~X, y = ~Y, z = ~Z,
        type = "scatter3d", mode = "lines")
