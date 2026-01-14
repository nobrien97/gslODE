library(deSolve)
library(tidyverse)
library(microbenchmark)
library(cowplot)
library(plotly)
library(paletteer)
library(ggh4x)
library(see)
library(latex2exp)



# NAR system
ODE_NAR <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- (t > Xstart && t < Xstop)
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

# Plot VDP for interval 0-10
# Setup parameters
dt <- 0.1
tmax <- 10

pars <- data.frame(mu = 1)

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

# Over all time points
ggplot(out_vdp_diff,
       aes(y = diff, x = solution)) +
  geom_boxplot(whisker.linewidth = 0.1,
               outliers = T) +
  theme_bw() +
  labs(y = "Difference in u \ndeSolve vs GSL (RK4)",
       x = "Solution") +
  theme(text = element_text(size = 12),
        legend.position = "bottom")

# No difference between them in these solutions
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

# Similarity
out_rob_gsl$id <- factor(out_rob_gsl$id)
out_rob_gsl$time <- factor(round(out_rob_gsl$time, digits = 1))
out_rob_desolve$id <- factor(out_rob_desolve$id)
out_rob_desolve$time <- factor(round(out_rob_desolve$time, digits = 1))

out_rob_diff <- full_join(out_rob_gsl, out_rob_desolve, by = c("id", "time"))

out_rob_diff <- out_rob_diff %>%
  rename(x_gsl = X.x,
         y_gsl = Y.x,
         z_gsl = Z.x,
         x_desolve = X.y,
         y_desolve = Y.y,
         z_desolve = Z.y) %>%
  group_by(id, time) %>%
  mutate(x_diff = x_gsl - x_desolve,
         y_diff = y_gsl - y_desolve,
         z_diff = z_gsl - z_desolve) %>%
  pivot_longer(cols = c(x_diff, y_diff, z_diff), names_to = "solution", values_to = "diff",
               names_pattern = "(.*)_diff")

# Boxplot
ggplot(out_rob_diff,
       aes(x = time, y = diff, colour = solution)) +
  geom_boxplot(whisker.linewidth = 0.1,
               outliers = T) +
  theme_bw() +
  labs(x = "Time",
       y = "Difference between \ndeSolve vs GSL (msbdf)",
       colour = "Variable") +
  theme(text = element_text(size = 12),
        legend.position = "bottom")

# Over all time points
ggplot(out_rob_diff,
       aes(y = diff, colour = solution)) +
  geom_boxplot(whisker.linewidth = 0.1,
               outliers = T) +
  theme_bw() +
  labs(y = "Difference between \ndeSolve vs GSL (msbdf)",
       colour = "Solution") +
  theme(text = element_text(size = 12),
        legend.position = "bottom")

# No difference between them in these solutions
t.test(out_rob_diff[out_rob_diff$solution == "x",]$x_gsl,
       out_rob_diff[out_rob_diff$solution == "x",]$x_desolve)

t.test(out_rob_diff[out_rob_diff$solution == "y",]$y_gsl,
       out_rob_diff[out_rob_diff$solution == "y",]$y_desolve)

t.test(out_rob_diff[out_rob_diff$solution == "z",]$z_gsl,
       out_rob_diff[out_rob_diff$solution == "z",]$z_desolve)


# Plot Lorenz system
dt <- 0.01
tmax <- 100

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

plotX_lor_desolve <- ggplot(out_lor_desolve) +
  geom_line(aes(as.numeric(time), X), linewidth = 1) +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "X") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")

plotY_lor_desolve <- ggplot(out_lor_desolve) +
  geom_line(aes(as.numeric(time), Y), linewidth = 1) +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "Y") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")


plotZ_lor_desolve <- ggplot(out_lor_desolve) +
  geom_line(aes(as.numeric(time), Z), linewidth = 1) +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "Z") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")

plot_grid(plotX_lor_desolve, plotY_lor_desolve, plotZ_lor_desolve,
          nrow = 3, labels = "AUTO")


# GSL

lor_base_call <- "./build/GSLODE -i './tests/test_lorenz_input.csv' -f 0.01 -t 100 -o Lorenz"
methods <- c("rkf45", "msbdf", "bsimp", "msadams")

out_lor_gsl <- data.frame(id = numeric(length(methods) * length(times)),
                          time = numeric(length(methods) * length(times)),
                          X = numeric(length(methods) * length(times)),
                          Y = numeric(length(methods) * length(times)),
                          Z = numeric(length(methods) * length(times)))

for (i in seq_along(methods)) {
  out_this_iter <- system(paste(lor_base_call, "-s", methods[i]),
                          intern = T, ignore.stderr = T)
  out_lor_gsl[((i-1)*length(times) + 1):(length(times) * i),] <- read_csv(paste0(out_this_iter, collapse = "\n"),
                                                                          col_names = c("id", "time", "X", "Y", "Z"))
  out_lor_gsl[((i-1)*length(times) + 1):(length(times) * i),]$id <- i
}

out_lor_gsl <- out_lor_gsl %>%
  mutate(method = methods[id])


# plot_ly(out_lor_gsl, x = ~X, y = ~Y, z = ~Z,
#         type = "scatter3d", mode = "lines")


plotX_lor_gsl <- ggplot(out_lor_gsl) +
  geom_line(aes(as.numeric(time), X, colour = method), linewidth = 1) +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "X", colour = "Stepper") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")

plotY_lor_gsl <- ggplot(out_lor_gsl) +
  geom_line(aes(as.numeric(time), Y, colour = method), linewidth = 1) +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "Y", colour = "Stepper") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")


plotZ_lor_gsl <- ggplot(out_lor_gsl) +
  geom_line(aes(as.numeric(time), Z, colour = method), linewidth = 1) +
  #scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "Z", colour = "Stepper") +
  theme_bw(base_size = 16) +
  theme(text = element_text(size = 12), legend.position = "bottom")

plot_grid(plotX_lor_gsl, plotY_lor_gsl, plotZ_lor_gsl,
          nrow = 3, labels = "AUTO")

# Similarity
out_lor_gsl$time <- factor(round(out_lor_gsl$time, digits = 2))
out_lor_desolve$time <- factor(round(out_lor_desolve$time, digits = 2))

out_lor_diff <- full_join(out_lor_gsl, out_lor_desolve, by = c("time"))

out_lor_diff <- out_lor_diff %>%
  rename(x_gsl = X.x,
         y_gsl = Y.x,
         z_gsl = Z.x,
         x_desolve = X.y,
         y_desolve = Y.y,
         z_desolve = Z.y) %>%
  group_by(time) %>%
  mutate(x_diff = x_gsl - x_desolve,
         y_diff = y_gsl - y_desolve,
         z_diff = z_gsl - z_desolve) %>%
  pivot_longer(cols = c(x_diff, y_diff, z_diff), names_to = "solution", values_to = "diff",
               names_pattern = "(.*)_diff")

# Boxplot
ggplot(out_lor_diff,
       aes(x = time, y = diff, colour = solution)) +
  facet_nested("Stepper" + method~.) +
  geom_boxplot(whisker.linewidth = 0.1,
               outliers = T) +
  theme_bw() +
  labs(x = "Time",
       y = "Difference between \ndeSolve vs GSL (msbdf)",
       colour = "Variable") +
  theme(text = element_text(size = 12),
        legend.position = "bottom")

# Over all time points
ggplot(out_lor_diff,
       aes(y = diff, colour = solution)) +
  facet_nested("Stepper" + method~.) +
  geom_boxplot(whisker.linewidth = 0.1,
               outliers = T) +
  theme_bw() +
  labs(y = "Difference between \ndeSolve vs GSL",
       colour = "Solution") +
  theme(text = element_text(size = 12),
        legend.position = "bottom")

# These ones are different in x and y - and by quite a lot in some cases!
# This could be down to precision? Or solver
t.test(out_lor_diff[out_lor_diff$solution == "x",]$x_gsl,
       out_lor_diff[out_lor_diff$solution == "x",]$x_desolve)

t.test(out_lor_diff[out_lor_diff$solution == "y",]$y_gsl,
       out_lor_diff[out_lor_diff$solution == "y",]$y_desolve)

t.test(out_lor_diff[out_lor_diff$solution == "z",]$z_gsl,
       out_lor_diff[out_lor_diff$solution == "z",]$z_desolve)
