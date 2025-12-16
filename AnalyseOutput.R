library(deSolve)
library(tidyverse)
library(microbenchmark)

# Seed
#seed <- sample(0:.Machine$integer.max, 1)
# [1] 1503521039
seed <- 1503521039
set.seed(seed)

# Generate example parameters
pars <- matrix(c(rep(c(1, 6), each = 10000),
                 abs(rnorm(50000))), nrow = 10000)

write_csv(as.data.frame(pars), "./tests/test_input.csv", col_names = F)

# Load example parameters
pars <- read_csv("./tests/test_input.csv",
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

plotZ <- ggplot(out_desolve %>% filter(id < 3)) +
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
out <- system("./build/GSLODE -i './tests/test_input.csv' -f 0.1 -s rk4",
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
out_diff <- out
out_diff$Z <- (out$Z - out_desolve$Z)
out_diff$id <- factor(out_diff$id)

ggplot(out_diff,
       aes(x = as.factor(time), y = Z)) +
  geom_boxplot(colour = "grey", outliers = F) +
  labs(x = "Time", y = "Difference in Z expression\ndeSolve vs GSL (RK4)") +
  theme_bw() +
  theme(text = element_text(size = 10))

# GSL is more likely to come up with solutions greater than deSolve than less than
# deSolve.

# What if we try a different solving method
# Run other examples
out_msadams <- system("./build/GSLODE -i './tests/test_input.csv' -f 0.1 -s msadams",
              intern = T, ignore.stderr = T)

out_msadams <- read_csv(paste0(out_msadams, collapse = "\n"),
                col_names = c("id", "time", "X", "Z"))

out_msadams_diff <- out_msadams
out_msadams_diff$Z <- (out_msadams$Z - out_desolve$Z)
out_msadams_diff$id <- factor(out_msadams_diff$id)

ggplot(out_msadams_diff,
       aes(x = as.factor(time), y = Z)) +
  geom_boxplot(colour = "grey", outliers = F) +
  labs(x = "Time", y = "Difference in Z expression\ndeSolve vs GSL (MS-Adams)") +
  theme_bw() +
  theme(text = element_text(size = 10))




microbenchmark(system("./build/GSLODE -i './tests/test_input.csv'",
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
