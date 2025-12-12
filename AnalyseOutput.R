library(deSolve)
library(tidyverse)
library(RColorBrewer)


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
dt <- 1
tmax <- 100

params <- split(pars, seq(nrow(pars)))
params <- lapply(params, function(x) {
  unlist(x)
})

iniState <- c(Z=0)
times <- seq(0,tmax,by=dt)
solution <- ode(iniState, times, ODE_NAR, params[[1]]) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(X = ifelse(time > params[[1]]["Xstart"] & time <= params[[1]]["Xstop"], 1, 0)) %>%
  select(time, X, Z)

# Plotting colours
colX   <- brewer.pal(6, "Paired")[2]
colZ   <- brewer.pal(6, "Paired")[4]

plotZ <- ggplot(solution) +
  annotate("rect", xmin = params[[1]]["Xstart"], xmax = params[[1]]["Xstop"], ymin = 0, ymax = 1.05,
           alpha = .2, fill = colX) +
  geom_line(aes(time, Z), color = colZ, linewidth = 1) +
  scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "Z expression") +
  theme_bw(base_size = 16)
plotZ

# Run other examples
out <- system("./build/GSLODE -i './tests/test_input.csv'",
       intern = T, ignore.stderr = T)

out <- read_csv(paste0(out, collapse = "\n"),
                col_names = c("id", "time", "X", "Z"))

plotZ_gsl <- ggplot(out) +
  annotate("rect", xmin = params[[1]]["Xstart"], xmax = params[[1]]["Xstop"], ymin = 0, ymax = 1.05,
           alpha = .2, fill = colX) +
  geom_line(aes(time, Z, colour = as.factor(id)), linewidth = 1) +
  scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "Z expression", colour = "Parameter combo") +
  theme_bw(base_size = 16)
plotZ_gsl

