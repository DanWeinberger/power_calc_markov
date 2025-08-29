library(lubridate)
library(msm)
library(tidyverse)

# Constants
set.seed(123)

n_individuals <- 100
n_weeks <- 52
n_days <- n_weeks * 7
p_infect_community_base <- 0.001  # Baseline example community infection rate

# Seasonal variation in community infection rates (sinusoidal)
community_infection_rate <- function(day) {
  p_infect_community_base * (1 + 0.5 * sin(2 * pi * (day / 365)))
}

# Mean duration of being infected
mean_duration_infected <- 14

simulate_transitions <- function(states, community_infection_rate, mean_duration_infected, n_days) {
  for (day in 2:n_days) {
    community_rate <- community_infection_rate(day)
    for (i in 1:nrow(states)) {
      current_state <- states[i, day - 1]
      if (current_state == 1) {  # S
        new_state <- ifelse(rbinom(1, 1, community_rate) == 1, 2, 1)  # Move to I or stay in S
      } else if (current_state == 2) {  # I
        new_state <- ifelse(rbinom(1, 1, 1 / mean_duration_infected) == 1, 3, 2)  # Move to R or stay in I
      } else if (current_state == 3) {  # R
        new_state <- 3  # Stay in R
      }
      states[i, day] <- new_state
    }
  }
  return(states)
}

# Initialize states: 1 = S, 2 = I, 3 = R
states <- matrix(1, nrow = n_individuals, ncol = n_days)
states <- simulate_transitions(states, community_infection_rate, mean_duration_infected, n_days)

# Convert to weekly observations
sample_days <- seq(1, n_days, by = 7)
weekly_states <- states[, seq(1, n_days, by = 7)]

# Convert data to long format for model fitting
simulation_data <- data.frame(
  individual = rep(1:n_individuals, times = n_weeks),
  week = rep(sample_days, each = n_individuals),
  state = as.vector(weekly_states)
) %>%
  arrange(individual, week) %>%
  group_by(individual) %>%
  filter(!(state == lag(state, default = first(state)) & state==3  )) %>% #censor if R
  ungroup()

head(simulation_data)

#head(simulation_data)


####MARKOV MODEL 


Q <- rbind(
  c(0, 1, 0),  # From S: S->I
  c(0, 0, 1 / mean_duration_infected),  # From I: I->R
  c(0, 0, 0)  # From R: no transitions
)

# Fit Markov transition model
model <- msm(state ~ week, subject = individual, data = simulation_data,
             qmatrix = Q, gen.inits = TRUE)

# Summary of the model
summary(model)

# Extract transition probabilities
pmatrix.msm(model)

