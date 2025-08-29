##This code generates data using a 2 state discrete-time Markov model
#We treat each year independently--this is for simplicity, otherwise need to move from R->S in season 2

library(lubridate)
library(msm)
library(tidyverse)
library(pbapply)
# Constants
set.seed(123)

nsim=500
n_individuals <- 240 #40 people per ward, on 3 wards (must be a multiple of 3!), in 2 years (each year considered independent)
n_weeks <- 52 #2 year (but just deal with 1 year of time because )
n_days <- n_weeks * 7
p_infect_community_base <- 0.00075  # Baseline community infection rate--yields ~23% of people getting infected using the seasonal parameters below
ward_epidemic_multiplier <- 2
ward_epidemics <- list(c(100,145),
                       c(500,545)) #when are there epidemics on the wward?

#For influenza, assume 23% of HCW infected during season: https://pmc.ncbi.nlm.nih.gov/articles/PMC2352563/

# Seasonal variation in community infection rates (sinusoidal)
community_infection_rate <- function(day, id) {
  if(
     (day>ward_epidemics[[1]][1] & day<ward_epidemics[[1]][2])|(day>ward_epidemics[[2]][1] & day<ward_epidemics[[2]][2])
    ){
    ward_epidemic_multiplier * p_infect_community_base * (1 + 0.5 * sin(2 * pi * (day / 365)))
  }else{
    p_infect_community_base * (1 + 0.5 * sin(2 * pi * (day / 365)))
    
  }
}

# Mean duration of being infected
mean_duration_infected <- 14

ward_individual_mapping <- cbind.data.frame(
  individual=1:n_individuals,
  ward = rep(1:3, each=n_individuals/3)
)

simulate_transitions <- function(states, community_infection_rate, mean_duration_infected, n_days) {
  for (day in 2:n_days) {
    for (i in 1:n_individuals) {
      community_rate <- community_infection_rate(day,i)
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

gen_data <- function(){
    # Initialize states: 1 = S, 2 = I, 3 = R
    states <- matrix(1, nrow = n_individuals, ncol = n_days)
    states <- simulate_transitions(states, community_infection_rate, mean_duration_infected, n_days)
    
    # Convert to weekly observations
    sample_days <- seq(1, n_days, by = 7)
    weekly_states <- states[, seq(1, n_days, by = 7)]
    
    
    # Convert data to long format for model fitting
    simulation_data <- data.frame(
      individual = rep(1:n_individuals, times = n_weeks),
      day = rep(sample_days, each = n_individuals),
      state = as.vector(weekly_states)
    ) %>%
      left_join(ward_individual_mapping, by='individual') %>%
      arrange(individual, day) %>%
      group_by(individual) %>%
      filter(!(state == lag(state, default = first(state)) & state==3  )) %>% #censor if R
      ungroup() %>%
      mutate(
             ward_epidemic = if_else(
               (day>ward_epidemics[[1]][1] & day<ward_epidemics[[1]][2])|(day>ward_epidemics[[2]][1] & day<ward_epidemics[[2]][2]),
                                     1,0),
             sin365 = sin(2*pi*day/365),
             cos365 = cos(2*pi*day/365)
             
             )
}



    ####MARKOV MODEL 
    
fit_data <- function(X){    
    Q <- rbind(
      c(0, 1, 0),  # From S: S->I
      c(0, 0, 1 / mean_duration_infected),  # From I: I->R
      c(0, 0, 0)  # From R: no transitions
    )
    

    # Fit Markov transition model
    model <- msm(state ~ day , 
                 subject = individual,
                 data = X,
                 qmatrix = Q, 
                 covariates = list("1-2" = ~ ward_epidemic + sin365 +cos365),
                 gen.inits = TRUE)
    
    res <- hazard.msm(model)
    estimate <- res
    
    return(estimate)
}

#Generate 500 simulations
sims <- pbreplicate(nsim,gen_data(), simplify = F)

#extract results
all_estimates <- t(pblapply(sims,fit_data))

covar_estimates <- t(sapply(all_estimates,function(x){
     x$ward_epidemic["State 1 - State 2",]
}))

#mean of HR estimates across all
mean(covar_estimates[,'HR'])

#Power
mean(covar_estimates[,'L']>1)


mean(covar_estimates[,'L']<ward_epidemic_multiplier & covar_estimates[,'U']>ward_epidemic_multiplier)


prev <- lapply(1:length(sims), function(x){ 
  ds <-  sims[[x]]
  ds$sim = x 
  return(ds) 
}
  ) %>%
  bind_rows() %>%
  mutate(season = if_else(day<=365,1,2)
         ) %>%
  group_by(sim,season,individual) %>%
  mutate(pos = if_else(max(state==2)==1, 1,0 )
         ) %>%
  summarize(pos=max(pos)) %>%
  ungroup() %>%
  group_by(sim,season) %>%
  summarize(prev=mean(pos))

hist(prev$prev)
quantile(prev$prev, probs=c(0.5, 0.25, 0.75))
