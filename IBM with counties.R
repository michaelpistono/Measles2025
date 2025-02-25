##########################################################################
## IBM SEIR Model for Measles with Age Groups
## Incorporating County-Specific Vaccine Exemption and 
## Age-Specific Contact Rates (Children vs. Adults; High-risk vs. Low-risk)
## Infectious individuals are assumed to stop contributing to transmission 
## 5 days after becoming infectious.
##
## This script simulates measles transmission in a heterogeneous population 
## using a stochastic SEIR individual-based model. Individuals are assigned 
## to counties (with county-specific K-12 proportions and vaccine exemption 
## rates) and categorized as either "child" (K-12 aged) or "adult." For each 
## age group, separate contact parameters are used for high-risk and low-risk 
## individuals. Once an individual becomes infectious, they contribute to 
## transmission only for the first 5 days.
##
## County Data:
##   Gaines County, TX:   Population = 22,500; 35% K-12; 17.6% vaccine exempt (children)
##   Dawson County, TX:   Population = 12,004; 25.6% K-12; 2.52% vaccine exempt (children)
##   Terry County, TX:    Population = 11,576; 28.4% K-12; 2.68% vaccine exempt (children)
##   Yoakum County, TX:   Population = 7,468;  32.1% K-12; 1.9% vaccine exempt (children)
##   Lea County, NM:      Population = 72,000; 28.8% K-12; 0.9% vaccine exempt (children)
##
## For adults, the assumed rate of non-MMR immunization is 13% for Gaines and 5% for others.
##
## Author: Mike Pistono
## Date: 24 Feb 2025
##########################################################################

## ---------------------------
## Initialize Environment
## ---------------------------
rm(list = ls())
setwd("~/Desktop/Measles Outbreak 2025")
set.seed(12345)

library(tidyverse)
library(doParallel)
library(foreach)
library(doRNG)      # For reproducible parallel RNG
library(gridExtra)

## ---------------------------
## Model Function: SEIR_ibm
## ---------------------------
# This function runs the individual-based SEIR model.
#
# Args:
#   initState: Named vector with initial counts (S: susceptible, I: infectious).
#   theta: Named vector of parameters, including:
#          beta           - baseline transmission probability per contact,
#          D              - infectious period (days),
#          latentPeriod   - latent period (days),
#          ve             - vaccine effectiveness,
#          cMean_k12      - average contacts per day for K-12 aged individuals,
#          cLow_k12       - contacts per day for low-risk (vaccinated) K-12 individuals,
#          cMean_adult    - average contacts per day for adults,
#          cLow_adult     - contacts per day for low-risk (vaccinated) adults.
#   numIter: Number of simulation time steps (days).
#
# Note:
#   County-specific parameters (population, K-12 proportion, exemption rates) are defined inside.
#   When an individual transitions from E → I, we initialize a counter "days_infectious."
#   Only individuals with days_infectious ≤ 4 contribute to transmission since individuals can transmit 4 days before rash.
#
# Returns:
#   A data frame with overall compartment counts, prevalence, incidence, and total cases over time.
SEIR_ibm <- function(initState, theta, numIter) {
  
  ## General transmission parameters --------------------------------
  beta         <- theta["beta"]
  ve           <- theta["ve"]
  latentMean   <- theta["latentPeriod"]
  gamma_shape  <- 2
  gamma_scale  <- theta["D"] / gamma_shape
  
  ## Age-specific contact parameters --------------------------------
  cMean_k12   <- theta["cMean_k12"]
  cLow_k12    <- theta["cLow_k12"]
  cMean_adult <- theta["cMean_adult"]
  cLow_adult  <- theta["cLow_adult"]
  
  ## County-Level Data --------------------------------
  countyData <- data.frame(
    county     = c("Gaines", "Dawson", "Terry", "Yoakum", "Lea"),
    population = c(22500, 12004, 11576, 7468, 72000),
    propSchool = c(0.35, 0.256, 0.284, 0.321, 0.288),  # proportion K-12 aged
    riskChild  = c(0.18, 0.025, 0.027, 0.02, 0.01), # % vaccine exemptions among k-12
    riskAdult  = c(0.13, 0.05, 0.05, 0.05, 0.05)   # estimated % unvaccinated adults
  )
  totalCountyPop <- sum(countyData$population)
  countyData$weight <- countyData$population / totalCountyPop
  
  ## Initialize Population --------------------------------
  S0 <- initState["S"]
  I0 <- initState["I"]
  N  <- S0 + I0
  
  # Create a list to store individual attributes.
  indiv <- vector(mode = "list", length = N)
  
  # Initialize individuals: first S0 as susceptible; next I0 as infectious.
  for (i in 1:S0) {
    indiv[[i]] <- list(state = "S", latentTime = NA, infectiousTime = NA)
  }
  for (i in (S0 + 1):N) {
    # Initialize infectious individuals with days_infectious = 1.
    indiv[[i]] <- list(state = "I", latentTime = NA,
                       infectiousTime = rgamma(1, shape = gamma_shape, scale = gamma_scale),
                       days_infectious = 1)
  }
  
  ## Assign County, Age, and Risk Status --------------------------------
  # For each individual, sample a county based on its population weight.
  # Then assign age: if a random draw is less than the county's propSchool then "child" (K-12) else "adult."
  # Risk status is determined using county-specific exemption rates.
  counties <- countyData$county
  weights <- countyData$weight
  
  for (i in 1:N) {
    assignedCounty <- sample(counties, size = 1, prob = weights)
    indiv[[i]]$county <- assignedCounty
    countyParams <- countyData[countyData$county == assignedCounty, ]
    
    if (runif(1) < countyParams$propSchool) {
      indiv[[i]]$age <- "child"
      indiv[[i]]$risk <- ifelse(runif(1) < countyParams$riskChild, "H", "L")
    } else {
      indiv[[i]]$age <- "adult"
      indiv[[i]]$risk <- ifelse(runif(1) < countyParams$riskAdult, "H", "L")
    }
  }
  
  ## Pre-allocate Storage for Time Series --------------------------------
  timeVec <- 1:numIter
  numS <- numeric(numIter); numE <- numeric(numIter)
  numI <- numeric(numIter); numR <- numeric(numIter)
  
  # Record initial overall compartment counts.
  stateVector <- sapply(indiv, function(x) x$state)
  numS[1] <- sum(stateVector == "S")
  numE[1] <- sum(stateVector == "E")
  numI[1] <- sum(stateVector == "I")
  numR[1] <- sum(stateVector == "R")
  
  ## Simulation Loop --------------------------------
  for (t in 1:(numIter - 1)) {
    
    # Update days_infectious for those in state "I"
    for (j in 1:N) {
      if (indiv[[j]]$state == "I") {
        if (is.null(indiv[[j]]$days_infectious)) {
          indiv[[j]]$days_infectious <- 1
        } else {
          indiv[[j]]$days_infectious <- indiv[[j]]$days_infectious + 1
        }
      }
    }
    
    # Subgroup counts and "active" infectious prevalence --------------------------------
    # only count infectious individuals with days_infectious <= 4.5
    total_k12 <- sum(sapply(indiv, function(x) x$age == "child"))
    total_adult <- sum(sapply(indiv, function(x) x$age == "adult"))
    
    total_k12_H <- sum(sapply(indiv, function(x) x$age == "child" && x$risk == "H"))
    total_k12_L <- total_k12 - total_k12_H
    total_adult_H <- sum(sapply(indiv, function(x) x$age == "adult" && x$risk == "H"))
    total_adult_L <- total_adult - total_adult_H
    
    active_infectious <- function(x) {
      return(x$state == "I" && !is.null(x$days_infectious) && x$days_infectious <= 4.5)
    }
    
    prev_k12_H <- if (total_k12_H > 0) {
      sum(sapply(indiv, function(x) active_infectious(x) && x$age == "child" && x$risk == "H")) / total_k12_H
    } else 0
    prev_k12_L <- if (total_k12_L > 0) {
      sum(sapply(indiv, function(x) active_infectious(x) && x$age == "child" && x$risk == "L")) / total_k12_L
    } else 0
    prev_adult_H <- if (total_adult_H > 0) {
      sum(sapply(indiv, function(x) active_infectious(x) && x$age == "adult" && x$risk == "H")) / total_adult_H
    } else 0
    prev_adult_L <- if (total_adult_L > 0) {
      sum(sapply(indiv, function(x) active_infectious(x) && x$age == "adult" && x$risk == "L")) / total_adult_L
    } else 0
    
    # Proportions of high-risk/unvaccinated within each age group
    p_k12_H <- if (total_k12 > 0) total_k12_H / total_k12 else 0
    p_adult_H <- if (total_adult > 0) total_adult_H / total_adult else 0
    
    # Compute effective contact rates for high-risk groups in each age group.
    cHigh_k12 <- if (p_k12_H > 0) (cMean_k12 - cLow_k12 * (1 - p_k12_H)) / p_k12_H else cMean_k12
    cHigh_adult <- if (p_adult_H > 0) (cMean_adult - cLow_adult * (1 - p_adult_H)) / p_adult_H else cMean_adult
    
    # Mixing weights and mixing prevalence for each age group:
    denom_k12 <- if (total_k12 > 0) {
      (cHigh_k12 * (total_k12_H/total_k12) + cLow_k12 * (total_k12_L/total_k12))
    } else 0
    gH_k12 <- if (denom_k12 > 0) (cHigh_k12 * (total_k12_H/total_k12)) / denom_k12 else 0
    pMix_k12 <- gH_k12 * prev_k12_H + (1 - gH_k12) * prev_k12_L
    
    denom_adult <- if (total_adult > 0) {
      (cHigh_adult * (total_adult_H/total_adult) + cLow_adult * (total_adult_L/total_adult))
    } else 0
    gH_adult <- if (denom_adult > 0) (cHigh_adult * (total_adult_H/total_adult)) / denom_adult else 0
    pMix_adult <- gH_adult * prev_adult_H + (1 - gH_adult) * prev_adult_L
    
    # Force of infection for each subgroup:
    lambda_k12_H <- cHigh_k12 * beta * pMix_k12
    lambda_k12_L <- cLow_k12  * beta * pMix_k12
    lambda_adult_H <- cHigh_adult * beta * pMix_adult
    lambda_adult_L <- cLow_adult  * beta * pMix_adult
    
    # Update each individual --------------------------------
    for (j in 1:N) {
      currentState <- indiv[[j]]$state
      if (currentState == "S") {
        if (indiv[[j]]$age == "child") {
          if (indiv[[j]]$risk == "H") {
            lambda <- lambda_k12_H
          } else {
            lambda <- (1 - ve) * lambda_k12_L
          }
        } else {  # adult
          if (indiv[[j]]$risk == "H") {
            lambda <- lambda_adult_H
          } else {
            lambda <- (1 - ve) * lambda_adult_L
          }
        }
        if (runif(1) < lambda) {
          indiv[[j]]$state <- "E"
          indiv[[j]]$latentTime <- rexp(1, rate = 1/latentMean)
        }
      }
      
      # Transition E -> I (latent period countdown)
      if (currentState == "E" && !is.na(indiv[[j]]$latentTime)) {
        indiv[[j]]$latentTime <- indiv[[j]]$latentTime - 1
        if (indiv[[j]]$latentTime <= 0) {
          indiv[[j]]$state <- "I"
          indiv[[j]]$latentTime <- NA
          indiv[[j]]$infectiousTime <- rgamma(1, shape = gamma_shape, scale = gamma_scale)
          indiv[[j]]$days_infectious <- 1  # start the counter
        }
      }
      
      # Transition I -> R (infectious period countdown)
      if (currentState == "I" && !is.na(indiv[[j]]$infectiousTime)) {
        indiv[[j]]$infectiousTime <- indiv[[j]]$infectiousTime - 1
        if (indiv[[j]]$infectiousTime <= 0) {
          indiv[[j]]$state <- "R"
          indiv[[j]]$infectiousTime <- NA
        }
      }
    }
    
    # Update overall compartment counts (aggregated over all groups)
    stateVector <- sapply(indiv, function(x) x$state)
    numS[t+1] <- sum(stateVector == "S")
    numE[t+1] <- sum(stateVector == "E")
    # Here, incidence is taken as the count of infectious individuals who are still "active" (days_infectious <= 4.5)
    numI[t+1] <- sum(sapply(indiv, function(x) x$state == "I" && !is.null(x$days_infectious) && x$days_infectious <= 4.5))
    numR[t+1] <- sum(stateVector == "R")
  }
  
  ## Post-Processing --------------------------------
  prevOverall <- numI / (numS + numE + numI + numR)
  inc <- numI  # Incidence approximated by the number of "active" infectious individuals
  totalCases <- S0 - numS[numIter]
  
  results <- data.frame(time = timeVec,
                        prevOverall = 100 * prevOverall,
                        inc = inc,
                        totalCases = totalCases)
  
  return(results)
}

## ---------------------------
## Simulation Parameters & Run
## ---------------------------
# Define initial state using the sum of county populations:
# Gaines: 22,500; Dawson: 12,004; Terry: 11,576; Yoakum: 7,468; Lea: 72000
initState <- c(S = 22500 + 12004 + 11576 + 7468 + 72000, I = 3) #3 initial infections seeding outbreak

# Parameter vector (general transmission parameters and age-specific contact rates)
theta <- c(beta = 0.06,           # Baseline transmission probability per contact
           D = 8,                 # Infectious period (days)
           latentPeriod = 10,     # Latent period (days)
           ve = 0.96,             # Vaccine effectiveness
           cMean_k12 = 250,        # Average contacts/day for K-12 aged individuals (avg school size is 409 students in Gaines County)
           cLow_k12 = 250,         
           cMean_adult = 15,      # Average contacts/day for adults
           cLow_adult = 15)       

numIter <- 160    # Simulate for 160 days
nSims   <- 50     # Number of simulation replicates

## Run simulations in parallel ---------------------------------
nCores <- parallel::detectCores() - 1
cl <- makeCluster(nCores)
registerDoParallel(cl)

simList <- foreach(i = 1:nSims, .packages = c("dplyr"), .options.RNG = 12345) %dorng% {
  SEIR_ibm(initState, theta, numIter)
}

stopCluster(cl)

avgSim <- bind_rows(simList) %>%
  group_by(time) %>%
  summarise(across(everything(), mean))

## ---------------------------
## Visualization
## ---------------------------


# Plot incidence over time
p2 <- ggplot(avgSim, aes(x = time)) +
  geom_line(aes(y = inc, color = "Incidence"), size = 1.2) +
  labs(x = "Time (days)", y = "Incidence", 
       title = "Model-Predicted Incidence in Measles Outbreak") +
  theme_minimal() +
  scale_color_manual(values = c("Incidence" = "blue"))

# Overlay with observed outbreak data if available
obsData <- read.csv("Measles.csv") %>%
  mutate(Week = seq(1, nrow(.)),
         Day  = Week * 7)

p3 <- ggplot() +
  geom_line(data = avgSim %>% filter(time <= 50),
            aes(x = time, y = inc, color = "Simulated Incidence"), size = 1.2) +
  geom_point(data = obsData,
             aes(x = Day, y = Total, color = "Observed Incidence"), size = 3) +
  labs(x = "Time (days)", y = "Incidence (new infections)", 
       title = "Simulated vs. Observed Incidence (First 5 Weeks)") +
  scale_color_manual(values = c("Simulated Incidence" = "blue", "Observed Incidence" = "red")) +
  theme_minimal()

grid.arrange(p2, p3, ncol = 1)


