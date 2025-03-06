##########################################################################
## IBM SEIR Model for Measles with 3 Age Groups (0-4, 5-17, 18+)
## using Age-Assortative Mixing (3x3 contact matrix)
##########################################################################

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
SEIR_ibm <- function(initState, theta, numIter) {
  
  ## General transmission parameters
  #beta in simulation loop allowing for varied transmission probability based on awareness of the outbreak
  ve           <- theta["ve"]           # Vaccine effectiveness (reduces force for low-risk)
  latentMean   <- theta["latentPeriod"] # Mean latent period (days)
  gamma_shape  <- 2
  gamma_scale  <- theta["D"] / gamma_shape
  
  ## County-Level Data (only the 5 counties from your original model)
  ## Columns:
  ##   - pctUnder5, pctUnder18: percentages, e.g. 7.3 => 7.3%
  ##   - kindExempt: "Kindergarten Vaccine Exempt" (for 0-4)
  ##   - k12Exempt: "Percent K12 Vaccine Exempt" (for 5-17)
  ##   - adultExempt: assumed fraction unvaccinated among 18+
  countyData <- data.frame(
    county      = c("Gaines"),
    population  = c(22500),
    pctUnder5   = c(10.4),
    pctUnder18  = c(36.0),
    kindExempt  = c(17.6),   # for 0-4
    k12Exempt   = c(13.6),   # for 5-17
    adultExempt = c(5.0)    # assumed 8% for 18+
  )
  
  totalCountyPop <- sum(countyData$population)
  countyData$weight <- countyData$population / totalCountyPop
  
  ## Initialize Population
  S0 <- initState["S"]
  I0 <- initState["I"]
  N  <- S0 + I0
  
  # Create a list to store individual attributes (state, times, county, age, risk).
  indiv <- vector(mode = "list", length = N)
  
  # Initialize states: first S0 as susceptible; next I0 as infectious.
  for (i in 1:S0) {
    indiv[[i]] <- list(state = "S", latentTime = NA, infectiousTime = NA)
  }
  for (i in (S0 + 1):N) {
    indiv[[i]] <- list(state = "I",
                       latentTime = NA,
                       infectiousTime = rgamma(1, shape = gamma_shape, scale = gamma_scale),
                       days_infectious = 1)
  }
  
  ## Assign County, Age Group, and Risk
  counties <- countyData$county
  weights  <- countyData$weight
  
  for (i in 1:N) {
    # 1) Pick county
    assignedCounty <- sample(counties, size = 1, prob = weights)
    indiv[[i]]$county <- assignedCounty
    countyParams <- countyData[countyData$county == assignedCounty, ]
    
    # 2) Decide which age group (0-4, 5-17, or 18+)
    r <- runif(1)
    frac_0_4  <- countyParams$pctUnder5 / 100
    frac_5_17 <- (countyParams$pctUnder18 - countyParams$pctUnder5) / 100
    # The remainder is 18+
    
    if (r < frac_0_4) {
      # 0-4
      indiv[[i]]$age <- "0-4"
      # Probability of being unvaccinated is 'kindExempt' from table
      if (runif(1) < (countyParams$kindExempt / 100)) {
        indiv[[i]]$risk <- "H"
      } else {
        indiv[[i]]$risk <- "L"
      }
      
    } else if (r < frac_0_4 + frac_5_17) {
      # 5-17
      indiv[[i]]$age <- "5-17"
      if (runif(1) < (countyParams$k12Exempt / 100)) {
        indiv[[i]]$risk <- "H"
      } else {
        indiv[[i]]$risk <- "L"
      }
      
    } else {
      # 18+
      indiv[[i]]$age <- "18+"
      if (runif(1) < (countyParams$adultExempt / 100)) {
        indiv[[i]]$risk <- "H"
      } else {
        indiv[[i]]$risk <- "L"
      }
    }
  }
  
  ## Pre-allocate Time Series
  timeVec <- 1:numIter
  numS <- numeric(numIter)
  numE <- numeric(numIter)
  numI <- numeric(numIter)
  numR <- numeric(numIter)
  
  stateVector <- sapply(indiv, function(x) x$state)
  numS[1] <- sum(stateVector == "S")
  numE[1] <- sum(stateVector == "E")
  numI[1] <- sum(stateVector == "I")
  numR[1] <- sum(stateVector == "R")
  
  ## Simulation Loop
  for (t in 1:(numIter - 1)) {
    
    
    if (t < 28.0) {
      beta <- theta["beta"] #transmission probability depending on time
    } else {
      beta <- 0.035
    }
    
    
    # --- Set contact rates based on day-of-week ---
    # t %% 7 == 6 or 0 => weekend
    if ((t %% 7 == 6) || (t %% 7 == 0)) {
      # On weekends, set all 3x3 contact rates to adjust for church services
      m_0_4_0_4  <- 30
      m_0_4_5_17 <- 30
      m_0_4_18   <- 30
      m_5_17_0_4 <- 30
      m_5_17_5_17<- 30
      m_5_17_18  <- 30
      m_18_0_4   <- 30
      m_18_5_17  <- 45
      m_18_18    <- 45
    } else {
      # Weekday contact rates from theta
      m_0_4_0_4  <- theta["m_0_4_0_4"]
      m_0_4_5_17 <- theta["m_0_4_5_17"]
      m_0_4_18   <- theta["m_0_4_18"]
      m_5_17_0_4 <- theta["m_5_17_0_4"]
      m_5_17_5_17<- theta["m_5_17_5_17"]
      m_5_17_18  <- theta["m_5_17_18"]
      m_18_0_4   <- theta["m_18_0_4"]
      m_18_5_17  <- theta["m_18_5_17"]
      m_18_18    <- theta["m_18_18"]
    }
    
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
    
    # Calculate prevalence by age group (actively infectious if days_infectious <= 4.5)
    total_0_4  <- sum(sapply(indiv, function(x) x$age == "0-4"))
    total_5_17 <- sum(sapply(indiv, function(x) x$age == "5-17"))
    total_18   <- sum(sapply(indiv, function(x) x$age == "18+"))
    
    active_infectious <- function(x) {
      x$state == "I" && !is.null(x$days_infectious) && x$days_infectious <= 4.5
    }
    
    prev_0_4  <- if (total_0_4  > 0) sum(sapply(indiv, function(x) x$age == "0-4"  && active_infectious(x))) / total_0_4  else 0
    prev_5_17 <- if (total_5_17 > 0) sum(sapply(indiv, function(x) x$age == "5-17" && active_infectious(x))) / total_5_17 else 0
    prev_18   <- if (total_18   > 0) sum(sapply(indiv, function(x) x$age == "18+"  && active_infectious(x))) / total_18   else 0
    
    # Update each individual
    for (j in 1:N) {
      currentState <- indiv[[j]]$state
      
      if (currentState == "S") {
        # Force of infection depends on individual's age group
        if (indiv[[j]]$age == "0-4") {
          lambda <- beta * (m_0_4_0_4  * prev_0_4 +
                              m_0_4_5_17 * prev_5_17 +
                              m_0_4_18   * prev_18)
        } else if (indiv[[j]]$age == "5-17") {
          lambda <- beta * (m_5_17_0_4  * prev_0_4 +
                              m_5_17_5_17 * prev_5_17 +
                              m_5_17_18   * prev_18)
        } else {  # "18+"
          lambda <- beta * (m_18_0_4   * prev_0_4 +
                              m_18_5_17  * prev_5_17 +
                              m_18_18    * prev_18)
        }
        
        # If "low-risk" (vaccinated), reduce force by vaccine effectiveness
        if (indiv[[j]]$risk == "L") {
          lambda <- (1 - ve) * lambda
        }
        
        # Infection event
        if (runif(1) < lambda) {
          indiv[[j]]$state <- "E"
          indiv[[j]]$latentTime <- rexp(1, rate = 1/latentMean)
        }
      }
      
      # E -> I (latent countdown)
      if (currentState == "E" && !is.na(indiv[[j]]$latentTime)) {
        indiv[[j]]$latentTime <- indiv[[j]]$latentTime - 1
        if (indiv[[j]]$latentTime <= 0) {
          indiv[[j]]$state <- "I"
          indiv[[j]]$latentTime <- NA
          indiv[[j]]$infectiousTime <- rgamma(1, shape = gamma_shape, scale = gamma_scale)
          indiv[[j]]$days_infectious <- 1
        }
      }
      
      # I -> R (infectious countdown)
      if (currentState == "I" && !is.na(indiv[[j]]$infectiousTime)) {
        indiv[[j]]$infectiousTime <- indiv[[j]]$infectiousTime - 1
        if (indiv[[j]]$infectiousTime <= 0) {
          indiv[[j]]$state <- "R"
          indiv[[j]]$infectiousTime <- NA
        }
      }
    }
    
    # Update compartment counts
    stateVector <- sapply(indiv, function(x) x$state)
    numS[t+1] <- sum(stateVector == "S")
    numE[t+1] <- sum(stateVector == "E")
    numI[t+1] <- sum(sapply(indiv, function(x) x$state == "I" &&
                              !is.null(x$days_infectious) &&
                              x$days_infectious <= 4.5))
    numR[t+1] <- sum(stateVector == "R")
  }
  
  ## Post-Processing
  prevOverall <- numI / (numS + numE + numI + numR)
  inc <- numI  # approximate incidence by number of actively infectious
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

# Initial state (S + I) from the new populations:
initState <- c(S = 22500, I = 2)

theta <- c(
  beta         = 0.0625,  # Baseline transmission probability per contact
  D            = 8,     # Infectious period (days)
  latentPeriod = 10,    # Latent period (days)
  ve           = 0.97,  # Vaccine effectiveness
  
  # 3x3 contact matrix (weekday rates):
  m_0_4_0_4   = 40.0,  # (0–4) to (0–4), avg school size with kindergarteners is 334 
  m_0_4_5_17  = 50.0,   # (0–4) to (5–17) #25 fits curve
  m_0_4_18    = 15.0,   # (0–4) to (18+)
  
  m_5_17_0_4  = 48.0,   # (5–17) to (0–4) avg school size with all other grades is 475
  m_5_17_5_17 = 95.0,  # (5–17) to (5–17) #65 fits curve
  m_5_17_18   = 25.0,   # (5–17) to (18+)
  
  m_18_0_4    = 5.0,    # (18+) to (0–4)
  m_18_5_17   = 20.0,   # (18+) to (5–17)
  m_18_18     = 15.0     # (18+) to (18+)
)

numIter <- 220
nSims   <- 50

## Run simulations in parallel
nCores <- parallel::detectCores() - 1
cl <- makeCluster(nCores)
registerDoParallel(cl)

simList <- foreach(i = 1:nSims, .packages = c("dplyr"), .options.RNG = 12345) %dorng% {
  SEIR_ibm(initState, theta, numIter)
}

stopCluster(cl)

## Summarize across replicates
avgSim <- bind_rows(simList) %>%
  group_by(time) %>%
  summarise(across(everything(), mean))

## ---------------------------
## Visualization
## ---------------------------


# If you have observed data to overlay:
obsData <- read.csv("Measles_gaines.csv") %>%
  mutate(Week = seq(1, nrow(.)),
         Day  = (Week * 7) - 7)

p2 <- ggplot() +
  geom_line(data = avgSim %>% filter(time <= 61),
            aes(x = time, y = inc, color = "Simulated Incidence"), size = 1.2) +
  geom_point(data = obsData,
             aes(x = Day, y = Total, color = "Observed Incidence"), size = 3) +
  labs(x = "Time (days)", y = "Incidence (new infections)",
       title = "Simulated vs. Observed Incidence in Gaines County, TX") +
  scale_color_manual(values = c("Simulated Incidence" = "blue", 
                                "Observed Incidence" = "red")) +
  theme_minimal()

p3 <- ggplot(avgSim, aes(x = time)) +
  geom_line(aes(y = inc, color = "Incidence"), size = 1.2) +
  labs(x = "Time (days)", y = "Incidence",
       title = "Model-Predicted Incidence in Measles Outbreak") +
  theme_minimal() +
  scale_color_manual(values = c("Incidence" = "blue"))

grid.arrange(p2, p3, ncol = 1)



