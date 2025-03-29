##########################################################################
## IBM SEIR Model for Measles with 4 Age Groups 
## (preK, kindergarten, 6-17, 18+) using Age-Assortative Mixing (4x4 contact matrix)
##########################################################################

rm(list = ls())
setwd("~/Desktop/Measles Outbreak 2025")
set.seed(12345)

library(tidyverse)
library(doParallel)
library(foreach)
library(doRNG)      # For reproducible parallel RNG
library(gridExtra)
library(purrr)

################################
## Model Function: SEIR_ibm
################################

SEIR_ibm <- function(initState, theta, numIter) {
  
  ## General transmission parameters
  ve <- theta["VE"]           # Vaccine effectiveness (reduces force for low-risk)
  
  
  # Latent period follows a gamma distribution where most density falls between 5 and 12 days, mean of 10
  latent_shape <- 3                     
  latent_scale <- theta["latentPeriod"] / latent_shape
  
  # Infectious period parameters (also a gamma distribution where infectiousness peaks around 4 days with a long right tail)
  infectious_shape  <- 2
  infectious_scale  <- theta["Duration"] / infectious_shape
  
  ## District-Level Data & the remainder of the county population
  countyData <- data.frame(
    district    = c("Loop", "Seagraves", "Seminole", "Mennonite", "Remainder"),
    population  = c(151,        519,        2961,         420,        18449),
    
    pctUnder5   = c(0.0,        0.0,        0.0,          0.0,         11.8),
    pctKinder   = c(1.5,        7.6,        3.2,          10.0,         0.0),
    pctUnder18  = c(100.0,      100.0,      100.0,        100.0,      21.95),
    
    preKExempt  = c(0.0,        0.0,        0.0,          0.0,         65.0),
    kinderExempt= c(53.85,      5.71,       7.86,         50.0,         0.0),
    k12Exempt   = c(47.95,      1.87,       13.8,         50.0,         9.8),
    adultExempt = c(0.0,        0.0,        0.0,          0.0,          7.0) 
  )
  
  totalCountyPop <- sum(countyData$population)
  countyData$weight <- countyData$population / totalCountyPop
  
  ## Initialize Population
  S0 <- initState["S"]
  I0 <- initState["I"]
  N  <- S0 + I0
  
  # Create a list to store individual attributes (state, times, district, age, risk).
  indiv <- vector(mode = "list", length = N)
  
  # Initialize states: first S0 as susceptible; next I0 as infectious.
  for (i in 1:S0) {
    indiv[[i]] <- list(state = "S", latentTime = NA, infectiousTime = NA)
  }
  for (i in (S0 + 1):N) {
    indiv[[i]] <- list(state = "I",
                       latentTime = NA,
                       infectiousTime = rgamma(1, shape = infectious_shape, scale = infectious_scale),
                       days_infectious = 1)
  }
  
  ## Assign District, Age Group, and Risk
  districts <- countyData$district
  weights  <- countyData$weight
  
  for (i in 1:N) {
    # 1) Pick district
    assignedDistrict <- sample(districts, size = 1, prob = weights)
    indiv[[i]]$district <- assignedDistrict
    countyParams <- countyData[countyData$district == assignedDistrict, ]
    
    # For the initial infectious individuals, force assignment to unvaccinated children.
    if (i > S0) {
      indiv[[i]]$age <- "5-17"
      indiv[[i]]$risk <- "H"
    } else {
      # 2) For susceptibles, decide which age group (preK, kinder, 5-17, or 18+)
      r <- runif(1)
      frac_preK  <- countyParams$pctUnder5 / 100
      frac_kinder <- countyParams$pctKinder / 100
      frac_k12   <- (countyParams$pctUnder18 - countyParams$pctUnder5 - countyParams$pctKinder) / 100
      
      if (r < frac_preK) {
        indiv[[i]]$age <- "preK"
        if (runif(1) < (countyParams$preKExempt / 100)) {
          indiv[[i]]$risk <- "H"
        } else {
          indiv[[i]]$risk <- "L"
        }
      } else if (r < frac_preK + frac_kinder) {
        indiv[[i]]$age <- "kinder"
        if (runif(1) < (countyParams$kinderExempt / 100)) {
          indiv[[i]]$risk <- "H"
        } else {
          indiv[[i]]$risk <- "L"
        }
      } else if (r < frac_preK + frac_kinder + frac_k12) {
        indiv[[i]]$age <- "5-17"
        if (runif(1) < (countyParams$k12Exempt / 100)) {
          indiv[[i]]$risk <- "H"
        } else {
          indiv[[i]]$risk <- "L"
        }
      } else {
        indiv[[i]]$age <- "18+"
        if (runif(1) < (countyParams$adultExempt / 100)) {
          indiv[[i]]$risk <- "H"
        } else {
          indiv[[i]]$risk <- "L"
        }
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
  for (t in 0:(numIter - 1)) {
    
    if (t < 55.0) {  # MMR vaccinations started ramping up around day 40, and vaccine takes 10-14 days to become effective
      beta <- theta["beta"]  # transmission probability depending on time
    } else {
      beta <- 0.045 # beta reduced by half after new wave of vaccinations    was .025!!!
    }
    
    # Set contact rates based on day-of-week 
    if ((t %% 7 == 6) || (t %% 7 == 0)) {
      
      m_preK_preK    <- 2.5;m_preK_kinder  <- 5;  m_preK_5_17   <- 5;  m_preK_18   <- 5;
      m_kinder_preK  <- 5; m_kinder_kinder <- 5;  m_kinder_5_17 <- 10; m_kinder_18 <- 10;
      m_5_17_preK    <- 10; m_5_17_kinder  <- 10; m_5_17_5_17   <- 20; m_5_17_18   <- 20;
      m_18_preK      <- 10; m_18_kinder    <- 10; m_18_5_17     <- 20; m_18_18     <- 20;
      
    } else {
      
      m_preK_preK    <- theta["m_preK_preK"]
      m_preK_kinder  <- theta["m_preK_kinder"]
      m_preK_5_17    <- theta["m_preK_5_17"]
      m_preK_18      <- theta["m_preK_18"]
      
      m_kinder_preK    <- theta["m_kinder_preK"]
      m_kinder_kinder  <- theta["m_kinder_kinder"]
      m_kinder_5_17    <- theta["m_kinder_5_17"]
      m_kinder_18      <- theta["m_kinder_18"]
      
      m_5_17_preK    <- theta["m_5_17_preK"]
      m_5_17_kinder  <- theta["m_5_17_kinder"]
      m_5_17_5_17    <- theta["m_5_17_5_17"]
      m_5_17_18      <- theta["m_5_17_18"]
      
      m_18_preK    <- theta["m_18_preK"]
      m_18_kinder  <- theta["m_18_kinder"]
      m_18_5_17    <- theta["m_18_5_17"]
      m_18_18      <- theta["m_18_18"]
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
    total_preK    <- sum(sapply(indiv, function(x) x$age == "preK"))
    total_kinder  <- sum(sapply(indiv, function(x) x$age == "kinder"))
    total_5_17    <- sum(sapply(indiv, function(x) x$age == "5-17"))
    total_18      <- sum(sapply(indiv, function(x) x$age == "18+"))
    
    active_infectious <- function(x) {
      x$state == "I" && !is.null(x$days_infectious) && x$days_infectious <= 4.5 #remove infected individuals from FOI calculation after rash appears
    }
    
    prev_preK    <- if (total_preK > 0)    sum(sapply(indiv, function(x) x$age == "preK"    && active_infectious(x))) / total_preK    else 0
    prev_kinder  <- if (total_kinder > 0)  sum(sapply(indiv, function(x) x$age == "kinder"  && active_infectious(x))) / total_kinder  else 0
    prev_5_17    <- if (total_5_17 > 0)    sum(sapply(indiv, function(x) x$age == "5-17"    && active_infectious(x))) / total_5_17    else 0
    prev_18      <- if (total_18 > 0)      sum(sapply(indiv, function(x) x$age == "18+"     && active_infectious(x))) / total_18      else 0
    
    # Update each individual
    for (j in 1:N) {
      currentState <- indiv[[j]]$state
      
      if (currentState == "S") {
        if (indiv[[j]]$age == "preK") {
          lambda <- beta * (m_preK_preK   * prev_preK +
                              m_preK_kinder * prev_kinder +
                              m_preK_5_17   * prev_5_17 +
                              m_preK_18     * prev_18)
        } else if (indiv[[j]]$age == "kinder") {
          lambda <- beta * (m_kinder_preK   * prev_preK +
                              m_kinder_kinder * prev_kinder +
                              m_kinder_5_17   * prev_5_17 +
                              m_kinder_18     * prev_18)
        } else if (indiv[[j]]$age == "5-17") {
          lambda <- beta * (m_5_17_preK   * prev_preK +
                              m_5_17_kinder * prev_kinder +
                              m_5_17_5_17   * prev_5_17 +
                              m_5_17_18     * prev_18)
        } else {  # "18+"
          lambda <- beta * (m_18_preK   * prev_preK +
                              m_18_kinder * prev_kinder +
                              m_18_5_17   * prev_5_17 +
                              m_18_18     * prev_18)
        }
        
        if (indiv[[j]]$risk == "L") {
          lambda <- (1 - ve) * lambda
        }
        
        if (runif(1) < lambda) {
          indiv[[j]]$state <- "E"
          # Draw the latent period from a gamma distribution
          indiv[[j]]$latentTime <- rgamma(1, shape = latent_shape, scale = latent_scale)
        }
      }
      
      if (currentState == "E" && !is.na(indiv[[j]]$latentTime)) {
        indiv[[j]]$latentTime <- indiv[[j]]$latentTime - 1
        
        if (indiv[[j]]$latentTime <= 0) {
          indiv[[j]]$state <- "I"
          indiv[[j]]$latentTime <- NA
          indiv[[j]]$infectiousTime <- rgamma(1, shape = infectious_shape, scale = infectious_scale)
          indiv[[j]]$days_infectious <- 1
        }
      }
      
      if (currentState == "I" && !is.na(indiv[[j]]$infectiousTime)) {
        indiv[[j]]$infectiousTime <- indiv[[j]]$infectiousTime - 1
        if (indiv[[j]]$infectiousTime <= 0) {
          indiv[[j]]$state <- "R"
          indiv[[j]]$infectiousTime <- NA
        }
      }
    }
    
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
  inc <- numI  # approximate incidence by number of actively infectious individuals
  cumCases <- S0 - numS
  
  results <- data.frame(time = timeVec,
                        prevOverall = 100 * prevOverall,
                        inc = inc,
                        cumCases = cumCases)
  
  return(results)
}

#########################################
## Simulation Parameters & Run
#########################################
initState <- c(S = 22500, I = 2)

theta <- c(
  beta = 0.05,
  Duration = 8,   # mean infectious duration
  latentPeriod = 10,  # mean latent period
  VE = 0.97, # VE after 2 doses of MMR vaccine
  
  ## Weekday contact matrix (4x4) parameters:
  m_preK_preK    = 1.8,
  m_preK_kinder  = 1.8,
  m_preK_5_17    = 1.8, # based on average household size
  m_preK_18      = 2.0,
  
  m_kinder_preK    = 2.0,
  m_kinder_kinder  = 17.0,
  m_kinder_5_17    = 5.0,
  m_kinder_18      = 10.0,
  
  m_5_17_preK    = 2.0,
  m_5_17_kinder  = 20.0, # average school size in Gaines County is ~ 400 students
  m_5_17_5_17    = 72.0,
  m_5_17_18      = 20.0,
  
  m_18_preK    = 2.0,
  m_18_kinder  = 5.0,
  m_18_5_17    = 5.0,
  m_18_18      = 20.0
)

numIter <- 365  # number of time steps/days per simulation
nSims   <- 1000   # number of replicates/simulations


nCores <- parallel::detectCores() - 1   #parallel processing so this doesn't take all year
cl <- makeCluster(nCores)
registerDoParallel(cl)

simList <- foreach(i = 1:nSims, .packages = c("dplyr"), .options.RNG = 12345) %dorng% {
  SEIR_ibm(initState, theta, numIter)
}

stopCluster(cl)


avgSim <- bind_rows(simList) %>%    # calculate averages across all simulations
  group_by(time) %>%
  summarise(across(everything(), mean))

simQuant <- bind_rows(simList) %>%    # calculate 2.5 and 97.5 quantiles for estimates
  group_by(time) %>%
  summarise(across(c(inc, cumCases),
                   list(lower = ~quantile(., probs = 0.025),
                        upper = ~quantile(., probs = 0.975))))

####################################
## Visualization
####################################

obsData <- read.csv("Measles_gaines_observed.csv") %>% # read in data from csv file
  mutate(Week = seq(1, nrow(.)),
         Day  = (Week * 7) - 7) # ensure days start at 0 instead of 7

simDF <- bind_rows(simList, .id = "sim") # move cumulative data into long form for plotting

p2 = ggplot(avgSim, aes(x = time)) +
  geom_line(data = simDF, aes(x = time, y = inc, group = sim),
            color = "gray", size = 0.6, alpha = 0.5) +
  geom_line(aes(y = inc, color = "Mean Simulated Incidence"), size = 1.2) +
  geom_point(data = obsData,
             aes(x = Day, y = Total/7, color = "Observed Incidence"), size = 3) +
  
  labs(x = "Time (days)", y = "Incidence (new infections)",
       title = "Simulated vs. Observed Incidence in Gaines County, TX") +
  scale_color_manual(values = c("Mean Simulated Incidence" = "blue",
                                "Observed Incidence" = "red")) +
  theme_test()

grid.arrange(p2, ncol = 1)



