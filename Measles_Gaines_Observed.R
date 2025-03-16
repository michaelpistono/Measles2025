##########################################################################
## IBM SEIR Model for Measles with 4 Age Groups 
## (preK, kinder, 6-17, 18+) using Age-Assortative Mixing (4x4 contact matrix)
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

## ---------------------------
## Model Function: SEIR_ibm
## ---------------------------
SEIR_ibm <- function(initState, theta, numIter) {
  
  ## General transmission parameters
  # beta in simulation loop allows for a change in transmission probability based on outbreak awareness  
  ve           <- theta["ve"]           # Vaccine effectiveness (reduces force for low-risk)
  latentMean   <- theta["latentPeriod"] # Mean latent period (days)
  gamma_shape  <- 2
  gamma_scale  <- theta["D"] / gamma_shape
  
  ## District-Level Data & the remainder of the county population

  countyData <- data.frame(
    district    = c("Loop", "Seagraves", "Seminole", "Mennonite", "Remainder"),
    
    population  = c(151,        519,        2961,         420,        18449),
    pctUnder5   = c(0.0,        0.0,        0.0,          0.0,         11.8), #2,185 under 5/not enrolled in kindergarten yet (11.8% of 18449)
    pctKinder   = c(1.5,        7.6,        3.2,          10.0,         0.0),
    pctUnder18  = c(100.0,      100.0,      100.0,        100.0,      21.95), #4,049 people under 18 not in school/homeschooled (21.95% of 18449)
    
    preKExempt  = c(0.0,        0.0,        0.0,          0.0,         50.0),
    kinderExempt= c(53.85,      5.71,       7.86,         50.0,         0.0),
    k12Exempt   = c(47.95,      1.87,       13.8,         50.0,         9.8),
    adultExempt = c(0.0,        0.0,        0.0,          0.0,          8.0) #assume 8% of adults in county not vaccinated since mennonites came to state after vax push
    
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
                       infectiousTime = rgamma(1, shape = gamma_shape, scale = gamma_scale),
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
    # Here we assign them to the "5-17" age group and set risk = "H"
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
    
    if (t < 55.0) {  #people started vaccinating around day 40, and vaccine takes 10-14 days to become effective
      beta <- theta["beta"]  # transmission probability depending on time
    } else {
      beta <- 0.025 #beta reduced by half after new wave of vaccinations
    }
    
    # --- Set contact rates based on day-of-week ---
    # On weekends (t %% 7 == 6 or 0) we use fixed contact values; otherwise use theta parameters.
    if ((t %% 7 == 6) || (t %% 7 == 0)) {
      m_preK_preK    <- 15; m_preK_kinder  <- 15; m_preK_5_17   <- 15; m_preK_18   <- 12;
      m_kinder_preK  <- 12; m_kinder_kinder<- 15; m_kinder_5_17 <- 15; m_kinder_18 <- 12;
      m_5_17_preK    <- 12; m_5_17_kinder  <- 12; m_5_17_5_17   <- 25; m_5_17_18   <- 25;
      m_18_preK      <- 12; m_18_kinder    <- 12; m_18_5_17     <- 25; m_18_18     <- 30;
      
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
      x$state == "I" && !is.null(x$days_infectious) && x$days_infectious <= 4.5
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
          indiv[[j]]$latentTime <- rexp(1, rate = 1/latentMean)
        }
      }
      
      if (currentState == "E" && !is.na(indiv[[j]]$latentTime)) {
        
        indiv[[j]]$latentTime <- indiv[[j]]$latentTime - 1
        
        if (indiv[[j]]$latentTime <= 0) {
          indiv[[j]]$state <- "I"
          indiv[[j]]$latentTime <- NA
          indiv[[j]]$infectiousTime <- rgamma(1, shape = gamma_shape, scale = gamma_scale)
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

## ---------------------------
## Simulation Parameters & Run
## ---------------------------
initState <- c(S = 22500, I = 2)

theta <- c(
  beta         = 0.051,
  D            = 8,
  latentPeriod = 10,
  ve           = 0.97,
  
  ## Weekday contact matrix (4x4) parameters:
  m_preK_preK    = 1.8,
  m_preK_kinder  = 1.8,
  m_preK_5_17    = 1.8, #based on average household size
  m_preK_18      = 2.0,
  
  m_kinder_preK    = 2.0,
  m_kinder_kinder  = 17.0,
  m_kinder_5_17    = 5.0,
  m_kinder_18      = 10.0,
  
  m_5_17_preK    = 2.0,
  m_5_17_kinder  = 25.0, #average school size in Gaines County is ~ 400 students
  m_5_17_5_17    = 79.0,
  m_5_17_18      = 25.0,
  
  m_18_preK    = 2.0,
  m_18_kinder  = 5.0,
  m_18_5_17    = 5.0,
  m_18_18      = 24.0
)

numIter <- 250
nSims   <- 100

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

simQuant <- bind_rows(simList) %>%
  group_by(time) %>%
  summarise(across(c(inc, cumCases),
                   list(lower = ~quantile(., probs = 0.025),
                        upper = ~quantile(., probs = 0.975))))


## ---------------------------
## Visualization
## ---------------------------



obsData <- read.csv("Measles_gaines_observed.csv") %>%
  mutate(Week = seq(1, nrow(.)),
         Day  = (Week * 7) - 7)

p2 = ggplot(avgSim, aes(x = time)) +
  geom_line(aes(y = inc), color = "blue", size = 1.2) +
  geom_point(data = obsData,
             aes(x = Day, y = Total/7, color = "Observed Incidence"), size = 3) +
  labs(x = "Time (days)", y = "Incidence (new infections)",
       title = "Simulated vs. Observed Incidence in Gaines County, TX") +
  scale_color_manual(values = c("Simulated Incidence" = "blue", 
                                "Observed Incidence" = "red")) +
  theme_test()

grid.arrange(p2, ncol = 1)




simDF <- bind_rows(simList, .id = "sim")

p3 <- ggplot() +
  geom_line(data = simDF, aes(x = time, y = cumCases, group = sim),
            color = "gray", size = 0.6, alpha = 0.5) +
  geom_line(data = avgSim, aes(x = time, y = cumCases, color = "Mean Cumulative Cases"),
            size = 1.2) +
  geom_point(data = obsData,
             aes(x = Day, y = Cumulative, color = "Observed Cases"), size = 3) +
  labs(x = "Time (days)", y = "Cumulative Cases",
       title = "Simulated Cumulative Cases in Gaines County, TX") +
  scale_color_manual(name = "", 
                     values = c("Mean Cumulative Cases" = "blue",
                                "Observed Cases" = "red")) +
  theme_test()

grid.arrange(p3, ncol = 1)



###############################
## MLE: Estimating beta via MLE
###############################

# For computational speed during estimation, use a reduced number of simulation replicates.
nSims_MLE <- 25  # adjust as needed

# Initialize a cluster for the MLE block.
nCores_MLE <- parallel::detectCores() - 1
cl_mle <- makeCluster(nCores_MLE)

# Export necessary objects to the worker nodes.
clusterExport(cl_mle, varlist = c("SEIR_ibm", "initState", "theta", "numIter", "obsData"))
registerDoParallel(cl_mle)

# Define the negative log-likelihood function.
neg_log_lik <- function(log_beta) {
  
  # Transform parameter back (we optimize on log-scale)
  beta_val <- exp(log_beta)
  
  # Update theta with the candidate beta
  theta_fit <- theta
  theta_fit["beta"] <- beta_val
  
  # Run the simulation replicates with the candidate parameter
  simList_temp <- foreach(i = 1:nSims_MLE, .packages = c("dplyr"), .options.RNG = 12345) %dorng% {
    SEIR_ibm(initState, theta_fit, numIter)
  }
  
  # Compute the average simulation output (here we use cumulative cases)
  avgSim_temp <- bind_rows(simList_temp) %>%
    group_by(time) %>%
    summarise(cumCases = mean(cumCases))
  
  # Extract model-predicted cumulative cases at the observation days.
  pred <- approx(x = avgSim_temp$time, y = avgSim_temp$cumCases,
                 xout = obsData$Day, rule = 2)$y
  
  # Compute the negative log-likelihood assuming a Poisson distribution.
  nll <- -sum(dpois(obsData$Cumulative, lambda = pred, log = TRUE))
  
  return(nll)
}

# Use optim to minimize the negative log-likelihood.
init_log_beta <- log(theta["beta"])
opt_res <- optim(par = init_log_beta, fn = neg_log_lik, method = "BFGS")

# Get the estimated beta value
beta_est <- exp(opt_res$par)
cat("Estimated beta =", beta_est, "\n")

# Stop the cluster for the MLE block
stopCluster(cl_mle)
