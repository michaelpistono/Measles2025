

This individual-based model (IBM) simulates a measles outbreak using a SEIR framework 
(Susceptible, Exposed, Infectious, Recovered) with four age groups—preK, kindergarten, 5–17, and 18+—and age-assortative mixing via a 4×4 contact matrix.
It incorporates district-level population data, vaccination effects, and day-of-week variations in contact rates.

---

## Table of Contents

- [Overview](#overview)
- [Model Description](#model-description)
- [Requirements](#requirements)
- [Installation](#installation)
- [Running the Model](#running-the-model)
- [Simulation Parameters](#simulation-parameters)
- [Contact](#contact)

---

## Overview

This IBM simulates measles transmission dynamics in a county by tracking individuals through the stages of:
- **S**: Susceptible
- **E**: Exposed (latent)
- **I**: Infectious
- **R**: Recovered

The model assigns each individual to a district and age group based on district-level demographic data and vaccination exemption rates. It incorporates a time-varying transmission probability (with adjustments for vaccination ramp-up) and distinct contact rates on weekdays versus weekends.

---

## Model Description

- **Disease Dynamics:**
  - **Latent Period:** Modeled with a gamma distribution (mean ≈ 10 days, shape = 3).
  - **Infectious Period:** Modeled with a gamma distribution (mean ≈ 8 days, shape = 2). Infectiousness is considered active only for a limited period (up to 4.5 days) to simulate the appearance of rash.
  - **Vaccination:** Vaccine effectiveness (VE) reduces the force of infection for low-risk individuals.

- **Age-Assortative Mixing:**
  - The model uses a 4×4 contact matrix where contact rates vary by age group.
  - Weekday contact rates are set by parameters in the `theta` vector.
  - Weekend contact rates are fixed at higher values, reflecting increased social mixing during these days.

- **District-Level Data:**
  - Population data and age-specific parameters (e.g., percentage of children in each age bracket, exemption rates) are defined for several districts.
  - Each individual is randomly assigned a district, which then influences their age group and risk status.

---

## Requirements

- **R** (version 3.5 or higher is recommended)
- **Required R packages:**
  - `tidyverse`
  - `doParallel`
  - `foreach`
  - `doRNG` (for reproducible parallel random number generation)
  - `gridExtra`
  - `purrr`

---

## Installation

1. Install R and RStudio** (optional but recommended).
2. Install the required packages** by running the following in an R console:
   ```r
   install.packages(c("tidyverse", "doParallel", "foreach", "doRNG", "gridExtra", "purrr"))


## Running the Model

The model will:

1. Initialize a population with 22,500 susceptible and 2 infectious individuals.
2. Run the simulation for 365 days and replicate it 100 times using parallel processing.
3. Output time series data for overall prevalence, daily incidence, and cumulative cases.
4. Generate a plot comparing the simulated incidence with observed incidence data from Measles_gaines_observed.csv.

## Simulation Parameters
1. Initial State:
- Susceptible (S): 22,500
- Infectious (I): 2
  
2. Key Parameters (theta vector):
- beta: Transmission probability (e.g., 0.05 during early days; reduced after day 55 due to vaccination efforts).
- Duration: Mean infectious duration (8 days, drawn from Gamma distribution) - individuals moved out of "active" infections after 4.5 days assuming rash leads to isolation.
- latentPeriod: Mean latent period (10 days, drawn from Gamma distribution).
- VE: Vaccine effectiveness (97% after 2 doses of MMR vaccine).

3. Contact Matrix Parameters
4. Simulation Duration: 365 days (time steps)
5. Number of Simulations: 1000 replicates

## Contact
For further information, questions, or contributions, please contact:

[michael_pistono@berkeley.edu]

