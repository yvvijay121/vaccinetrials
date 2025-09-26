#######################################################################
######################## SINGLE RUN OF PREDICTEE ALGORITHM #############
#######################################################################

# Load required libraries
library(magrittr)
library(survival)
library(survminer)
library(dplyr)
library(pwr)
library(tidyverse)
library(randomForestSRC)

# Set working directory to current directory (adjust path as needed)
# Note: Original script used Windows path, adjusted for current directory
setwd("~/Documents/Programming_Projects/vaccinetrials")

#######################################################################
######################## LOAD DATA AND MODELS #########################
#######################################################################

# Load pre-trained models
cox_model <- readRDS("universal_cox.rds")
rsf_model <- readRDS("universal_rsf.rds")

# Check what variables the Cox model expects
print("Cox model formula:")
print(formula(cox_model))
print("Cox model variables:")
print(names(cox_model$coefficients))

# Load CNEP data for target demographics
cnep <- read.csv("data/cnep_plus_all_2018.02.13.csv")
cnep_susceptible <- cnep[cnep$HCV == "susceptible",]

# Stratify age categories in CNEP data
cnep_susceptible$Age_Category <- NA
cnep_susceptible[cnep_susceptible$Age < 20,]$Age_Category <- "<20"
cnep_susceptible[cnep_susceptible$Age >= 20 & cnep_susceptible$Age < 30,]$Age_Category <- "20-29"
cnep_susceptible[cnep_susceptible$Age >= 30 & cnep_susceptible$Age < 40,]$Age_Category <- "30-39"
cnep_susceptible[cnep_susceptible$Age >= 40 & cnep_susceptible$Age < 50,]$Age_Category <- "40-49"
cnep_susceptible[cnep_susceptible$Age >= 50,]$Age_Category <- "49+"

# Calculate target demographic proportions
target_pop <- cnep_susceptible
target_gender_comp <- matrix(table(target_pop$Gender)/nrow(target_pop), nrow = 1)
colnames(target_gender_comp) <- c("Female", "Male")
target_race_comp <- matrix(table(target_pop$Race)/nrow(target_pop), nrow = 1)
colnames(target_race_comp) <- c("Hispanic", "NHBlack", "NHWhite", "Other")
target_age_comp <- matrix(table(target_pop$Age_Category)/nrow(target_pop), nrow = 1)
colnames(target_age_comp) <- c("<20", "20-29", "30-39", "40-49", "49+")

# Combine target demographics (gender, race, and age)
target_demo <- cbind(target_gender_comp, target_race_comp, target_age_comp)

#######################################################################
###################### CREATE RECRUITMENT POOL ######################
#######################################################################

# Note: This assumes you have events_pool.csv file with simulation data
# If you don't have this file, you'll need to generate or obtain it
simulation_data <- read.csv("data/events_pool.csv")

# Create recruitment pool following the same process as in recruitment_pool.R
# Label data with infection times
agents_infected_total <- simulation_data[simulation_data$Event == "infected",]
agents_infected_total <- agents_infected_total[!duplicated(agents_infected_total$Agent), c("Time", "Agent")]
colnames(agents_infected_total) <- c("infected_time", "Agent")

# Create recruitment pool
recruitment_pool <- merge(simulation_data, agents_infected_total, by = "Agent", all.x = TRUE)
recruitment_pool <- recruitment_pool[recruitment_pool$Event == "activated",]

# Label infection status within trial period
recruitment_pool$infected_time_after_start <- recruitment_pool$infected_time - recruitment_pool$Time
recruitment_pool[is.na(recruitment_pool$infected_time_after_start),]$infected_time_after_start <- 999
recruitment_pool$infected_by_trialend <- 0
recruitment_pool[recruitment_pool$infected_time_after_start <= 1.5, ]$infected_by_trialend <- 1

# Add age categories
recruitment_pool$Age_Category <- NA
recruitment_pool[recruitment_pool$Age < 20,]$Age_Category <- "<20"
recruitment_pool[recruitment_pool$Age >= 20 & recruitment_pool$Age < 30,]$Age_Category <- "20-29"
recruitment_pool[recruitment_pool$Age >= 30 & recruitment_pool$Age < 40,]$Age_Category <- "30-39"
recruitment_pool[recruitment_pool$Age >= 40 & recruitment_pool$Age < 50,]$Age_Category <- "40-49"
recruitment_pool[recruitment_pool$Age >= 50,]$Age_Category <- "49+"

# Convert variables to factors with correct levels that match the Cox model
recruitment_pool$Gender <- factor(recruitment_pool$Gender, levels = c("Female", "Male"))
recruitment_pool$Race <- factor(recruitment_pool$Race)
recruitment_pool$Syringe_source <- factor(recruitment_pool$Syringe_source, levels = c("HR", "nonHR"))
recruitment_pool$chicago_community_name <- factor(recruitment_pool$chicago_community_name)
recruitment_pool$Age_Category <- factor(recruitment_pool$Age_Category)

# Label susceptible status
recruitment_pool$susceptible <- 0
recruitment_pool$susceptible[recruitment_pool$HCV == "susceptible"] <- 1

# Check what variables we have in recruitment_pool
print("Variables in recruitment_pool:")
print(colnames(recruitment_pool))

# Ensure recruitment_pool has all variables needed by the Cox model
# Based on the actual Cox model, it does NOT include chicago_community_name
required_vars <- c("Age", "Gender", "Syringe_source", "Drug_in_degree", "Drug_out_degree", 
                   "current_total_network_size", "Daily_injection_intensity", 
                   "Fraction_recept_sharing")

missing_vars <- required_vars[!required_vars %in% colnames(recruitment_pool)]
if(length(missing_vars) > 0) {
  print("Missing variables for Cox model:")
  print(missing_vars)
  stop("Required variables missing from recruitment pool")
}

print("All required variables present!")

# CRITICAL FIX: The Cox model has an environment issue where it's looking for 'unique_agents'
# We need to create a workaround for the predict function
print("Creating workaround for Cox model environment issue...")

# Create a minimal unique_agents object with the right structure and correct factor levels
unique_agents_sample <- recruitment_pool[1:10, required_vars]
unique_agents_sample$Gender <- factor(unique_agents_sample$Gender, levels = c("Female", "Male"))
unique_agents_sample$Syringe_source <- factor(unique_agents_sample$Syringe_source, levels = c("HR", "nonHR"))
unique_agents_sample$status <- c(rep(0, 5), rep(1, 5)) 
unique_agents_sample$survival_time <- rep(1.5, 10)

# Assign to global environment so the Cox model can find it
assign("unique_agents", unique_agents_sample, envir = .GlobalEnv)

print("Workaround applied!")

#######################################################################
####################### LOAD ALGORITHM FUNCTION #####################
#######################################################################

# Source the algorithm function from batch_algorithm.R
source("batch_algorithm.R")

#######################################################################
########################## RUN ALGORITHM ONCE ########################
#######################################################################

# Set algorithm parameters
recruitment_per_batch <- 20    # B: number of PWID screened per batch
recruited_per_batch <- 10      # R: number of PWID recruited per batch
trial_followup_years <- 1.5    # Follow-up period for trial
req_sample_size <- 500          # Initial estimated sample size
work_constraint <- 2000         # Maximum number of PWID that can be screened

# Define target demographics for matching
target_demographics <- c("Gender", "Race", "Age_Category")

print("Starting PREDICTEE algorithm with Cox model...")
print(paste("Target sample size:", req_sample_size))
print(paste("Recruitment per batch:", recruitment_per_batch))
print(paste("Recruited per batch:", recruited_per_batch))
print(paste("Work constraint:", work_constraint))

# Run algorithm with Cox model
cox_results <- algorithm(
  recruitment_dataset = recruitment_pool,
  model_used = cox_model,
  target_demographics = target_demographics,
  target_props = target_demo,
  recruitment_per_batch = recruitment_per_batch,
  recruited_per_batch = recruited_per_batch,
  trial_followup_years = trial_followup_years,
  req_sample_size = req_sample_size,
  work_constraint = work_constraint,
  incidence_weight_min = 25,
  attrition_prob = 0,
  high_demo_error_adjustment = FALSE,
  ssmethod = "cohen",
  print_diagnostics = TRUE
)

print("Cox model algorithm completed!")
print(paste("Final recruited sample size:", nrow(cox_results$algorithm_output)))
print(paste("Total agents screened:", cox_results$agent_count))

# Display demographic composition of recruited sample
recruited_sample <- cox_results$algorithm_output

print("\n=== DEMOGRAPHIC COMPOSITION OF RECRUITED SAMPLE ===")
print("Gender distribution:")
print(table(recruited_sample$Gender))

print("\nRace distribution:")
print(table(recruited_sample$Race))

print("\nAge category distribution:")
print(table(recruited_sample$Age_Category))

print("\n=== COMPARISON WITH TARGET DEMOGRAPHICS ===")
print("Target gender proportions:")
print(target_gender_comp)
print("Actual gender proportions:")
print(table(recruited_sample$Gender)/nrow(recruited_sample))

print("\nTarget race proportions:")
print(target_race_comp)
print("Actual race proportions:")
print(table(recruited_sample$Race)/nrow(recruited_sample))

print("\nTarget age proportions:")
print(target_age_comp)
print("Actual age proportions:")
print(table(recruited_sample$Age_Category)/nrow(recruited_sample))

# Calculate mean predicted infection probability
mean_infection_prob <- mean(recruited_sample$infected_probability)
print(paste("\nMean predicted infection probability:", round(mean_infection_prob, 4)))

# Optional: Run with RSF model as well (uncomment if desired)
# print("\n\nStarting PREDICTEE algorithm with RSF model...")
# rsf_results <- algorithm(
#   recruitment_dataset = recruitment_pool,
#   model_used = rsf_model,
#   target_demographics = target_demographics,
#   target_props = target_demo,
#   recruitment_per_batch = recruitment_per_batch,
#   recruited_per_batch = recruited_per_batch,
#   trial_followup_years = trial_followup_years,
#   req_sample_size = req_sample_size,
#   work_constraint = work_constraint,
#   print_diagnostics = TRUE
# )
# 
# print("RSF model algorithm completed!")
# print(paste("Final recruited sample size:", nrow(rsf_results$algorithm_output)))
# print(paste("Total agents screened:", rsf_results$agent_count))

print("\n=== ALGORITHM EXECUTION COMPLETE ===")