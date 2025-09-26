# Set the working directory to the directory with the CSV of results data
setwd("~/Documents/Programming_Projects/vaccinetrials")

# Load additional packages used in the algorithm
library(pwr)
library(tidyverse)

# Read CNEP data (real PWID survey data from Chicago) for use in demographic comparisons and as target population
cnep <- read.csv("data/cnep_plus_all_2018.02.13.csv")
cnep_susceptible <- cnep[cnep$HCV == "susceptible",]

# Stratify age categories in CNEP data for demographic matching
cnep_susceptible$Age_Category <- NA
cnep_susceptible[cnep_susceptible$Age < 20,]$Age_Category <- "<20"
cnep_susceptible[cnep_susceptible$Age >= 20 & cnep_susceptible$Age < 30,]$Age_Category <- "20-29"
cnep_susceptible[cnep_susceptible$Age >= 30 & cnep_susceptible$Age < 40,]$Age_Category <- "30-39"
cnep_susceptible[cnep_susceptible$Age >= 40 & cnep_susceptible$Age < 50,]$Age_Category <- "40-49"
cnep_susceptible[cnep_susceptible$Age >= 50,]$Age_Category <- "49+"

###########################################################################
## Find the demographic composition of the target population
###########################################################################
# Target pop set to susceptible PWID
target_pop <- cnep_susceptible
target_gender_comp <- matrix(table(target_pop$Gender)/nrow(target_pop), nrow = 1)
colnames(target_gender_comp) <- c("Female", "Male")
target_race_comp <- matrix(table(target_pop$Race)/nrow(target_pop), nrow = 1)
colnames(target_race_comp) <- c("Hispanic", "NHBlack", "NHWhite", "Other")
target_age_comp <- matrix(table(target_pop$Age_Category)/nrow(target_pop), nrow = 1)
colnames(target_age_comp) <- c("<20", "20-29", "30-39", "40-49", "49+")

# # Target demographic only includes gender and race
# target_demo <- cbind(target_gender_comp, target_race_comp)

# Target demographic includes gender, race, and age
target_demo <- cbind(target_gender_comp, target_race_comp, target_age_comp)


#########################################################################
########################### RANDOM ALGORITHM FUNCTION ###################
#########################################################################
# Random recruitment algorithm that goes through patients one at a time and randomly decides recruitment
# Parameters:
# recruitment_dataset - the recruitment pool from which PWID are sampled
# model_used - the survival analysis (Cox or RSF) model used to predict probability of HCV infection
# trial_followup_years - the follow up period for the simulated trial
# req_sample_size - the required sample size for the trial
# recruitment_probability - probability of recruiting each eligible patient (default 0.5)
# print_diagnostics - T/F value to print simulation details

random_algorithm <- function(recruitment_dataset, model_used, trial_followup_years, req_sample_size, recruitment_probability = 0.5, print_diagnostics = FALSE) {
  
  # Load prediction model functions
  source("prediction_models.R")
  
  # Record start time of the simulation
  start_time <- proc.time()
  
  # Initialize variables
  algorithm_output <- data.frame()
  agent_count <- 0
  
  # Randomly shuffle the recruitment dataset to simulate random patient encounters
  shuffled_dataset <- recruitment_dataset[sample(nrow(recruitment_dataset)),]
  
  # Go through patients one at a time
  for(i in 1:nrow(shuffled_dataset)) {
    # Break if we've reached the required sample size
    if(nrow(algorithm_output) >= req_sample_size) {
      break
    }
    
    # Get current patient
    current_patient <- shuffled_dataset[i, , drop = FALSE]
    agent_count <- agent_count + 1
    
    # Apply the predictive model to the current patient
    if(class(model_used)[1] == "coxph"){
      current_patient_predicted <- apply_cox(current_patient, model_used, trial_followup_years)
    }
    
    if(class(model_used)[1] == "rfsrc"){
      current_patient_predicted <- apply_rsf(current_patient, model_used, trial_followup_years)
    }
    
    # Check if patient is susceptible (eligible)
    if(current_patient_predicted$susceptible == 1) {
      # Randomly decide whether to recruit this patient
      random_decision <- runif(1)
      
      if(random_decision <= recruitment_probability) {
        # Recruit this patient
        algorithm_output <- rbind(algorithm_output, current_patient_predicted)
        
        if(print_diagnostics) {
          print(paste("Patient", i, "recruited. Total recruited:", nrow(algorithm_output)))
        }
      }
    }
    
    # Print progress periodically
    if(print_diagnostics && agent_count %% 100 == 0) {
      print(paste("Screened", agent_count, "patients, recruited", nrow(algorithm_output)))
    }
    
    # Safety check to prevent infinite loops
    if(agent_count >= nrow(recruitment_dataset)) {
      if(print_diagnostics) {
        print("Warning: Exhausted all patients in recruitment dataset before reaching target sample size")
      }
      break
    }
  }
  
  # Print elapsed time for algorithm completion
  if(print_diagnostics) {
    print(paste("Final results: Screened", agent_count, "patients, recruited", nrow(algorithm_output)))
  }
  print(proc.time()-start_time)
  
  # Construct list to be returned
  output_list <- list("algorithm_output" = algorithm_output, "agent_count" = agent_count)
  return(output_list)
}

### This random algorithm function represents a simple random recruitment approach where each eligible patient 
### has an equal probability of being recruited, regardless of their predicted infection risk or demographic characteristics.
### This serves as a baseline comparison to the more sophisticated PREDICTEE batch algorithm.