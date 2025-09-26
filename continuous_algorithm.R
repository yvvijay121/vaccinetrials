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
########################### CONTINUOUS PREDICTEE FUNCTION ##############
#########################################################################
# CONTINUOUS PREDICTEE recruitment algorithm is defined as a function with multiple parameters defining the recruitment population, predictive models, and other PREDICTEE parameters. An explanation of each parameter is as follows:
# recruitment_dataset - the recruitment pool from which PWID are sampled from
# model_used - the survival analysis (Cox or RSF) model that is used to predict probability of HCV infection by trial_followup_years
# target_demographics - a vector delineating the demographic groups being matched
# target_props - a vector delineating the target proportion for each demographic group
# trial_followup_years - the follow up period for the simulated trial, which is used for prediction of infection probability with the survival analysis model
# req_sample_size - an estimated value for required sample size based on conventional recruitment methods (such as in-network recruitment)
# high_risk_percentile - percentile threshold for high-risk candidates (default: top 5%)
# medium_risk_percentile - percentile threshold for medium-risk candidates (default: top 15%) 
# low_risk_percentile - percentile threshold for low-risk candidates (default: top 50%)
# min_clinical_threshold - minimum clinically meaningful risk threshold (default: 1% annual risk)
# print_diagnostics - T/F value used only for simulations purposes to print simulation details
continuous_algorithm <- function(recruitment_dataset, model_used, target_demographics, target_props, trial_followup_years, req_sample_size, high_risk_percentile = 0.05, medium_risk_percentile = 0.15, low_risk_percentile = 0.5, min_clinical_threshold = 0.01, print_diagnostics = FALSE) {

  # Throw error if number of target groups does not equal number of matched groups
  if((length(unlist(sapply(recruitment_dataset[,target_demographics], levels))) != length(target_props)) & (!is.null(target_demographics))){
    stop("Number of target groups does not equal number of matched groups")
  }
  
  # Create all the functions corresponding to each step in the algorithm:
  # function uses the Cox model to predict the probability of HCV infection for PWID in a provided data frame at a certain follow_up_time 
  apply_cox <- function(data, cox, follow_up_time) {
    data$status <- NA
    data$survival_time <- follow_up_time
    probabilities <- as.data.frame(predict(cox, data, type = "survival"))
    colnames(probabilities) <- "survival_probability"
    output <- cbind(data, probabilities)
    output$infected_probability <- 1-output$survival_probability
    
    return(output)
  }
  
  # function uses the RSF model to predict the probability of HCV infection for PWID in a provided data frame at a certain follow_up_time 
  apply_rsf <- function(data, model, follow_up_time) {
    prediction_data <- predict.rfsrc(model, data, na.action = "na.impute")
    time_interest <- which(abs(prediction_data$time.interest-follow_up_time)==min(abs(prediction_data$time.interest-follow_up_time)))
    probabilities <- as.data.frame(prediction_data$survival[,time_interest])
    colnames(probabilities) <- "survival_probability"
    output <- cbind(data, probabilities)
    output$infected_probability <- 1-output$survival_probability
    
    return(output)
  }
  
  #########################################################################
  ########################### DEMOGRAPHIC COUNTER SYSTEM #################
  #########################################################################
  
  # Initialize demographic counters if target_demographics is provided
  demographic_counters <- NULL
  target_proportions <- NULL
  
  if(!is.null(target_demographics) && !is.null(target_props)) {
    # Create a list to store counters for each demographic category
    demographic_counters <- list()
    target_proportions <- list()
    
    # Extract target proportions for each demographic category
    current_col <- 1
    for(demo_var in target_demographics) {
      # Get the levels for this demographic variable
      demo_levels <- levels(recruitment_dataset[,demo_var])
      num_levels <- length(demo_levels)
      
      # Extract the target proportions for this demographic variable
      target_props_for_var <- target_props[current_col:(current_col + num_levels - 1)]
      names(target_props_for_var) <- demo_levels
      
      # Initialize counters to zero for each level
      current_counts <- rep(0, num_levels)
      names(current_counts) <- demo_levels
      
      # Store in lists
      demographic_counters[[demo_var]] <- current_counts
      target_proportions[[demo_var]] <- target_props_for_var
      
      # Update column index
      current_col <- current_col + num_levels
    }
  }
  
  # Function to update demographic counters when a candidate is recruited
  update_demographic_counters <- function(recruited_candidate, counters) {
    if(!is.null(counters) && !is.null(target_demographics)) {
      for(demo_var in target_demographics) {
        candidate_value <- as.character(recruited_candidate[,demo_var])
        counters[[demo_var]][candidate_value] <- counters[[demo_var]][candidate_value] + 1
      }
    }
    return(counters)
  }
  
  # Function to calculate demographic priority scores
  calculate_demographic_priority <- function(candidate, counters, target_props, total_recruited) {
    if(is.null(counters) || is.null(target_props) || total_recruited == 0) {
      return(0)  # No demographic priority if no targets set or no one recruited yet
    }
    
    priority_score <- 0
    
    for(demo_var in target_demographics) {
      candidate_value <- as.character(candidate[,demo_var])
      
      # Calculate current proportion for this demographic value
      current_count <- counters[[demo_var]][candidate_value]
      current_prop <- current_count / total_recruited
      
      # Get target proportion
      target_prop <- target_props[[demo_var]][candidate_value]
      
      # Calculate difference (positive means we need more of this demographic)
      prop_difference <- target_prop - current_prop
      
      # Add to priority score (weight by the magnitude of the difference)
      priority_score <- priority_score + prop_difference
    }
    
    return(priority_score)
  }
    
  #########################################################################
  ########################### ADAPTIVE THRESHOLD CALCULATION #############
  #########################################################################
  
  # Calculate adaptive thresholds based on recruitment goals and clinical constraints
  # Apply model to entire recruitment dataset to understand probability distribution  
  recruitment_with_probs <- recruitment_dataset
  if(class(model_used)[1] == "coxph"){
    recruitment_with_probs <- apply_cox(recruitment_with_probs, model_used, trial_followup_years)
  } else if(class(model_used)[1] == "rfsrc"){
    recruitment_with_probs <- apply_rsf(recruitment_with_probs, model_used, trial_followup_years)
  }
  
  # Calculate percentile-based thresholds from the current recruitment pool
  prob_quantiles <- quantile(recruitment_with_probs$infected_probability, 
                            probs = c(1 - low_risk_percentile, 1 - medium_risk_percentile, 1 - high_risk_percentile), 
                            na.rm = TRUE)
  
  # Set adaptive thresholds with clinical floor constraints
  high_risk_threshold <- max(prob_quantiles[3], min_clinical_threshold * 3)    # At least 3% or top 5%
  medium_risk_threshold <- max(prob_quantiles[2], min_clinical_threshold * 2)  # At least 2% or top 15%  
  low_risk_threshold <- max(prob_quantiles[1], min_clinical_threshold)         # At least 1% or top 50%
  
  if(print_diagnostics) {
    cat("=== ADAPTIVE THRESHOLDS CALCULATED ===\n")
    cat(sprintf("High-risk threshold (top %0.0f%%): %0.4f\n", high_risk_percentile*100, high_risk_threshold))
    cat(sprintf("Medium-risk threshold (top %0.0f%%): %0.4f\n", medium_risk_percentile*100, medium_risk_threshold))
    cat(sprintf("Low-risk threshold (top %0.0f%%): %0.4f\n", low_risk_percentile*100, low_risk_threshold))
    cat(sprintf("Clinical floor: %0.4f\n", min_clinical_threshold))
    cat("=========================================\n")
  }
  
  #########################################################################
  ########################### ALGORITHM VARIABLES ########################
  #########################################################################
  
  # Initiating variables
  algorithm_output <- data.frame()
  candidate_pool <- data.frame()
  agent_count <- 0
  time_step <- 0
  
  #########################################################################
  ########################### MAIN ALGORITHM LOOP ########################
  #########################################################################
  
  while(nrow(algorithm_output) < req_sample_size) {
    time_step <- time_step + 1
    
    # Sample a candidate from the recruitment pool
    if(nrow(recruitment_dataset) > 0) {
      candidate <- recruitment_dataset[sample(nrow(recruitment_dataset), 1), ]
      agent_count <- agent_count + 1
      
      # Apply predictive model to candidate
      if(class(model_used)[1] == "coxph"){
        candidate <- apply_cox(candidate, model_used, trial_followup_years)
      } else if(class(model_used)[1] == "rfsrc"){
        candidate <- apply_rsf(candidate, model_used, trial_followup_years)
      }
      
      # Calculate demographic priority score
      total_recruited <- nrow(algorithm_output)
      demographic_priority <- calculate_demographic_priority(candidate, demographic_counters, target_proportions, total_recruited)
      
      # Enhanced recruitment decision based on adaptive infected probability criteria
      # Recruit if candidate's infected_probability is above high_risk_threshold, or
      # If demographic priority is positive and infected_probability is above medium_risk_threshold, or
      # Random chance (10%) if infected_probability is above low_risk_threshold

      recruit_candidate <- FALSE

      if(candidate$infected_probability >= high_risk_threshold) {
        recruit_candidate <- TRUE
      } else if(candidate$infected_probability >= medium_risk_threshold && demographic_priority > 0) {
        recruit_candidate <- TRUE
      } else if(candidate$infected_probability >= low_risk_threshold && runif(1) < 0.1) {
        recruit_candidate <- TRUE
      }
      
      # Recruit the candidate if selected
      if(recruit_candidate) {
        algorithm_output <- rbind(algorithm_output, candidate)
        
        # Update demographic counters
        demographic_counters <- update_demographic_counters(candidate, demographic_counters)
        
        # Print diagnostic information
        if(print_diagnostics && nrow(algorithm_output) %% 50 == 0) {
          cat(paste("Recruited", nrow(algorithm_output), "candidates after", time_step, "time steps\n"))
        }
      }
    } else {
      # No more candidates in recruitment dataset
      cat("WARNING: Recruitment dataset exhausted\n")
      break
    }
    
    # Safety break to prevent infinite loops
    if(time_step > 50000) {
      cat("WARNING: Maximum time steps reached. Breaking loop.\n")
      break
    }
  }
  
  # Construct list to be returned
  output_list <- list("algorithm_output" = algorithm_output, "agent_count" = agent_count, "time_steps" = time_step)
  return(output_list)
}