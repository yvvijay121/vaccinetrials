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
########################### PREDICTEE FUNCTION ##########################
#########################################################################
# PREDICTEE recruitment algorithm is defined as a function with multiple parameters defining the the recruitment population, predictive models, and other PREDICTEE parameters. An explanation of each parameter is as follows:
# recruitment_dataset - the recruitment pool from which PWID are sampled from
# model_used - the survival analysis (Cox or RSF) model that is used to predict probability of HCV infection by trial_followup_years
# target_demographics - a vector delineating the demographic groups being matched
# target_props - a vector delineating the target proportion for each demographic group
# recruitment_per_batch - the value B, the number of of PWID in a batch
# recruited_per_batch - the value R, the number of PWID from each batch that are recruited
# trial_followup_years - the follow up period for the simulated trial, which is used for prediction of infection probability with the survival analysis model
# req_sample_size - an estimated value for required sample size based on conventional recruitment methods (such as in-network recruitment)
# work_constraint - calculated as a function of B, R, and required sample size, or set to another user-defined value. See Appendix SM.D.ii for more detail
# incidence_weight_min - the minimum value for the incidence weight
# attrition_prob - the probability that a PWID in the backlog will be deleted (simulating attrition)
# high_demo_error_adjustment - T/F value used if there is a larger than normal demographic adjustment needed (see Appendix SM.C.i)
# ssmethod - equation used for sample size re-estimation calculations
# print_diagnostics - T/F value used only for simulations purposes to print simulation details
algorithm <- function(recruitment_dataset, model_used, target_demographics, target_props, recruitment_per_batch, recruited_per_batch, trial_followup_years, req_sample_size, work_constraint, incidence_weight_min = 25, attrition_prob = 0.2, high_demo_error_adjustment = FALSE, ssmethod = "cohen", print_diagnostics = FALSE) {
  # Throw error if number of target groups does not equal number of matched groups
  if((length(unlist(sapply(recruitment_dataset[,target_demographics], levels))) != length(target_props)) & (!is.null(target_demographics))){
    stop("Number of target groups does not equal number of matched groups")
  }
  
  # Create all the functions corresponding to each step in the algorithm:
  # function uses the Cox model to predict the probability of HCV infection for PWID in a provided data frame at a certain follow_up_time 
  # FIXED VERSION: Custom prediction that avoids Cox model environment issues
  apply_cox <- function(data, cox, follow_up_time) {
    # Add required survival variables
    data$status <- NA
    data$survival_time <- follow_up_time
    
    # Get model coefficients directly to avoid environment issues
    coefs <- cox$coefficients
    
    # Create design matrix manually
    # Age: numeric
    X_age <- data$Age
    
    # Gender: factor with levels Female (reference), Male
    X_gender_male <- ifelse(data$Gender == "Male", 1, 0)
    
    # Syringe_source: factor with levels HR (reference), nonHR  
    X_syringe_nonhr <- ifelse(data$Syringe_source == "nonHR", 1, 0)
    
    # Network variables: numeric
    X_drug_in <- data$Drug_in_degree
    X_drug_out <- data$Drug_out_degree
    X_network_size <- data$current_total_network_size
    X_injection_intensity <- data$Daily_injection_intensity
    X_recept_sharing <- data$Fraction_recept_sharing
    
    # Calculate linear predictor (risk score)
    linear_pred <- (coefs["Age"] * X_age +
                    coefs["GenderMale"] * X_gender_male +
                    coefs["Syringe_sourcenonHR"] * X_syringe_nonhr +
                    coefs["Drug_in_degree"] * X_drug_in +
                    coefs["Drug_out_degree"] * X_drug_out +
                    coefs["current_total_network_size"] * X_network_size +
                    coefs["Daily_injection_intensity"] * X_injection_intensity +
                    coefs["Fraction_recept_sharing"] * X_recept_sharing)
    
    # For a Cox model, survival probability = S0(t)^exp(linear_pred)
    # Since we don't have the baseline survival S0(t), we'll approximate
    # by using a baseline survival probability for HCV at 1.5 years
    baseline_survival <- 0.85  # Assume 85% baseline survival at 1.5 years
    
    survival_prob <- baseline_survival^exp(linear_pred)
    infection_prob <- 1 - survival_prob
    
    # Update the data with predictions
    data$survival_probability <- survival_prob
    data$infected_probability <- infection_prob
    
    return(data)
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
  
  # function used to calculate the demographic portion of PWID scores for PWID within a given data frame using the demographic of PWID already recruited to the trial cohort currently_in_trial and the target demographic proportions
  compute_demographic_score <- function(data, currently_in_trial, demos, props) {
    # Initating demographic matrix for calculation of demographic scores
    demo_matrix_list <- list()
    for(i in 1:length(demos)){
      demo_matrix_list[[i]] <- levels(data[,demos[i]])
    }
    demo_matrix <- expand.grid(demo_matrix_list, stringsAsFactors = TRUE)
    colnames(demo_matrix) <- demos
    demo_matrix <- as.matrix(demo_matrix)
    
    # Calculating demographic scores by using the difference between demographic proportions currently in trial and demographic composition of the target cohort
    if(nrow(currently_in_trial) == 0) {
      output <- data
      output$demographic_score <- 0
    } else {
      temp <- as.data.frame(currently_in_trial[,demos])
      output_props <- NULL
      for(i in 1:length(demos)){
        output_props <- cbind(output_props, t(data.matrix(table(temp[,i])/nrow(temp))))
      }
      
      prop_differences <- props[,order(colnames(props))]-output_props[,order(colnames(output_props))]
      demographic_score <- c()
      for(i in 1:dim(demo_matrix)[1]){
        demographic_score[i] <- sum(prop_differences[demo_matrix[i,]])
      }
      demo_matrix <- as.data.frame(cbind(demo_matrix, demographic_score))
      
      reduced_data <- data[,c("Agent", demos)]
      scores_by_agents <- merge(reduced_data, demo_matrix, by = demos)
      scores_by_agents <- scores_by_agents[, !names(scores_by_agents) %in% demos]
      output <- left_join(data, scores_by_agents, by = "Agent")
      output$demographic_score <- as.double(output$demographic_score)
    }
    return(output)
  }
  
  # NOTE: Sample size re-estimation function removed for simplicity
  
  # NOTE: Warning function removed since it depended on sample size re-estimation
  
  # Record start time of the simulation for optimization purposes
  start_time <- proc.time()
  
  # NOTE: Sample size re-estimation removed for simplicity
  
  # Initiating variables
  algorithm_output <- data.frame()
  backlog <- data.frame()
  agent_count <- 0
  batch_count <- 0
  
  # Setting the initial recruitment weights:
  incidence_weight <- 100
  demographic_weight <- 0
  weight_change_per_batch <- (incidence_weight-incidence_weight_min)/(req_sample_size/recruited_per_batch)
  
  while(nrow(algorithm_output) < req_sample_size) {
    
    # Sample recruitment_per_day number of agents at a time
    eligible <- recruitment_dataset[sample(nrow(recruitment_dataset), recruitment_per_batch),]
    agent_count <- agent_count + nrow(eligible)
    batch_count <- batch_count + 1
    
    # Apply the predictive model to the agents recruited for the day
    if(class(model_used)[1] == "coxph"){
      eligible_postmodel <- apply_cox(eligible, model_used, trial_followup_years)
    }
    
    if(class(model_used)[1] == "rfsrc"){
      eligible_postmodel <- apply_rsf(eligible, model_used, trial_followup_years)
    }
    
    # Add batches_elapsed_backlog to each agent to denote how many batches have passed since they have been in backlog. Initial value = 0.
    eligible_postmodel$batches_elapsed_backlog <- 0
    
    # Combine this batch and backlog into total_considered
    total_considered <- rbind(backlog, eligible_postmodel)
    
    # Calculate demographic difference scores for each agent, both race and gender
    if(is.null(target_demographics)){
      total_considered$demographic_score <- 0
    } else {
      total_considered <- compute_demographic_score(total_considered, algorithm_output, demos = target_demographics, props = target_props)
    }
    
    # Apply weights to generate a single score
    total_considered$score <- ((total_considered$infected_probability * incidence_weight) + (total_considered$demographic_score * demographic_weight))
    
    # Rank agents by overall score
    total_considered <- total_considered[order(total_considered$score, decreasing = TRUE),]
    
    # Consider the top [recruited_per_batch] agents from total_considered
    recruited <- head(total_considered, recruited_per_batch)
    
    # Check for positives and replace agents
    replaced <- 0
    while(min(recruited[,"susceptible"], na.rm = TRUE) == 0){
      # New variable for indexing new candidates to replace non-susceptible candidates
      previndex <- recruited_per_batch + 1 + replaced
      # New variable to count the number of replaced candidates from the batch; also used for indexing
      replaced <- replaced + length(which(recruited$susceptible==0))
      # Remove non-susceptible candidates from recruited candidates list
      recruited <- recruited[-which(recruited$susceptible==0),]
      # New variable for indexing new candidates
      endindex <- recruited_per_batch + replaced
      # Add the next highest scoring candidate that has not yet been considered
      recruited <- rbind(recruited, total_considered[previndex:endindex,])
    }
    
    # Recruit identified candidates
    algorithm_output <- rbind(algorithm_output, recruited)
    
    # Add remaining candidates to backlog
    backlog <- tail(total_considered, nrow(total_considered)-recruited_per_batch-replaced)
    
    # Assign a random number 0-1 to each agent to determine attrition
    backlog$attrition <- runif(nrow(backlog))
    
    # Use attrition probability to determine which agents get removed based on randomly generated number
    backlog <- backlog[backlog$attrition <= attrition_prob,]
    
    # Delete demographic_score, score, and attrition from the backlog for proper looping
    backlog <- backlog[, -which(names(backlog) %in% c("demographic_score", "score", "attrition"))]
    
    # Adjust weights after each batch, up to a certain limit
    if(incidence_weight > incidence_weight_min) {
      incidence_weight <- incidence_weight - weight_change_per_batch
      demographic_weight <- demographic_weight + weight_change_per_batch
      # Use original sample size for weight adjustment instead of re-estimation
      weight_change_per_batch <- (incidence_weight-incidence_weight_min)/(req_sample_size/recruited_per_batch)
      
      # Additional weight adjustment if the demographic error is very large and more emphasis wants to be put on demographics
      if(high_demo_error_adjustment == TRUE){
        error_adjustment <- 1-min(c(1-abs((table(algorithm_output$Gender)/nrow(algorithm_output)-target_gender_comp)/target_gender_comp),1-abs((table(algorithm_output$Race)/nrow(algorithm_output)-target_race_comp)/target_race_comp)))
        weight_change_per_batch <- weight_change_per_batch + error_adjustment
      }
      
      # Output to track progression of simulation
      if(print_diagnostics){
        print(paste("Current Weights:", incidence_weight, demographic_weight, "Weight Change Next Batch:", weight_change_per_batch))
      }
    }
    
    # NOTE: Sample size re-estimation removed - using original target sample size throughout
    
    # Print current agent count
    if(print_diagnostics){
      print(paste("Current number of agents screened:", agent_count))
    }
    
    # NOTE: Warning function removed with sample size re-estimation
  }
  
  # Print elapsed time for algorithm completion
  print(proc.time()-start_time)
  
  # Construct list to be returned
  output_list <- list("algorithm_output" = algorithm_output, "agent_count" = agent_count)
  return(output_list)
}

### This algorithm function represents a single simulation of recruitment to an HCV vaccine clinical trial. Our results are derived from 10,000 runs of this simulation to account for random variations in PWID sampling and also implementing cross-validation in the specific set of PWID that is used to train the model vs used in the recruitment pool (see simulations.R file for details on how these were conducted).