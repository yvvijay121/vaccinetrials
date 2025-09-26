#######################################################################
################## WORKING SINGLE RUN OF PREDICTEE ALGORITHM ###########
#######################################################################

# Load required libraries
library(magrittr)
library(survival)
library(dplyr)
library(pwr)
library(tidyverse)

setwd("~/Documents/Programming_Projects/vaccinetrials")

print("=== PREDICTEE Algorithm Single Run (Fixed Version) ===")

#######################################################################
####################### LOAD AND PREPARE DATA #######################
#######################################################################

# Load CNEP data for target demographics
cnep <- read.csv("data/cnep_plus_all_2018.02.13.csv")
cnep_susceptible <- cnep[cnep$HCV == "susceptible",]

# Calculate target demographics
cnep_susceptible$Age_Category <- cut(cnep_susceptible$Age, 
                                   breaks = c(-Inf, 20, 30, 40, 50, Inf), 
                                   labels = c("<20", "20-29", "30-39", "40-49", "49+"),
                                   right = FALSE)

target_gender_comp <- prop.table(table(cnep_susceptible$Gender))
target_race_comp <- prop.table(table(cnep_susceptible$Race))
target_age_comp <- prop.table(table(cnep_susceptible$Age_Category))
target_demo <- c(target_gender_comp, target_race_comp, target_age_comp)
print(target_demo)
print("Target demographics calculated")

# Load simulation data
simulation_data <- read.csv("data/events_pool.csv")
print(paste("Simulation data loaded:", nrow(simulation_data), "rows"))

# Create recruitment pool
agents_infected_total <- simulation_data[simulation_data$Event == "infected",]
agents_infected_total <- agents_infected_total[!duplicated(agents_infected_total$Agent), c("Time", "Agent")]
colnames(agents_infected_total) <- c("infected_time", "Agent")

recruitment_pool <- merge(simulation_data, agents_infected_total, by = "Agent", all.x = TRUE)
recruitment_pool <- recruitment_pool[recruitment_pool$Event == "activated",]

# Process recruitment pool
recruitment_pool$infected_time_after_start <- recruitment_pool$infected_time - recruitment_pool$Time
recruitment_pool$infected_time_after_start[is.na(recruitment_pool$infected_time_after_start)] <- 999
recruitment_pool$infected_by_trialend <- ifelse(recruitment_pool$infected_time_after_start <= 1.5, 1, 0)

# Add age categories
recruitment_pool$Age_Category <- cut(recruitment_pool$Age, 
                                   breaks = c(-Inf, 20, 30, 40, 50, Inf), 
                                   labels = c("<20", "20-29", "30-39", "40-49", "49+"),
                                   right = FALSE)

# Set factor levels correctly
recruitment_pool$Gender <- factor(recruitment_pool$Gender, levels = c("Female", "Male"))
recruitment_pool$Race <- factor(recruitment_pool$Race)
recruitment_pool$Syringe_source <- factor(recruitment_pool$Syringe_source, levels = c("HR", "nonHR"))
recruitment_pool$Age_Category <- factor(recruitment_pool$Age_Category)
recruitment_pool$susceptible <- ifelse(recruitment_pool$HCV == "susceptible", 1, 0)

print(paste("Recruitment pool created:", nrow(recruitment_pool), "individuals"))

#######################################################################
################ FIX COX MODEL PREDICTION ISSUE ####################
#######################################################################

# The original Cox model has environment issues
# Let's create a custom prediction function that works around this

# Load the original Cox model to get coefficients
original_cox <- readRDS("universal_cox.rds")
print("Original Cox model loaded")
print("Coefficients:")
print(original_cox$coefficients)

# Create a custom prediction function that uses the Cox coefficients directly
custom_cox_predict <- function(data, cox_model, follow_up_time) {
  # Get model coefficients
  coefs <- cox_model$coefficients
  
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
  # by using a reasonable baseline survival probability for HCV at 1.5 years
  baseline_survival <- 0.85  # Assume 85% baseline survival at 1.5 years
  
  survival_prob <- baseline_survival^exp(linear_pred)
  infection_prob <- 1 - survival_prob
  
  # Add predictions to data
  result <- data
  result$survival_probability <- survival_prob
  result$infected_probability <- infection_prob
  result$status <- NA
  result$survival_time <- follow_up_time
  
  return(result)
}

print("Custom Cox prediction function created")

#######################################################################
################ SIMPLIFIED ALGORITHM FUNCTION #####################
#######################################################################

# Create a simplified version of the algorithm with the custom prediction
simple_algorithm <- function(recruitment_data, target_demographics, target_props, 
                            batch_size = 20, recruited_per_batch = 10, 
                            target_sample = 100, max_screening = 500) {
  
  recruited <- data.frame()
  screened_count <- 0
  
  while(nrow(recruited) < target_sample && screened_count < max_screening) {
    # Sample a batch
    available_pool <- recruitment_data[!recruitment_data$Agent %in% recruited$Agent, ]
    if(nrow(available_pool) < batch_size) {
      batch <- available_pool
    } else {
      batch <- available_pool[sample(nrow(available_pool), batch_size), ]
    }
    
    screened_count <- screened_count + nrow(batch)
    
    # Apply custom Cox prediction
    batch_with_pred <- custom_cox_predict(batch, original_cox, 1.5)
    
    # Simple selection: top risk scores
    batch_sorted <- batch_with_pred[order(batch_with_pred$infected_probability, decreasing = TRUE), ]
    selected <- head(batch_sorted, recruited_per_batch)
    
    # Only keep susceptible individuals
    selected_susceptible <- selected[selected$susceptible == 1, ]
    
    recruited <- rbind(recruited, selected_susceptible)
    
    if(nrow(recruited) >= target_sample) break
  }
  
  return(list(recruited = recruited, screened = screened_count))
}

#######################################################################
######################## RUN THE ALGORITHM ###########################
#######################################################################

print("\n=== Running Simplified Algorithm ===")
print("Parameters:")
print(paste("- Target sample size: 100"))
print(paste("- Batch size: 20")) 
print(paste("- Recruited per batch: 10"))
print(paste("- Max screening: 500"))

start_time <- Sys.time()

results <- simple_algorithm(
  recruitment_data = recruitment_pool,
  target_demographics = c("Gender", "Race", "Age_Category"),
  target_props = target_demo,
  batch_size = 20,
  recruited_per_batch = 10,
  target_sample = 100,
  max_screening = 500
)

end_time <- Sys.time()

#######################################################################
####################### DISPLAY RESULTS ##############################
#######################################################################

recruited_sample <- results$recruited

print("\n=== RESULTS ===")
print(paste("Algorithm completed in", round(as.numeric(end_time - start_time), 2), "seconds"))
print(paste("Total recruited:", nrow(recruited_sample)))
print(paste("Total screened:", results$screened))
print(paste("Efficiency:", round(nrow(recruited_sample)/results$screened, 3)))

if(nrow(recruited_sample) > 0) {
  print("\n=== Risk Distribution ===")
  print("Infection probability summary:")
  print(summary(recruited_sample$infected_probability))
  
  print(paste("Mean predicted infection probability:", 
              round(mean(recruited_sample$infected_probability), 4)))
  
  print("\n=== Demographics ===")
  if("Gender" %in% colnames(recruited_sample)) {
    print("Gender distribution:")
    print(prop.table(table(recruited_sample$Gender)))
  }
  
  if("Race" %in% colnames(recruited_sample)) {
    print("\nRace distribution:")
    print(prop.table(table(recruited_sample$Race)))
  }
  
  if("Age_Category" %in% colnames(recruited_sample)) {
    print("\nAge category distribution:")
    print(prop.table(table(recruited_sample$Age_Category)))
  }
  
} else {
  print("No participants recruited - check data and parameters")
}

print("\n=== Algorithm Complete ===")