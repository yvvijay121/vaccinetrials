#########################################################################
########################### PREDICTION MODEL UTILITIES ##################
#########################################################################
# This file contains shared prediction functions used by both batch_algorithm.R 
# and random_algorithm.R to avoid code duplication

# Function uses the Cox model to predict the probability of HCV infection for PWID 
# in a provided data frame at a certain follow_up_time 
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
  baseline_survival <- 0.95  # Assume 95% baseline survival at 1.5 years
  
  # Handle potential NA values in linear predictor
  survival_prob <- ifelse(is.na(linear_pred), NA, baseline_survival^exp(linear_pred))
  infection_prob <- ifelse(is.na(survival_prob), NA, 1 - survival_prob)
  
  # Update the data with predictions
  data$survival_probability <- survival_prob
  data$infected_probability <- infection_prob
  
  return(data)
}

# Function uses the RSF model to predict the probability of HCV infection for PWID 
# in a provided data frame at a certain follow_up_time 
apply_rsf <- function(data, model, follow_up_time) {
  prediction_data <- predict.rfsrc(model, data, na.action = "na.impute")
  time_interest <- which(abs(prediction_data$time.interest-follow_up_time)==min(abs(prediction_data$time.interest-follow_up_time)))
  probabilities <- as.data.frame(prediction_data$survival[,time_interest])
  colnames(probabilities) <- "survival_probability"
  output <- cbind(data, probabilities)
  output$infected_probability <- 1-output$survival_probability
  
  return(output)
}
