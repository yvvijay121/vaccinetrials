#######################################################################
############################### ANALYSIS ##############################
############################### COX MODEL #############################
############ 100 RANDOM DBLABEL ASSIGNMENTS, 100 RUNS EACH ############
######################### WITH AGE DEMOGRAPHIC ########################
library(magrittr)
library(survival)
library(survminer)
library(dplyr)
library(pwr)
library(tidyverse)

# Set the working directory to the directory with the CSV of data
setwd("~/Documents/Programming_Projects/vaccinetrials")

source("batch_algorithm.R")
source("continuous_algorithm.R")

# Read the CSV file as a data frame and create the training set using 20% of the data.
simulation_data <- read.csv("data/events_pool.csv")

# Generate new DBLabel split
all_DBLabel <- 0:12013
DBLabel_train <- sample(all_DBLabel, size=2403)
DBLabel_test <- all_DBLabel[-(DBLabel_train+1)]

# Train the Cox model
# Train Cox model based on random DBLabel assignment
train <- simulation_data[simulation_data$DBLabel %in% DBLabel_train,]

# Create data set of unique agents when they are first activated
agents_activated <- train[train$Event == "activated",]

# Create dataset specifically of activated agents who are susceptible
agents_activated_susceptible <- agents_activated[agents_activated$HCV == "susceptible",]

# Create dataset of all the events for the susceptible agents
agents_susceptible <- train[train$Agent %in% agents_activated_susceptible$Agent,]

# Identify when agents were deactivated
agents_deactivated_time <- agents_susceptible[agents_susceptible$Event == "deactivated",]

# Identify when agents were infected for the first time
agents_infected_time <- agents_susceptible[agents_susceptible$Event == "infected",]
agents_infected_time <- agents_infected_time[!(duplicated(agents_infected_time$Agent)),]

# Clean agents_susceptible to remove unnecessary columns and duplicate records, and convert reclassify factor variables as factor
unique_agents <- agents_susceptible[agents_susceptible$Event == "activated", c("Time", "Agent", "Age", "Gender", "Race", "Syringe_source", "Age_Started", "Drug_in_degree", "Drug_out_degree", "current_total_network_size", "Daily_injection_intensity", "Fraction_recept_sharing", "DBLabel", "chicago_community_name")]
unique_agents$Gender <- factor(unique_agents$Gender)
unique_agents$Race <- factor(unique_agents$Race)
unique_agents$Syringe_source <- factor(unique_agents$Syringe_source)
unique_agents$chicago_community_name <- factor(unique_agents$chicago_community_name)


# Label censored (0) if never infected, or status (1) if infected
unique_agents$status <- 0
unique_agents$status[unique_agents$Agent %in% agents_infected_time$Agent] <- 1

# Record the time of the relevant event (deactivated or infected) for each agent
agents_deactivated_time <- agents_deactivated_time[, c("Time","Agent")]
colnames(agents_deactivated_time) <- c("deactivated_time", "Agent")

agents_infected_time <- agents_infected_time[, c("Time","Agent")]
colnames(agents_infected_time) <- c("infected_time", "Agent")

unique_agents <- merge(unique_agents, agents_deactivated_time, by = "Agent", all = TRUE)
unique_agents <- merge(unique_agents, agents_infected_time, by = "Agent", all = TRUE)

# Select the correct time for the "event_time"
# Hierarchy: infected > deactivated
unique_agents$event_time <- NA
unique_agents[is.na(unique_agents$event_time) == TRUE,]$event_time <- unique_agents[is.na(unique_agents$event_time) == TRUE,]$infected_time
unique_agents[is.na(unique_agents$event_time) == TRUE,]$event_time <- unique_agents[is.na(unique_agents$event_time) == TRUE,]$deactivated_time

# Some agents never reach any of these outcomes (are activated but nothing else happens). Thus they remain susceptible through the entire simulation, and their event time is assigned as 9.858 (last tick of simulation)

unique_agents[is.na(unique_agents$event_time) == TRUE,]$event_time <- 9.858
unique_agents$survival_time <- unique_agents$event_time - unique_agents$Time

# Create Cox Model
cox_model <- coxph(Surv(time = survival_time, event = status) ~ Age + Gender + Syringe_source + Drug_in_degree + Drug_out_degree + current_total_network_size + Daily_injection_intensity + Fraction_recept_sharing + chicago_community_name, data = unique_agents, x = TRUE)

# Create the recruitment pool

# Create the same test data as seen in the Cox Model validation
test <- simulation_data[simulation_data$DBLabel %in% DBLabel_test,]

# Label test data with the time each agent has their first infection, if applicable.
agents_infected_total <- simulation_data[simulation_data$Event == "infected",]
agents_infected_total <- agents_infected_total[!duplicated(agents_infected_total$Agent), c("Time", "Agent")]
colnames(agents_infected_total) <- c("infected_time", "Agent")

# Label test data with the time each agent becomes chronic for potential future analysis.
agents_chronic_total <- simulation_data[simulation_data$Event == "chronic",]
agents_chronic_total <- agents_chronic_total[!duplicated(agents_chronic_total$Agent), c("Time", "Agent")]
colnames(agents_chronic_total) <- c("chronic_time", "Agent")

# Create recruitment_pool data frame with time of first infection and time of chronicity
recruitment_pool <- merge(test, agents_infected_total, by = "Agent", all.x = TRUE)
recruitment_pool <- merge(recruitment_pool, agents_chronic_total, by = "Agent", all.x = TRUE)

# recruitment_pool$status <- NA
# recruitment_pool$survival_time <- 1.5
# Set recruitment pool only to unique agents when they are first activated
recruitment_pool <- recruitment_pool[recruitment_pool$Event == "activated",]

# Label if the agent was infected within 1.5 years of activation
recruitment_pool$infected_time_after_start <- recruitment_pool$infected_time - recruitment_pool$Time
recruitment_pool[is.na(recruitment_pool$infected_time_after_start),]$infected_time_after_start <- 999
recruitment_pool$infected_by_trialend <- 0
recruitment_pool[recruitment_pool$infected_time_after_start <= 1.5, ]$infected_by_trialend <- 1

# Label if agent becomes chronic within 1.5 years of activation (for potential future analysis)
recruitment_pool$chronic_time_after_start <- recruitment_pool$chronic_time - recruitment_pool$Time
recruitment_pool[is.na(recruitment_pool$chronic_time_after_start),]$chronic_time_after_start <- 999
recruitment_pool$chronic_by_trialend <- 0
recruitment_pool[recruitment_pool$chronic_time_after_start <= 1.5, ]$chronic_by_trialend <- 1

# Stratify age categories for demographic matching
recruitment_pool$Age_Category <- NA
recruitment_pool[recruitment_pool$Age < 20,]$Age_Category <- "<20"
recruitment_pool[recruitment_pool$Age >= 20 & recruitment_pool$Age < 30,]$Age_Category <- "20-29"
recruitment_pool[recruitment_pool$Age >= 30 & recruitment_pool$Age < 40,]$Age_Category <- "30-39"
recruitment_pool[recruitment_pool$Age >= 40 & recruitment_pool$Age < 50,]$Age_Category <- "40-49"
recruitment_pool[recruitment_pool$Age >= 50,]$Age_Category <- "49+"


# Reclassify recruitment pool variables as factors if applicable
recruitment_pool$Gender <- factor(recruitment_pool$Gender)
recruitment_pool$Race <- factor(recruitment_pool$Race)
recruitment_pool$Syringe_source <- factor(recruitment_pool$Syringe_source)
recruitment_pool$chicago_community_name <- factor(recruitment_pool$chicago_community_name)
recruitment_pool$Age_Category <- factor(recruitment_pool$Age_Category)

# Label if each agent is susceptible or not for use in the algorithm. Realistic algorithm will not know if they are susceptible or not until they undergo a "test," after which they will be removed and replaced with another agent who is not susceptible if they are determined to not be susceptible.
recruitment_pool$susceptible <- 0
recruitment_pool$susceptible[recruitment_pool$HCV == "susceptible"] <- 1

# Initiation of output data frames:
# Batch algorithm results
batch_expectation_vs_reality <- data.frame()
batch_race_analysis <- data.frame()
batch_gender_analysis <- data.frame()
batch_age_analysis <- data.frame()
batch_demographic_incidence <- data.frame()

# Continuous algorithm results
continuous_expectation_vs_reality <- data.frame()
continuous_race_analysis <- data.frame()
continuous_gender_analysis <- data.frame()
continuous_age_analysis <- data.frame()
continuous_demographic_incidence <- data.frame()

number_of_runs <- 1

while(nrow(batch_expectation_vs_reality) < number_of_runs) {
  
  print(paste("Starting run:", nrow(batch_expectation_vs_reality) + 1))
  
  #########################################################################
  ######################### BATCH ALGORITHM ##############################
  #########################################################################
  # Run the batch algorithm and store raw results in batch_algorithm_output
  batch_algorithm_list <- algorithm(recruitment_dataset = recruitment_pool, model_used = cox_model, target_demographics = c("Gender", "Race", "Age_Category"), target_props = target_demo, recruitment_per_batch = 50, recruited_per_batch = 5, trial_followup_years = 1.5, req_sample_size = 800, work_constraint = 8000, attrition_prob = 0.2)
  
  batch_algorithm_output <- batch_algorithm_list$algorithm_output
  
  #########################################################################
  ######################### CONTINUOUS ALGORITHM #########################
  #########################################################################
  # Run the continuous algorithm with equivalent parameters
  # Convert batch parameters to continuous equivalents:
  # screening_rate = recruitment_per_batch = 50 per time step
  # recruitment_rate = recruited_per_batch = 5 per time step
  # Use adaptive thresholds with clinically meaningful constraints:
  # - High risk: top 5% (or minimum 3% annual risk)  
  # - Medium risk: top 15% (or minimum 2% annual risk)
  # - Low risk: top 50% (or minimum 1% annual risk)
  continuous_algorithm_list <- continuous_algorithm(recruitment_dataset = recruitment_pool, model_used = cox_model, target_demographics = c("Gender", "Race", "Age_Category"), target_props = target_demo, trial_followup_years = 1.5, req_sample_size = 800, high_risk_percentile = 0.05, medium_risk_percentile = 0.15, low_risk_percentile = 0.5, min_clinical_threshold = 0.01)
  
  continuous_algorithm_output <- continuous_algorithm_list$algorithm_output
  
  #########################################################################
  ######################### BATCH ALGORITHM RESULTS ######################
  #########################################################################
  # Calculate incidence results and agent count for batch algorithm
  batch_newrow <- matrix(c(table(batch_algorithm_output$infected_by_trialend)[2]/sum(table(batch_algorithm_output$infected_by_trialend)), mean(batch_algorithm_output$infected_probability), table(batch_algorithm_output$chronic_by_trialend)[2]/sum(table(batch_algorithm_output$chronic_by_trialend)), batch_algorithm_list$agent_count), nrow = 1)
  colnames(batch_newrow) <- c("actual", "expected", "chronic", "agent_count")
  batch_expectation_vs_reality <- rbind(batch_expectation_vs_reality, batch_newrow)
  
  # Calculate demographic results for batch algorithm
  batch_racerow <- matrix(table(batch_algorithm_output$Race)/nrow(batch_algorithm_output), nrow = 1)
  colnames(batch_racerow) <- c("Hispanic", "NHBlack", "NHWhite", "Other")
  batch_genderrow <- matrix(table(batch_algorithm_output$Gender)/nrow(batch_algorithm_output), nrow = 1)
  colnames(batch_genderrow) <- c("Female", "Male")
  batch_agerow <- matrix(table(batch_algorithm_output$Age_Category)/nrow(batch_algorithm_output), nrow = 1)
  colnames(batch_agerow) <- c("<20", "20-29", "30-39", "40-49", "49+")
  
  batch_race_analysis <- rbind(batch_race_analysis, batch_racerow)
  batch_gender_analysis <- rbind(batch_gender_analysis, batch_genderrow)
  batch_age_analysis <- rbind(batch_age_analysis, batch_agerow)
  
  #########################################################################
  ######################### CONTINUOUS ALGORITHM RESULTS #################
  #########################################################################
  # Calculate incidence results and agent count for continuous algorithm
  if(nrow(continuous_algorithm_output) > 0) {
    continuous_newrow <- matrix(c(table(continuous_algorithm_output$infected_by_trialend)[2]/sum(table(continuous_algorithm_output$infected_by_trialend)), mean(continuous_algorithm_output$infected_probability), table(continuous_algorithm_output$chronic_by_trialend)[2]/sum(table(continuous_algorithm_output$chronic_by_trialend)), continuous_algorithm_list$agent_count), nrow = 1)
    colnames(continuous_newrow) <- c("actual", "expected", "chronic", "agent_count")
    continuous_expectation_vs_reality <- rbind(continuous_expectation_vs_reality, continuous_newrow)
    
    # Calculate demographic results for continuous algorithm
    continuous_racerow <- matrix(table(continuous_algorithm_output$Race)/nrow(continuous_algorithm_output), nrow = 1)
    colnames(continuous_racerow) <- c("Hispanic", "NHBlack", "NHWhite", "Other")
    continuous_genderrow <- matrix(table(continuous_algorithm_output$Gender)/nrow(continuous_algorithm_output), nrow = 1)
    colnames(continuous_genderrow) <- c("Female", "Male")
    continuous_agerow <- matrix(table(continuous_algorithm_output$Age_Category)/nrow(continuous_algorithm_output), nrow = 1)
    colnames(continuous_agerow) <- c("<20", "20-29", "30-39", "40-49", "49+")
    
    continuous_race_analysis <- rbind(continuous_race_analysis, continuous_racerow)
    continuous_gender_analysis <- rbind(continuous_gender_analysis, continuous_genderrow)
    continuous_age_analysis <- rbind(continuous_age_analysis, continuous_agerow)
  } else {
    print("WARNING: Continuous algorithm produced no results for this run")
  }
  
  #########################################################################
  ######################### DEMOGRAPHIC INCIDENCE ANALYSIS ################
  #########################################################################
  # Calculate incidence results within each demographic group for batch algorithm
  batch_race_incidence_actual <- matrix(tapply(batch_algorithm_output$infected_by_trialend, batch_algorithm_output$Race, sum)/table(batch_algorithm_output$Race), nrow = 1)
  colnames(batch_race_incidence_actual) <- c("hispanic_actual", "nhblack_actual", "nhwhite_actual", "other_actual")
  batch_race_incidence_expected <- matrix(tapply(batch_algorithm_output$infected_probability, batch_algorithm_output$Race, mean), nrow = 1)
  colnames(batch_race_incidence_expected) <- c("hispanic_expected", "nhblack_expected", "nhwhite_expected", "other_expected")
  batch_race_incidence <- cbind(batch_race_incidence_actual, batch_race_incidence_expected)
  
  batch_gender_incidence_actual <- matrix(tapply(batch_algorithm_output$infected_by_trialend, batch_algorithm_output$Gender, sum)/table(batch_algorithm_output$Gender), nrow = 1)
  colnames(batch_gender_incidence_actual) <- c("female_actual", "male_actual")
  batch_gender_incidence_expected <- matrix(tapply(batch_algorithm_output$infected_probability, batch_algorithm_output$Gender, mean), nrow = 1)
  colnames(batch_gender_incidence_expected) <- c("female_expected", "male_expected")
  batch_gender_incidence <- cbind(batch_gender_incidence_actual, batch_gender_incidence_expected)
  
  batch_age_incidence_actual <- matrix(tapply(batch_algorithm_output$infected_by_trialend, batch_algorithm_output$Age_Category, sum)/table(batch_algorithm_output$Age_Category), nrow = 1)
  colnames(batch_age_incidence_actual) <- c("less20_actual", "20to29_actual", "30to39_actual", "40to49_actual", "greater49_actual")
  batch_age_incidence_expected <- matrix(tapply(batch_algorithm_output$infected_probability, batch_algorithm_output$Age_Category, mean), nrow = 1)
  colnames(batch_age_incidence_expected) <- c("less20_expected", "20to29_expected", "30to39_expected", "40to49_expected", "greater49_expected")
  batch_age_incidence <- cbind(batch_age_incidence_actual, batch_age_incidence_expected)
  
  batch_demographic_incidence_row <- cbind(batch_race_incidence, batch_gender_incidence, batch_age_incidence)
  batch_demographic_incidence <- rbind(batch_demographic_incidence, batch_demographic_incidence_row)
  
  # Calculate incidence results within each demographic group for continuous algorithm
  if(nrow(continuous_algorithm_output) > 0) {
    continuous_race_incidence_actual <- matrix(tapply(continuous_algorithm_output$infected_by_trialend, continuous_algorithm_output$Race, sum)/table(continuous_algorithm_output$Race), nrow = 1)
    colnames(continuous_race_incidence_actual) <- c("hispanic_actual", "nhblack_actual", "nhwhite_actual", "other_actual")
    continuous_race_incidence_expected <- matrix(tapply(continuous_algorithm_output$infected_probability, continuous_algorithm_output$Race, mean), nrow = 1)
    colnames(continuous_race_incidence_expected) <- c("hispanic_expected", "nhblack_expected", "nhwhite_expected", "other_expected")
    continuous_race_incidence <- cbind(continuous_race_incidence_actual, continuous_race_incidence_expected)
    
    continuous_gender_incidence_actual <- matrix(tapply(continuous_algorithm_output$infected_by_trialend, continuous_algorithm_output$Gender, sum)/table(continuous_algorithm_output$Gender), nrow = 1)
    colnames(continuous_gender_incidence_actual) <- c("female_actual", "male_actual")
    continuous_gender_incidence_expected <- matrix(tapply(continuous_algorithm_output$infected_probability, continuous_algorithm_output$Gender, mean), nrow = 1)
    colnames(continuous_gender_incidence_expected) <- c("female_expected", "male_expected")
    continuous_gender_incidence <- cbind(continuous_gender_incidence_actual, continuous_gender_incidence_expected)
    
    continuous_age_incidence_actual <- matrix(tapply(continuous_algorithm_output$infected_by_trialend, continuous_algorithm_output$Age_Category, sum)/table(continuous_algorithm_output$Age_Category), nrow = 1)
    colnames(continuous_age_incidence_actual) <- c("less20_actual", "20to29_actual", "30to39_actual", "40to49_actual", "greater49_actual")
    continuous_age_incidence_expected <- matrix(tapply(continuous_algorithm_output$infected_probability, continuous_algorithm_output$Age_Category, mean), nrow = 1)
    colnames(continuous_age_incidence_expected) <- c("less20_expected", "20to29_expected", "30to39_expected", "40to49_expected", "greater49_expected")
    continuous_age_incidence <- cbind(continuous_age_incidence_actual, continuous_age_incidence_expected)
    
    continuous_demographic_incidence_row <- cbind(continuous_race_incidence, continuous_gender_incidence, continuous_age_incidence)
    continuous_demographic_incidence <- rbind(continuous_demographic_incidence, continuous_demographic_incidence_row)
  }
  
  print(paste("Completed runs:", nrow(batch_expectation_vs_reality)))
}

# Write output data to csv for future analysis
print("=== BATCH ALGORITHM RESULTS ===")
batch_fulldata <- cbind(batch_expectation_vs_reality, batch_race_analysis, batch_gender_analysis, batch_age_analysis, batch_demographic_incidence)
batch_fulldata_means <- colMeans(batch_fulldata)
print("Batch Algorithm - Means of all runs:")
print(batch_fulldata_means)

print("=== CONTINUOUS ALGORITHM RESULTS ===")
if(nrow(continuous_expectation_vs_reality) > 0) {
  continuous_fulldata <- cbind(continuous_expectation_vs_reality, continuous_race_analysis, continuous_gender_analysis, continuous_age_analysis, continuous_demographic_incidence)
  continuous_fulldata_means <- colMeans(continuous_fulldata)
  print("Continuous Algorithm - Means of all runs:")
  print(continuous_fulldata_means)
  
  print("=== ALGORITHM COMPARISON ===")
  print("Batch vs Continuous - Key Metrics:")
  print(paste("Batch actual incidence:", round(batch_fulldata_means["actual"], 4)))
  print(paste("Continuous actual incidence:", round(continuous_fulldata_means["actual"], 4)))
  print(paste("Batch agent count:", round(batch_fulldata_means["agent_count"], 0)))
  print(paste("Continuous agent count:", round(continuous_fulldata_means["agent_count"], 0)))
} else {
  print("No continuous algorithm results to analyze")
}

# Ensure results directory exists
if (!dir.exists("results")) {
  dir.create("results")
}

# save results to CSV files
write.csv(batch_fulldata, "results/batch_algorithm_results.csv", row.names = FALSE)
if(nrow(continuous_expectation_vs_reality) > 0) {
  write.csv(continuous_fulldata, "results/continuous_algorithm_results.csv", row.names = FALSE)
}