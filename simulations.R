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
setwd("C:/Users/richa/OneDrive - University of Illinois at Chicago/Stats/HCV")

# Read the CSV file as a data frame and create the training set using 20% of the data.
simulation_data <- read.csv("events_pool.csv")

# # Create results data frame - RUN ONCE
# fulldata_dblabel100 <- data.frame()

while(nrow(fulldata_dblabel100) < 100){
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
  expectation_vs_reality <- data.frame()
  race_analysis <- data.frame()
  gender_analysis <- data.frame()
  age_analysis <- data.frame()
  demographic_incidence <- data.frame()
  
  number_of_runs <- 100
  
  while(nrow(expectation_vs_reality) < number_of_runs) {
    
    # Run the actual algorithm and store raw results in algorithm_output
    algorithm_list <- algorithm(recruitment_dataset = recruitment_pool, model_used = cox_model, target_demographics = c("Gender", "Race", "Age_Category"), target_props = target_demo, recruitment_per_batch = 50, recruited_per_batch = 5, trial_followup_years = 1.5, req_sample_size = 800, work_constraint = 8000, attrition_prob = 0.2)
    # algorithm_list <- algorithm(recruitment_dataset = recruitment_pool, model_used = cox_model, target_demographics = NULL, target_props = target_demo, recruitment_per_batch = 50, recruited_per_batch = 5, trial_followup_years = 1.5, req_sample_size = 800, work_constraint = 8000)
    
    algorithm_output <- algorithm_list$algorithm_output
    
    # Calculate incidence results and agent count, and store them in a matrix
    newrow <- matrix(c(table(algorithm_output$infected_by_trialend)[2]/sum(table(algorithm_output$infected_by_trialend)), mean(algorithm_output$infected_probability), table(algorithm_output$chronic_by_trialend)[2]/sum(table(algorithm_output$chronic_by_trialend)), algorithm_list$agent_count), nrow = 1)
    colnames(newrow) <- c("actual", "expected", "chronic", "agent_count")
    expectation_vs_reality <- rbind(expectation_vs_reality, newrow)
    
    # Calculate demographic results and store them in a matrix
    racerow <- matrix(table(algorithm_output$Race)/nrow(algorithm_output), nrow = 1)
    colnames(racerow) <- c("Hispanic", "NHBlack", "NHWhite", "Other")
    genderrow <- matrix(table(algorithm_output$Gender)/nrow(algorithm_output), nrow = 1)
    colnames(genderrow) <- c("Female", "Male")
    agerow <- matrix(table(algorithm_output$Age_Category)/nrow(algorithm_output), nrow = 1)
    colnames(agerow) <- c("<20", "20-29", "30-39", "40-49", "49+")
    
    race_analysis <- rbind(race_analysis, racerow)
    gender_analysis <- rbind(gender_analysis, genderrow)
    age_analysis <- rbind(age_analysis, agerow)
    
    # Calculate incidence results within each demographic group and store them in a matrix
    race_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Race, sum)/table(algorithm_output$Race), nrow = 1)
    colnames(race_incidence_actual) <- c("hispanic_actual", "nhblack_actual", "nhwhite_actual", "other_actual")
    race_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Race, mean), nrow = 1)
    colnames(race_incidence_expected) <- c("hispanic_expected", "nhblack_expected", "nhwhite_expected", "other_expected")
    race_incidence <- cbind(race_incidence_actual, race_incidence_expected)
    
    gender_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Gender, sum)/table(algorithm_output$Gender), nrow = 1)
    colnames(gender_incidence_actual) <- c("female_actual", "male_actual")
    gender_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Gender, mean), nrow = 1)
    colnames(gender_incidence_expected) <- c("female_expected", "male_expected")
    gender_incidence <- cbind(gender_incidence_actual, gender_incidence_expected)
    
    age_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Age_Category, sum)/table(algorithm_output$Age_Category), nrow = 1)
    colnames(age_incidence_actual) <- c("less20_actual", "20to29_actual", "30to39_actual", "40to49_actual", "greater49_actual")
    age_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Age_Category, mean), nrow = 1)
    colnames(age_incidence_expected) <- c("less20_expected", "20to29_expected", "30to39_expected", "40to49_expected", "greater49_expected")
    age_incidence <- cbind(age_incidence_actual, age_incidence_expected)
    
    demographic_incidence_row <- cbind(race_incidence, gender_incidence, age_incidence)
    demographic_incidence <- rbind(demographic_incidence, demographic_incidence_row)
    
    print(paste("Completed runs:", nrow(expectation_vs_reality)))
  }
  
  # Write output data to csv for future analysis
  fulldata <- cbind(expectation_vs_reality, race_analysis, gender_analysis, age_analysis, demographic_incidence)
  fulldata_means <- colMeans(fulldata)
  fulldata_dblabel100 <- rbind(fulldata_dblabel100, fulldata_means)
}

colnames(fulldata_dblabel100) <- colnames(fulldata)
write.csv(fulldata_dblabel100, "cox_results_dblabel100-100_final.csv") # input file name












#######################################################################
############################### ANALYSIS ##############################
############################### RSF MODEL #############################
############ 100 RANDOM DBLABEL ASSIGNMENTS, 100 RUNS EACH ############
######################### WITH AGE DEMOGRAPHIC ########################
library(magrittr)
library(survival)
library(survminer)
library(dplyr)
library(pwr)
library(tidyverse)
library(randomForestSRC)

# Set the working directory to the directory with the CSV of data
setwd("C:/Users/richa/OneDrive - University of Illinois at Chicago/Stats/HCV")

# Read the CSV file as a data frame and create the training set using 20% of the data.
simulation_data <- read.csv("events_pool.csv")

# # Create results data frame
# fulldata_dblabel100 <- data.frame()
# fulldata_dblabel100 <- read.csv("rsf_results_dblabel100-100.csv")

while(nrow(fulldata_dblabel100) < 100){
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
  
  #### Create Random Survival Forest Model
  
  # Best model supported by Qiu et al. (A Comparison Study of Machine Learning (Random Survival Forest) and Classic Statistic (Cox Proportional Hazards) for Predicting Progression in High-Grade Glioma after Proton and Carbon Ion Radiotherapy)
  rsf_model <- rfsrc(Surv(time = survival_time, event = status) ~ Age + Gender + Syringe_source + Drug_in_degree + Drug_out_degree + current_total_network_size + Daily_injection_intensity + Fraction_recept_sharing + chicago_community_name, data = unique_agents, ntree = 100, mtry = 3)
  
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
  expectation_vs_reality <- data.frame()
  race_analysis <- data.frame()
  gender_analysis <- data.frame()
  age_analysis <- data.frame()
  demographic_incidence <- data.frame()
  
  number_of_runs <- 100
  
  while(nrow(expectation_vs_reality) < number_of_runs) {
    
    # Run the actual algorithm and store raw results in algorithm_output
    algorithm_list <- algorithm(recruitment_dataset = recruitment_pool, model_used = rsf_model, target_demographics = c("Gender", "Race", "Age_Category"), target_props = target_demo, recruitment_per_batch = 50, recruited_per_batch = 5, trial_followup_years = 1.5, req_sample_size = 800, work_constraint = 8000, attrition_prob = 0.2)
    # algorithm_list <- algorithm(recruitment_dataset = recruitment_pool, model_used = cox_model, target_demographics = NULL, target_props = target_demo, recruitment_per_batch = 50, recruited_per_batch = 5, trial_followup_years = 1.5, req_sample_size = 800, work_constraint = 8000)
    
    algorithm_output <- algorithm_list$algorithm_output
    
    # Calculate incidence results and agent count, and store them in a matrix
    newrow <- matrix(c(table(algorithm_output$infected_by_trialend)[2]/sum(table(algorithm_output$infected_by_trialend)), mean(algorithm_output$infected_probability), table(algorithm_output$chronic_by_trialend)[2]/sum(table(algorithm_output$chronic_by_trialend)), algorithm_list$agent_count), nrow = 1)
    colnames(newrow) <- c("actual", "expected", "chronic", "agent_count")
    expectation_vs_reality <- rbind(expectation_vs_reality, newrow)
    
    # Calculate demographic results and store them in a matrix
    racerow <- matrix(table(algorithm_output$Race)/nrow(algorithm_output), nrow = 1)
    colnames(racerow) <- c("Hispanic", "NHBlack", "NHWhite", "Other")
    genderrow <- matrix(table(algorithm_output$Gender)/nrow(algorithm_output), nrow = 1)
    colnames(genderrow) <- c("Female", "Male")
    agerow <- matrix(table(algorithm_output$Age_Category)/nrow(algorithm_output), nrow = 1)
    colnames(agerow) <- c("<20", "20-29", "30-39", "40-49", "49+")
    
    race_analysis <- rbind(race_analysis, racerow)
    gender_analysis <- rbind(gender_analysis, genderrow)
    age_analysis <- rbind(age_analysis, agerow)
    
    # Calculate incidence results within each demographic group and store them in a matrix
    race_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Race, sum)/table(algorithm_output$Race), nrow = 1)
    colnames(race_incidence_actual) <- c("hispanic_actual", "nhblack_actual", "nhwhite_actual", "other_actual")
    race_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Race, mean), nrow = 1)
    colnames(race_incidence_expected) <- c("hispanic_expected", "nhblack_expected", "nhwhite_expected", "other_expected")
    race_incidence <- cbind(race_incidence_actual, race_incidence_expected)
    
    gender_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Gender, sum)/table(algorithm_output$Gender), nrow = 1)
    colnames(gender_incidence_actual) <- c("female_actual", "male_actual")
    gender_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Gender, mean), nrow = 1)
    colnames(gender_incidence_expected) <- c("female_expected", "male_expected")
    gender_incidence <- cbind(gender_incidence_actual, gender_incidence_expected)
    
    age_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Age_Category, sum)/table(algorithm_output$Age_Category), nrow = 1)
    colnames(age_incidence_actual) <- c("less20_actual", "20to29_actual", "30to39_actual", "40to49_actual", "greater49_actual")
    age_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Age_Category, mean), nrow = 1)
    colnames(age_incidence_expected) <- c("less20_expected", "20to29_expected", "30to39_expected", "40to49_expected", "greater49_expected")
    age_incidence <- cbind(age_incidence_actual, age_incidence_expected)
    
    demographic_incidence_row <- cbind(race_incidence, gender_incidence, age_incidence)
    demographic_incidence <- rbind(demographic_incidence, demographic_incidence_row)
    
    print(paste("Completed runs:", nrow(expectation_vs_reality)))
  }
  
  # Write output data to csv for future analysis
  fulldata <- cbind(expectation_vs_reality, race_analysis, gender_analysis, age_analysis, demographic_incidence)
  fulldata_means <- colMeans(fulldata)
  fulldata_dblabel100 <- rbind(fulldata_dblabel100, fulldata_means)
}

colnames(fulldata_dblabel100) <- colnames(fulldata)
write.csv(fulldata_dblabel100, "rsf_results_dblabel100-100_final.csv") # input file name












#######################################################################
############################### ANALYSIS ##############################
############################### COX MODEL #############################
############ 100 RANDOM DBLABEL ASSIGNMENTS, 100 RUNS EACH ############
########################## EXTREME POPULATION #########################
library(magrittr)
library(survival)
library(survminer)
library(dplyr)
library(pwr)
library(tidyverse)

# Set the working directory to the directory with the CSV of data
setwd("C:/Users/richa/OneDrive - University of Illinois at Chicago/Stats/HCV")

# Read the CSV file as a data frame and create the training set using 20% of the data.
simulation_data <- read.csv("events_pool.csv")

# # Create results data frame
# fulldata_dblabel100 <- data.frame()
# fulldata_dblabel100 <- read.csv("cox_results_dblabel100-100_extreme.csv")

# Create the target demographic
target_demo <- matrix(c(0.5, 0.5, 0.33, 0.33, 0.33, 0.01), nrow = 1)
colnames(target_demo) <- c("Female", "Male", "Hispanic", "NHBlack", "NHWhite", "Other")

while(nrow(fulldata_dblabel100) < 100){
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
  expectation_vs_reality <- data.frame()
  race_analysis <- data.frame()
  gender_analysis <- data.frame()
  demographic_incidence <- data.frame()
  
  number_of_runs <- 100
  
  while(nrow(expectation_vs_reality) < number_of_runs) {
    
    # Run the actual algorithm and store raw results in algorithm_output
    algorithm_list <- algorithm(recruitment_dataset = recruitment_pool, model_used = cox_model, target_demographics = c("Gender", "Race"), target_props = target_demo, recruitment_per_batch = 50, recruited_per_batch = 5, trial_followup_years = 1.5, req_sample_size = 800, work_constraint = 8000, attrition_prob = 0.2)
    # algorithm_list <- algorithm(recruitment_dataset = recruitment_pool, model_used = cox_model, target_demographics = NULL, target_props = target_demo, recruitment_per_batch = 50, recruited_per_batch = 5, trial_followup_years = 1.5, req_sample_size = 800, work_constraint = 8000)
    
    algorithm_output <- algorithm_list$algorithm_output
    
    # Calculate incidence results and agent count, and store them in a matrix
    newrow <- matrix(c(table(algorithm_output$infected_by_trialend)[2]/sum(table(algorithm_output$infected_by_trialend)), mean(algorithm_output$infected_probability), table(algorithm_output$chronic_by_trialend)[2]/sum(table(algorithm_output$chronic_by_trialend)), algorithm_list$agent_count), nrow = 1)
    colnames(newrow) <- c("actual", "expected", "chronic", "agent_count")
    expectation_vs_reality <- rbind(expectation_vs_reality, newrow)
    
    # Calculate demographic results and store them in a matrix
    racerow <- matrix(table(algorithm_output$Race)/nrow(algorithm_output), nrow = 1)
    colnames(racerow) <- c("Hispanic", "NHBlack", "NHWhite", "Other")
    genderrow <- matrix(table(algorithm_output$Gender)/nrow(algorithm_output), nrow = 1)
    colnames(genderrow) <- c("Female", "Male")
    
    race_analysis <- rbind(race_analysis, racerow)
    gender_analysis <- rbind(gender_analysis, genderrow)
    
    # Calculate incidence results within each demographic group and store them in a matrix
    race_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Race, sum)/table(algorithm_output$Race), nrow = 1)
    colnames(race_incidence_actual) <- c("hispanic_actual", "nhblack_actual", "nhwhite_actual", "other_actual")
    race_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Race, mean), nrow = 1)
    colnames(race_incidence_expected) <- c("hispanic_expected", "nhblack_expected", "nhwhite_expected", "other_expected")
    race_incidence <- cbind(race_incidence_actual, race_incidence_expected)
    
    gender_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Gender, sum)/table(algorithm_output$Gender), nrow = 1)
    colnames(gender_incidence_actual) <- c("female_actual", "male_actual")
    gender_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Gender, mean), nrow = 1)
    colnames(gender_incidence_expected) <- c("female_expected", "male_expected")
    gender_incidence <- cbind(gender_incidence_actual, gender_incidence_expected)
    
    demographic_incidence_row <- cbind(race_incidence, gender_incidence)
    demographic_incidence <- rbind(demographic_incidence, demographic_incidence_row)
    
    print(paste("Completed runs:", nrow(expectation_vs_reality)))
  }
  
  # Write output data to csv for future analysis
  fulldata <- cbind(expectation_vs_reality, race_analysis, gender_analysis, demographic_incidence)
  fulldata_means <- colMeans(fulldata)
  fulldata_dblabel100 <- rbind(fulldata_dblabel100, fulldata_means)
}

colnames(fulldata_dblabel100) <- colnames(fulldata)
write.csv(fulldata_dblabel100, "cox_results_dblabel100-100_extreme_final.csv") # input file name












#######################################################################
############################### ANALYSIS ##############################
############################### RSF MODEL #############################
############ 100 RANDOM DBLABEL ASSIGNMENTS, 100 RUNS EACH ############
########################## EXTREME POPULATION #########################
library(magrittr)
library(survival)
library(survminer)
library(dplyr)
library(pwr)
library(tidyverse)
library(randomForestSRC)

# Set the working directory to the directory with the CSV of data
setwd("C:/Users/richa/OneDrive - University of Illinois at Chicago/Stats/HCV")

# Read the CSV file as a data frame and create the training set using 20% of the data.
simulation_data <- read.csv("events_pool.csv")

# # Create results data frame
# fulldata_dblabel100 <- data.frame()

# Create the target demographic
target_demo <- matrix(c(0.5, 0.5, 0.33, 0.33, 0.33, 0.01), nrow = 1)
colnames(target_demo) <- c("Female", "Male", "Hispanic", "NHBlack", "NHWhite", "Other")

while(nrow(fulldata_dblabel100) < 100){
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
  
  #### Create Random Survival Forest Model
  
  # Best model supported by Qiu et al. (A Comparison Study of Machine Learning (Random Survival Forest) and Classic Statistic (Cox Proportional Hazards) for Predicting Progression in High-Grade Glioma after Proton and Carbon Ion Radiotherapy)
  rsf_model <- rfsrc(Surv(time = survival_time, event = status) ~ Age + Gender + Syringe_source + Drug_in_degree + Drug_out_degree + current_total_network_size + Daily_injection_intensity + Fraction_recept_sharing + chicago_community_name, data = unique_agents, ntree = 100, mtry = 3)
  
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
  expectation_vs_reality <- data.frame()
  race_analysis <- data.frame()
  gender_analysis <- data.frame()
  demographic_incidence <- data.frame()
  
  number_of_runs <- 100
  
  while(nrow(expectation_vs_reality) < number_of_runs) {
    
    # Run the actual algorithm and store raw results in algorithm_output
    algorithm_list <- algorithm(recruitment_dataset = recruitment_pool, model_used = rsf_model, target_demographics = c("Gender", "Race"), target_props = target_demo, recruitment_per_batch = 50, recruited_per_batch = 5, trial_followup_years = 1.5, req_sample_size = 800, work_constraint = 8000, attrition_prob = 0.2)
    
    algorithm_output <- algorithm_list$algorithm_output
    
    # Calculate incidence results and agent count, and store them in a matrix
    newrow <- matrix(c(table(algorithm_output$infected_by_trialend)[2]/sum(table(algorithm_output$infected_by_trialend)), mean(algorithm_output$infected_probability), table(algorithm_output$chronic_by_trialend)[2]/sum(table(algorithm_output$chronic_by_trialend)), algorithm_list$agent_count), nrow = 1)
    colnames(newrow) <- c("actual", "expected", "chronic", "agent_count")
    expectation_vs_reality <- rbind(expectation_vs_reality, newrow)
    
    # Calculate demographic results and store them in a matrix
    racerow <- matrix(table(algorithm_output$Race)/nrow(algorithm_output), nrow = 1)
    colnames(racerow) <- c("Hispanic", "NHBlack", "NHWhite", "Other")
    genderrow <- matrix(table(algorithm_output$Gender)/nrow(algorithm_output), nrow = 1)
    colnames(genderrow) <- c("Female", "Male")
    
    race_analysis <- rbind(race_analysis, racerow)
    gender_analysis <- rbind(gender_analysis, genderrow)
    
    # Calculate incidence results within each demographic group and store them in a matrix
    race_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Race, sum)/table(algorithm_output$Race), nrow = 1)
    colnames(race_incidence_actual) <- c("hispanic_actual", "nhblack_actual", "nhwhite_actual", "other_actual")
    race_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Race, mean), nrow = 1)
    colnames(race_incidence_expected) <- c("hispanic_expected", "nhblack_expected", "nhwhite_expected", "other_expected")
    race_incidence <- cbind(race_incidence_actual, race_incidence_expected)
    
    gender_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Gender, sum)/table(algorithm_output$Gender), nrow = 1)
    colnames(gender_incidence_actual) <- c("female_actual", "male_actual")
    gender_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Gender, mean), nrow = 1)
    colnames(gender_incidence_expected) <- c("female_expected", "male_expected")
    gender_incidence <- cbind(gender_incidence_actual, gender_incidence_expected)
    
    demographic_incidence_row <- cbind(race_incidence, gender_incidence)
    demographic_incidence <- rbind(demographic_incidence, demographic_incidence_row)
    
    print(paste("Completed runs:", nrow(expectation_vs_reality)))
  }
  
  # Write output data to csv for future analysis
  fulldata <- cbind(expectation_vs_reality, race_analysis, gender_analysis, demographic_incidence)
  fulldata_means <- colMeans(fulldata)
  fulldata_dblabel100 <- rbind(fulldata_dblabel100, fulldata_means)
}

colnames(fulldata_dblabel100) <- colnames(fulldata)
write.csv(fulldata_dblabel100, "rsf_results_dblabel100-100_extreme_final.csv") # input file name









################# Unweighted Simulations subgroup incidence DBLabel Methodology - Cox #################
library(magrittr)
library(survival)
library(survminer)
library(dplyr)
library(pwr)
library(tidyverse)

# Set the working directory to the directory with the CSV of data
setwd("C:/Users/richa/OneDrive - University of Illinois at Chicago/Stats/HCV")

# Read the CSV file as a data frame and create the training set using 20% of the data.
simulation_data <- read.csv("events_pool.csv")

# Create results data frame
# fulldata_dblabel100 <- data.frame()

while(nrow(fulldata_dblabel100) < 100){
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
  expectation_vs_reality <- data.frame()
  race_analysis <- data.frame()
  gender_analysis <- data.frame()
  age_analysis <- data.frame()
  demographic_incidence <- data.frame()
  
  number_of_runs <- 100
  
  while(nrow(expectation_vs_reality) < number_of_runs) {
    
    # Run the actual algorithm and store raw results in algorithm_output
    algorithm_list <- algorithm(recruitment_dataset = recruitment_pool, model_used = cox_model, target_demographics = NULL, target_props = target_demo, recruitment_per_batch = 50, recruited_per_batch = 5, trial_followup_years = 1.5, req_sample_size = 800, work_constraint = 8000, attrition_prob = 0.2)
    
    algorithm_output <- algorithm_list$algorithm_output
    
    # Calculate incidence results and agent count, and store them in a matrix
    newrow <- matrix(c(table(algorithm_output$infected_by_trialend)[2]/sum(table(algorithm_output$infected_by_trialend)), mean(algorithm_output$infected_probability), table(algorithm_output$chronic_by_trialend)[2]/sum(table(algorithm_output$chronic_by_trialend)), algorithm_list$agent_count), nrow = 1)
    colnames(newrow) <- c("actual", "expected", "chronic", "agent_count")
    expectation_vs_reality <- rbind(expectation_vs_reality, newrow)
    
    # Calculate demographic results and store them in a matrix
    racerow <- matrix(table(algorithm_output$Race)/nrow(algorithm_output), nrow = 1)
    colnames(racerow) <- c("Hispanic", "NHBlack", "NHWhite", "Other")
    genderrow <- matrix(table(algorithm_output$Gender)/nrow(algorithm_output), nrow = 1)
    colnames(genderrow) <- c("Female", "Male")
    agerow <- matrix(table(algorithm_output$Age_Category)/nrow(algorithm_output), nrow = 1)
    colnames(agerow) <- c("<20", "20-29", "30-39", "40-49", "49+")
    
    race_analysis <- rbind(race_analysis, racerow)
    gender_analysis <- rbind(gender_analysis, genderrow)
    age_analysis <- rbind(age_analysis, agerow)
    
    # Calculate incidence results within each demographic group and store them in a matrix
    race_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Race, sum)/table(algorithm_output$Race), nrow = 1)
    colnames(race_incidence_actual) <- c("hispanic_actual", "nhblack_actual", "nhwhite_actual", "other_actual")
    race_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Race, mean), nrow = 1)
    colnames(race_incidence_expected) <- c("hispanic_expected", "nhblack_expected", "nhwhite_expected", "other_expected")
    race_incidence <- cbind(race_incidence_actual, race_incidence_expected)
    
    gender_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Gender, sum)/table(algorithm_output$Gender), nrow = 1)
    colnames(gender_incidence_actual) <- c("female_actual", "male_actual")
    gender_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Gender, mean), nrow = 1)
    colnames(gender_incidence_expected) <- c("female_expected", "male_expected")
    gender_incidence <- cbind(gender_incidence_actual, gender_incidence_expected)
    
    age_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Age_Category, sum)/table(algorithm_output$Age_Category), nrow = 1)
    colnames(age_incidence_actual) <- c("less20_actual", "20to29_actual", "30to39_actual", "40to49_actual", "greater49_actual")
    age_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Age_Category, mean), nrow = 1)
    colnames(age_incidence_expected) <- c("less20_expected", "20to29_expected", "30to39_expected", "40to49_expected", "greater49_expected")
    age_incidence <- cbind(age_incidence_actual, age_incidence_expected)
    
    demographic_incidence_row <- cbind(race_incidence, gender_incidence, age_incidence)
    demographic_incidence <- rbind(demographic_incidence, demographic_incidence_row)
    
    print(paste("Completed runs:", nrow(expectation_vs_reality)))
  }
  
  # Write output data to csv for future analysis
  fulldata <- cbind(expectation_vs_reality, race_analysis, gender_analysis, age_analysis, demographic_incidence)
  fulldata_means <- colMeans(fulldata)
  fulldata_dblabel100 <- rbind(fulldata_dblabel100, fulldata_means)
}

colnames(fulldata_dblabel100) <- colnames(fulldata)
write.csv(fulldata_dblabel100, "cox_subgroup_unweighted_dblabel_final.csv") # input file name












################# Unweighted Simulations subgroup incidence DBLabel Methodology - RSF #################
library(magrittr)
library(survival)
library(survminer)
library(dplyr)
library(pwr)
library(tidyverse)
library(randomForestSRC)

# Set the working directory to the directory with the CSV of data
setwd("C:/Users/richa/OneDrive - University of Illinois at Chicago/Stats/HCV")

# Read the CSV file as a data frame and create the training set using 20% of the data.
simulation_data <- read.csv("events_pool.csv")

# # Create results data frame
# fulldata_dblabel100 <- data.frame()

while(nrow(fulldata_dblabel100) < 100){
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
  
  #### Create Random Survival Forest Model
  
  # Best model supported by Qiu et al. (A Comparison Study of Machine Learning (Random Survival Forest) and Classic Statistic (Cox Proportional Hazards) for Predicting Progression in High-Grade Glioma after Proton and Carbon Ion Radiotherapy)
  rsf_model <- rfsrc(Surv(time = survival_time, event = status) ~ Age + Gender + Syringe_source + Drug_in_degree + Drug_out_degree + current_total_network_size + Daily_injection_intensity + Fraction_recept_sharing + chicago_community_name, data = unique_agents, ntree = 100, mtry = 3)
  
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
  expectation_vs_reality <- data.frame()
  race_analysis <- data.frame()
  gender_analysis <- data.frame()
  age_analysis <- data.frame()
  demographic_incidence <- data.frame()
  
  number_of_runs <- 100
  
  while(nrow(expectation_vs_reality) < number_of_runs) {
    
    # Run the actual algorithm and store raw results in algorithm_output
    algorithm_list <- algorithm(recruitment_dataset = recruitment_pool, model_used = rsf_model, target_demographics = NULL, target_props = target_demo, recruitment_per_batch = 50, recruited_per_batch = 5, trial_followup_years = 1.5, req_sample_size = 800, work_constraint = 8000, attrition_prob = 0.2)
    
    algorithm_output <- algorithm_list$algorithm_output
    
    # Calculate incidence results and agent count, and store them in a matrix
    newrow <- matrix(c(table(algorithm_output$infected_by_trialend)[2]/sum(table(algorithm_output$infected_by_trialend)), mean(algorithm_output$infected_probability), table(algorithm_output$chronic_by_trialend)[2]/sum(table(algorithm_output$chronic_by_trialend)), algorithm_list$agent_count), nrow = 1)
    colnames(newrow) <- c("actual", "expected", "chronic", "agent_count")
    expectation_vs_reality <- rbind(expectation_vs_reality, newrow)
    
    # Calculate demographic results and store them in a matrix
    racerow <- matrix(table(algorithm_output$Race)/nrow(algorithm_output), nrow = 1)
    colnames(racerow) <- c("Hispanic", "NHBlack", "NHWhite", "Other")
    genderrow <- matrix(table(algorithm_output$Gender)/nrow(algorithm_output), nrow = 1)
    colnames(genderrow) <- c("Female", "Male")
    agerow <- matrix(table(algorithm_output$Age_Category)/nrow(algorithm_output), nrow = 1)
    colnames(agerow) <- c("<20", "20-29", "30-39", "40-49", "49+")
    
    race_analysis <- rbind(race_analysis, racerow)
    gender_analysis <- rbind(gender_analysis, genderrow)
    age_analysis <- rbind(age_analysis, agerow)
    
    # Calculate incidence results within each demographic group and store them in a matrix
    race_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Race, sum)/table(algorithm_output$Race), nrow = 1)
    colnames(race_incidence_actual) <- c("hispanic_actual", "nhblack_actual", "nhwhite_actual", "other_actual")
    race_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Race, mean), nrow = 1)
    colnames(race_incidence_expected) <- c("hispanic_expected", "nhblack_expected", "nhwhite_expected", "other_expected")
    race_incidence <- cbind(race_incidence_actual, race_incidence_expected)
    
    gender_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Gender, sum)/table(algorithm_output$Gender), nrow = 1)
    colnames(gender_incidence_actual) <- c("female_actual", "male_actual")
    gender_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Gender, mean), nrow = 1)
    colnames(gender_incidence_expected) <- c("female_expected", "male_expected")
    gender_incidence <- cbind(gender_incidence_actual, gender_incidence_expected)
    
    age_incidence_actual <- matrix(tapply(algorithm_output$infected_by_trialend, algorithm_output$Age_Category, sum)/table(algorithm_output$Age_Category), nrow = 1)
    colnames(age_incidence_actual) <- c("less20_actual", "20to29_actual", "30to39_actual", "40to49_actual", "greater49_actual")
    age_incidence_expected <- matrix(tapply(algorithm_output$infected_probability, algorithm_output$Age_Category, mean), nrow = 1)
    colnames(age_incidence_expected) <- c("less20_expected", "20to29_expected", "30to39_expected", "40to49_expected", "greater49_expected")
    age_incidence <- cbind(age_incidence_actual, age_incidence_expected)
    
    demographic_incidence_row <- cbind(race_incidence, gender_incidence, age_incidence)
    demographic_incidence <- rbind(demographic_incidence, demographic_incidence_row)
    
    print(paste("Completed runs:", nrow(expectation_vs_reality)))
  }
  
  # Write output data to csv for future analysis
  fulldata <- cbind(expectation_vs_reality, race_analysis, gender_analysis, age_analysis, demographic_incidence)
  fulldata_means <- colMeans(fulldata)
  fulldata_dblabel100 <- rbind(fulldata_dblabel100, fulldata_means)
}

colnames(fulldata_dblabel100) <- colnames(fulldata)
write.csv(fulldata_dblabel100, "rsf_subgroup_unweighted_dblabel_final.csv") # input file name









#######################################################################
############################### ANALYSIS ##############################
############################### COX MODEL #############################
############ 100 RANDOM DBLABEL ASSIGNMENTS, 100 RUNS EACH ############
######################### C-INDEX CALCULATIONS ########################
library(magrittr)
library(survival)
library(survminer)
library(dplyr)
library(pwr)
library(tidyverse)
library(pec)

# Set the working directory to the directory with the CSV of data
setwd("C:/Users/richa/OneDrive - University of Illinois at Chicago/Stats/HCV")

# Read the CSV file as a data frame and create the training set using 20% of the data.
simulation_data <- read.csv("events_pool.csv")

# Create results data frame
# cox_cindex_dblabel100 <- data.frame()

while(nrow(cox_cindex_dblabel100) < 10000){
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
  
  # Reset recruitment_pool from prior run
  recruitment_pool <- NA
  
  # Create recruitment_pool data frame with time of first infection and time of chronicity
  recruitment_pool <- merge(test, agents_infected_total, by = "Agent", all.x = TRUE)
  recruitment_pool <- merge(recruitment_pool, agents_chronic_total, by = "Agent", all.x = TRUE)
  
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
  
  # Set status to infection by trial end
  recruitment_pool$status <- recruitment_pool$infected_by_trialend
  recruitment_pool$survival_time <- recruitment_pool$infected_time_after_start
  
  # Find concordance of model using training data
  train_cindex <- cindex(cox_model)$AppCindex
  
  # Find concordance of model using test data
  test_cindex <- cindex(cox_model, data = recruitment_pool)$AppCindex
  
  # Write output data to csv for future analysis
  cox_cindex_dblabel100 <- rbind(cox_cindex_dblabel100, c(train_cindex, test_cindex))
  
  print(paste("Completed runs:", nrow(cox_cindex_dblabel100)))
  write.csv(cox_cindex_dblabel100, "cox_cindex_10000_final.csv") # input file name
}














#######################################################################
############################### ANALYSIS ##############################
############################### RSF MODEL #############################
############ 100 RANDOM DBLABEL ASSIGNMENTS, 100 RUNS EACH ############
######################### C-INDEX CALCULATIONS ########################
library(magrittr)
library(survival)
library(survminer)
library(dplyr)
library(pwr)
library(tidyverse)
library(randomForestSRC)
library(pec)

# Set the working directory to the directory with the CSV of data
setwd("C:/Users/richa/OneDrive - University of Illinois at Chicago/Stats/HCV")

# Read the CSV file as a data frame and create the training set using 20% of the data.
simulation_data <- read.csv("events_pool.csv")

# Create results data frame
# rsf_cindex_dblabel100 <- data.frame()

while(nrow(rsf_cindex_dblabel100) < 10000){
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
  
  #### Create Random Survival Forest Model
  
  # Best model supported by Qiu et al. (A Comparison Study of Machine Learning (Random Survival Forest) and Classic Statistic (Cox Proportional Hazards) for Predicting Progression in High-Grade Glioma after Proton and Carbon Ion Radiotherapy)
  rsf_model <- rfsrc(Surv(time = survival_time, event = status) ~ Age + Gender + Syringe_source + Drug_in_degree + Drug_out_degree + current_total_network_size + Daily_injection_intensity + Fraction_recept_sharing + chicago_community_name, data = unique_agents, ntree = 100, mtry = 3)
  
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
  
  # Reset recruitment_pool from prior run
  recruitment_pool <- NA
  
  # Create recruitment_pool data frame with time of first infection and time of chronicity
  recruitment_pool <- merge(test, agents_infected_total, by = "Agent", all.x = TRUE)
  recruitment_pool <- merge(recruitment_pool, agents_chronic_total, by = "Agent", all.x = TRUE)
  
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
  
  recruitment_pool$status <- recruitment_pool$infected_by_trialend
  recruitment_pool$survival_time <- recruitment_pool$infected_time_after_start
  
  # Find concordance of model using training data
  train_cindex <- cindex(rsf_model)$AppCindex
  
  # Find concordance of model using test data
  test_cindex <- cindex(rsf_model, data = recruitment_pool)$AppCindex
  
  # Write output data to csv for future analysis
  rsf_cindex_dblabel100 <- rbind(rsf_cindex_dblabel100, c(train_cindex, test_cindex))
  
  print(paste("Completed runs:", nrow(rsf_cindex_dblabel100)))
  
  if(nrow(rsf_cindex_dblabel100) %% 100 == 0){
    write.csv(rsf_cindex_dblabel100, "rsf_cindex_dblabel10000_final.csv") # input file name
  }
}














#######################################################################