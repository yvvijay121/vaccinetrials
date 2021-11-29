library(randomForestSRC)
library(magrittr)
library(survival)
library(pROC)
library(pwr)

# Set the working directory to the directory with the CSV of data
setwd("C:/Users/richa/OneDrive - University of Illinois at Chicago/Stats/HCV")

# Read the CSV file as a data frame and create the training set using 20% of the data.
simulation_data <- read.csv("events_pool.csv")
train <- simulation_data[mod(simulation_data$DBLabel, 5) == 0,]

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
# Standard: nodesize set to 10% of the training data sample size
rf <- rfsrc(Surv(time = survival_time, event = status) ~ Age + Gender + Race + Syringe_source + Drug_in_degree + Drug_out_degree + current_total_network_size + Daily_injection_intensity + Fraction_recept_sharing + chicago_community_name, data = unique_agents, nodesize = round(nrow(unique_agents)/10, digits = 0))

# Best model: Not supported by literature, but produces best results for algorithm
rf_best <- rfsrc(Surv(time = survival_time, event = status) ~ Age + Gender + Race + Syringe_source + Drug_in_degree + Drug_out_degree + current_total_network_size + Daily_injection_intensity + Fraction_recept_sharing + chicago_community_name, data = unique_agents, nodesize = 150)









########################################################################
########################## Test for RSF Model ##########################
########################################################################
# Create test data frame using the remaining 80% of data not used to train the Cox Model
test <- simulation_data[mod(simulation_data$DBLabel, 5) != 0,]

# Label test data with the time each agent has their first infection, if applicable.
agents_infected_total <- simulation_data[simulation_data$Event == "infected",]
agents_infected_total <- agents_infected_total[!duplicated(agents_infected_total$Agent), c("Time", "Agent")]
colnames(agents_infected_total) <- c("infected_time", "Agent")
test <- merge(test, agents_infected_total, by = "Agent", all.x = TRUE)

# Reduce test data to unique agents at the time they were activated and only those who are HCV susceptible
test <- test[test$Event == "activated",]
test <- test[test$HCV == "susceptible",]

# Calculate new column for the test dataset with how long after agent activation they become infected; if they never become infected, set as 999 as placeholder value since working with NA values is inconvenient
test$infected_time_after_start <- test$infected_time - test$Time
test[is.na(test$infected_time_after_start),]$infected_time_after_start <- 999

# Using the infected time, determine if the agent would be infected by the end of a trial with a 1.5 year follow-up.
test$infected_by_trialend <- 0
test[test$infected_time_after_start <= 1.5, ]$infected_by_trialend <- 1

# Reclassify test data variables as factors if applicable
test$Gender <- factor(test$Gender)
test$Race <- factor(test$Race)
test$Syringe_source <- factor(test$Syringe_source)
test$chicago_community_name <- factor(test$chicago_community_name)

# Apply the trained Cox Model to the test set using the predict function, producing a dataframe of survival probabilities (the probability they do NOT become infected)
prediction_data <- predict.rfsrc(rf_best, test, na.action = "na.impute")
time_interest <- which(abs(prediction_data$time.interest-1.5)==min(abs(prediction_data$time.interest-1.5)))
probabilities <- 1 - prediction_data$survival[,time_interest]

# Combine the data frame of infection probabilities with the test data frame
test <- cbind(test, probabilities)
colnames(test)[colnames(test) == "probabilities"] <- "infected_probability"

# RESULTS to show the accuracy of the Cox Model in predicting:
# Calculates the mean predicted incidence for the test data using the Cox Model as well as the actual incidence based on the simulation.
mean(test$infected_probability)
table(test$infected_by_trialend)[2] / (table(test$infected_by_trialend)[1] + table(test$infected_by_trialend)[2])

# Evaluation of the RSF model using AUC
roc(data = test, response = infected_by_trialend, predictor = infected_probability)












########################################################################
##################### Algorithm incorporating RSF ######################
########################################################################

###########################################################################
## Find the demographic composition of the target population
###########################################################################
# Target pop set to susceptible PWID
target_pop <- recruitment_pool[recruitment_pool$susceptible == 1,]
target_gender_comp <- table(target_pop$Gender)/nrow(target_pop)
target_race_comp <- table(target_pop$Race)/nrow(target_pop)

# Function to create custom demographic targets
create_demographic_target <- function(race_numbers, gender_numbers){
  hispanic_raw <- race_numbers[1]
  nhblack_raw <- race_numbers[2]
  nhwhite_raw <- race_numbers[3]
  other_raw <- race_numbers[4]
  
  female_raw <- gender_numbers[1]
  male_raw <- gender_numbers[2]
  
  race_df <- as.data.frame(matrix(data = c(rep.int("Hispanic", hispanic_raw), rep.int("NHBlack", nhblack_raw), rep.int("NHWhite", nhwhite_raw), rep.int("Other", other_raw)), ncol = 1))
  
  gender_df <- as.data.frame(matrix(data = c(rep.int("Female", female_raw), rep.int("Male", male_raw)), ncol = 1))
  
  output_race_table <- table(race_df)
  output_gender_table <- table(gender_df)
  
  return_list <- list("race" = output_race_table, "gender" = output_gender_table)
  
  return(return_list)
}

## Target pop set to CDC 2019 Surveillance Report data
# target_pop <- create_demographic_target(c(350, 267, 2683, 119), c(1653, 2471))
# target_gender_comp <- target_pop$gender/sum(target_pop$gender)
# target_race_comp <- target_pop$race/sum(target_pop$race)

## Target pop set to test case for extreme differences between target population and available population
# target_pop <- create_demographic_target(c(40, 30, 20, 10), c(50, 50))
# target_gender_comp <- target_pop$gender/sum(target_pop$gender)
# target_race_comp <- target_pop$race/sum(target_pop$race)

######################### Setting initial parameters #########################
recruitment_per_batch <- 50
recruited_per_batch <- 5
trial_followup_years <- 1.5
backlog_batch_limit <- 5
recruitment_dataset <- recruitment_pool
rsf_model_used <- rf_best
algorithm_output <- data.frame()
backlog <- data.frame()
agent_count <- 0
batch_count <- 0

# Setting the initial recruitment weights:
incidence_weight <- 100
demographic_weight <- 0
incidence_weight_min <- 25
weight_change_per_batch <- (incidence_weight-incidence_weight_min)/(req_sample_size/recruited_per_batch)
high_demo_error_adjustment <- FALSE

# Sample size constraints:
req_sample_size <- 800
work_constraint <- 8000
# Reestimation point will be set by default at half the work constraint, but will be able to be adjusted:
reestimation_point <- work_constraint/3
reestimate_completion <- 0

# Create all the functions corresponding to each step in the algorithm:
apply_rsf <- function(data, model, follow_up_time) {
  prediction_data <- predict.rfsrc(model, data, na.action = "na.impute")
  time_interest <- which(abs(prediction_data$time.interest-follow_up_time)==min(abs(prediction_data$time.interest-follow_up_time)))
  probabilities <- as.data.frame(prediction_data$survival[,time_interest])
  colnames(probabilities) <- "survival_probability"
  output <- cbind(data, probabilities)
  output$infected_probability <- 1-output$survival_probability
  
  return(output)
}

compute_demographic_score <- function(data, currently_in_trial, target_race, target_gender) {
  output_gender_comp <- table(currently_in_trial$Gender)/nrow(currently_in_trial)
  output_race_comp <- table(currently_in_trial$Race)/nrow(currently_in_trial)
  if(nrow(currently_in_trial) == 0){
    gender_diff <- target_gender
    race_diff <- target_race
  } else {
    gender_diff <- target_gender - output_gender_comp
    race_diff <- target_race - output_race_comp
  }
  gender_diff <- as.data.frame(gender_diff)
  colnames(gender_diff) <- c("Gender", "gender_diff")
  race_diff <- as.data.frame(race_diff)
  colnames(race_diff) <- c("Race", "race_diff")
  
  output <- merge(data, gender_diff, by = "Gender")
  output <- merge(output, race_diff, by = "Race")
  output$demographic_score <- output$gender_diff + output$race_diff
  
  return(output)
}

reestimate_n <- function(currently_in_trial) {
  p1 <- mean(currently_in_trial$infected_probability)
  p2 <- 0.4 * p1
  h <- 2*asin(sqrt(p1))-2*asin(sqrt(p2))
  new_n <- 2*ceiling(pwr.2p.test(h = h, sig.level = 0.05, power = .80, alternative="greater")$n)
  
  return(new_n)
}

print_warning <- function(work_constraint, agent_count, recruitment_per_batch, algorithm_output){
  work_remaining <- work_constraint - agent_count
  remaining_agents_needed <- reestimate_n(algorithm_output) - nrow(algorithm_output)
  remaining_batches_needed <- remaining_agents_needed/recruited_per_batch
  if((remaining_batches_needed * recruitment_per_batch) > work_remaining){
    print("WARNING: With current predicted incidence, work constraint is destined to be exceeded.")
  }
}

#########################################################################
############################### ALGORITHM ###############################
#########################################################################

start_time <- proc.time()

while(nrow(algorithm_output) < req_sample_size) {
  # Sample recruitment_per_day number of agents at a time
  eligible <- recruitment_dataset[sample(nrow(recruitment_dataset), recruitment_per_batch),]
  agent_count <- agent_count + nrow(eligible)
  batch_count <- batch_count + 1
  
  # Apply the Cox model to the agents recruited for the day
  eligible_postmodel <- apply_rsf(eligible, rsf_model_used, trial_followup_years)
  
  # Eliminate agents who are not susceptible
  eligible_postmodel <- eligible_postmodel[eligible_postmodel$susceptible == 1,]
  
  # Add batches_elapsed_backlog to each agent to denote how many batches have passed since they have been in backlog. Initial value = 0.
  eligible_postmodel$batches_elapsed_backlog <- 0
  
  # Combine this batch and backlog into total_considered
  total_considered <- rbind(backlog, eligible_postmodel)
  
  # Calculate demographic difference scores for each agent, both race and gender
  total_considered <- compute_demographic_score(total_considered, algorithm_output, target_race_comp, target_gender_comp)
  
  # Apply weights to generate a single score
  total_considered$score <- ((total_considered$infected_probability * incidence_weight) + (total_considered$demographic_score * demographic_weight))
  
  # Rank agents by overall score
  total_considered <- total_considered[order(total_considered$score, decreasing = TRUE),]
  
  # Recruit the top [recruited_per_batch] agents from total_considered, and add the rest back to the backlog
  recruited <- head(total_considered, recruited_per_batch)
  algorithm_output <- rbind(algorithm_output, recruited)
  
  backlog <- tail(total_considered, nrow(total_considered)-recruited_per_batch)
  
  # Delete gender_diff, race_diff, demographic_score, and score from the backlog for proper looping
  backlog <- backlog[, -which(names(backlog) %in% c("gender_diff", "race_diff", "demographic_score", "score"))]
  
  # Increase batches_elapsed_backlog for agents left in backlog
  backlog$batches_elapsed_backlog <- backlog$batches_elapsed_backlog + 1
  
  # Eliminate agents who have been in backlog for more than backlog_batch_limit batches
  backlog <- backlog[backlog$batches_elapsed_backlog <= backlog_batch_limit,]
  
  # Adjust weights after each batch, up to a certain limit
  if(incidence_weight > incidence_weight_min) {
    incidence_weight <- incidence_weight - weight_change_per_batch
    demographic_weight <- demographic_weight + weight_change_per_batch
    blinded_n <- reestimate_n(algorithm_output)
    weight_change_per_batch <- (incidence_weight-incidence_weight_min)/(blinded_n/5)
    
    if(high_demo_error_adjustment == TRUE){
      error_adjustment <- 1-min(c(1-abs((table(algorithm_output$Gender)/nrow(algorithm_output)-target_gender_comp)/target_gender_comp),1-abs((table(algorithm_output$Race)/nrow(algorithm_output)-target_race_comp)/target_race_comp)))
      weight_change_per_batch <- weight_change_per_batch + error_adjustment
    }
    
    print(paste("Current Weights:", incidence_weight, demographic_weight, "Weight Change Next Batch:", weight_change_per_batch))
  }
  
  # Calculate new estimated sample size requirement at the prespecified reestimation point
  if(agent_count >= reestimation_point && reestimate_completion == 0) {
    req_sample_size <- reestimate_n(algorithm_output)
    print(paste("Sample size reestimated to be:", req_sample_size))
    reestimate_completion <- 1
  }
  
  # Print current agent count
  print(paste("Current number of agents screened:", agent_count))
  
  # Print warning statement in console if work constraint is destined to be violated
  print_warning(work_constraint, agent_count, recruitment_per_batch, algorithm_output)
}

# Print elapsed time for algorithm completion
print(proc.time()-start_time)

#########################################################################
#########################################################################
#########################################################################




# After running algorithm, calculate actual incidence based on simulation data and predicted incidence
table(algorithm_output$infected_by_trialend)[2]/nrow(algorithm_output)
mean(algorithm_output$infected_probability)

## GRAPHS - Gender, Race
# Gender
gender <- rbind(table(algorithm_output$Gender)/nrow(algorithm_output), target_gender_comp)
rownames(gender) = c("Algorithm", "Susceptible")
barplot(gender, xlab = "Gender", ylab = "Percent", beside=TRUE)
legend(x = "topleft", legend = c("Algorithm", "Susceptible"), fill = c("#4D4D4D", "#E6E6E6"), cex = 0.60)

# Race
race <- rbind(table(algorithm_output$Race)/nrow(algorithm_output), target_race_comp)
rownames(race) = c("Algorithm", "Susceptible")
barplot(race,xlab = "Race", ylab = "Percent", beside=TRUE)
legend(x = "topleft", legend = c("Algorithm", "Susceptible"), fill = c("#4D4D4D", "#E6E6E6"), cex = 0.60)