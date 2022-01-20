# BEFORE running algorithm, run the code in cox_model.R to load the Cox Model or algorithm_ML_model.R to load the RSF Model, and recruitment_pool.R to load the recruitment dataset and target population used in the algorithm.

# Load additional packages used in the algorithm
library(pwr)
library(tidyverse)

###########################################################################
## Find the demographic composition of the target population
###########################################################################
# Target pop set to susceptible PWID
target_pop <- cnep_susceptible
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


###########################################################################
##################### FINAL VERSION OF BATCH ALGORITHM ####################
###########################################################################

recruitment_per_batch <- 50
recruited_per_batch <- 5
trial_followup_years <- 1.5
attrition_prob <- 0.2
recruitment_dataset <- recruitment_pool
#model_used <- cox_model
model_used <- rsf_model

# Sample size constraints:
req_sample_size <- 800
work_constraint <- 8000
ssmethod <- "cohen"

# Setting recruitment parameters:
incidence_weight_min <- 25
high_demo_error_adjustment <- FALSE

# Create all the functions corresponding to each step in the algorithm:
apply_cox <- function(data, cox, follow_up_time) {
  data$status <- NA
  data$survival_time <- follow_up_time
  probabilities <- as.data.frame(predict(cox, data, type = "survival"))
  colnames(probabilities) <- "survival_probability"
  output <- cbind(data, probabilities)
  output$infected_probability <- 1-output$survival_probability
  
  return(output)
}

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

reestimate_n <- function(currently_in_trial, sig.level = 0.05, power = .80, method = "cohen") {
  if(method == "cohen"){
    p1 <- mean(currently_in_trial$infected_probability)
    p2 <- 0.4 * p1
    h <- 2*asin(sqrt(p1))-2*asin(sqrt(p2))
    new_n <- 2*ceiling(pwr.2p.test(h = h, sig.level = sig.level, power = power, alternative="greater")$n)
  }
  
  if(method == "kelsey"){
    z_a <- qnorm(sig.level/2, lower.tail = F)
    z_b <- qnorm(1-power, lower.tail = F)
    p1 <- mean(currently_in_trial$infected_probability)
    p2 <- 0.4 * p1
    p <- (p1+p2)/2
    new_n <- 2*ceiling((((z_a+z_b)^2)*(p)*(1-p)*2)/((p1-p2)^2))
  }
  
  if(method == "schoenfeld"){
    z_a <- qnorm(sig.level/2, lower.tail = F)
    z_b <- qnorm(1-power, lower.tail = F)
    theta <- 0.4
    d <- 4*(((z_a+z_b)^2)/(log(theta)^2))
    p1 <- mean(currently_in_trial$infected_probability)
    p2 <- 0.4 * p1
    p <- (p1+p2)/2
    new_n <- ceiling(d/p)
  }
  
  return(new_n)
}

print_warning <- function(work_constraint, agent_count, recruitment_per_batch, recruited_per_batch, algorithm_output){
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

# Reestimation point will be set by default at 1/2 work constraint for Cox, 1/3 for ML, but can be adjusted:
if(class(model_used) == "coxph"){
  reestimation_point <- work_constraint/2
}

if(class(model_used)[1] == "rfsrc"){
  reestimation_point <- work_constraint/3
}
reestimate_completion <- 0

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
  
  # Apply the Cox model to the agents recruited for the day
  if(class(model_used) == "coxph"){
    eligible_postmodel <- apply_cox(eligible, model_used, trial_followup_years)
  }
  
  if(class(model_used)[1] == "rfsrc"){
    eligible_postmodel <- apply_rsf(eligible, model_used, trial_followup_years)
  }
  
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
  
  # Assign a random number 0-100 to each agent to determine attrition
  backlog$attrition <- runif(nrow(backlog))
  
  # Use attrition probability to determine which agents get removed based on randomly generated number
  backlog <- backlog[backlog$attrition <= attrition_prob,]
  
  # Delete gender_diff, race_diff, demographic_score, score, and attrition from the backlog for proper looping
  backlog <- backlog[, -which(names(backlog) %in% c("gender_diff", "race_diff", "demographic_score", "score", "attrition"))]
  
  # Adjust weights after each batch, up to a certain limit
  if(incidence_weight > incidence_weight_min) {
    incidence_weight <- incidence_weight - weight_change_per_batch
    demographic_weight <- demographic_weight + weight_change_per_batch
    blinded_n <- reestimate_n(algorithm_output, method = ssmethod)
    weight_change_per_batch <- (incidence_weight-incidence_weight_min)/(blinded_n/recruited_per_batch)
    
    if(high_demo_error_adjustment == TRUE){
      error_adjustment <- 1-min(c(1-abs((table(algorithm_output$Gender)/nrow(algorithm_output)-target_gender_comp)/target_gender_comp),1-abs((table(algorithm_output$Race)/nrow(algorithm_output)-target_race_comp)/target_race_comp)))
      weight_change_per_batch <- weight_change_per_batch + error_adjustment
    }
    
    print(paste("Current Weights:", incidence_weight, demographic_weight, "Weight Change Next Batch:", weight_change_per_batch))
  }
  
  # Calculate new estimated sample size requirement at the prespecified reestimation point
  if(agent_count >= reestimation_point && reestimate_completion == 0) {
    req_sample_size <- reestimate_n(algorithm_output, method = ssmethod)
    print(paste("Sample size reestimated to be:", req_sample_size))
    reestimate_completion <- 1
  }
  
  # Print current agent count
  print(paste("Current number of agents screened:", agent_count))
  
  # Print warning statement in console if work constraint is destined to be violated
  print_warning(work_constraint, agent_count, recruitment_per_batch, recruited_par_batch, algorithm_output)
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

# Representativeness Max Error Calculation:
1-abs((table(algorithm_output$Gender)/nrow(algorithm_output)-target_gender_comp)/target_gender_comp)
1-abs((table(algorithm_output$Race)/nrow(algorithm_output)-target_race_comp)/target_race_comp)











#########################################################################
################################ FUNCTION ###############################
#########################################################################
algorithm <- function(recruitment_dataset, model_used, recruitment_per_batch, recruited_per_batch, trial_followup_years, req_sample_size, work_constraint, incidence_weight_min = 25, attrition_prob = 0.2, high_demo_error_adjustment = FALSE, ssmethod = "cohen") {
  start_time <- proc.time()
  
  # Reestimation point will be set by default at 1/2 work constraint for Cox, 1/3 for ML, but can be adjusted:
  if(class(model_used) == "coxph"){
    reestimation_point <- work_constraint/2
  }
  
  if(class(model_used)[1] == "rfsrc"){
    reestimation_point <- work_constraint/3
  }
  reestimate_completion <- 0
  
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
    
    # Apply the Cox model to the agents recruited for the day
    if(class(model_used)[1] == "coxph"){
      eligible_postmodel <- apply_cox(eligible, model_used, trial_followup_years)
    }
    
    if(class(model_used)[1] == "rfsrc"){
      eligible_postmodel <- apply_rsf(eligible, model_used, trial_followup_years)
    }
    
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
    
    # Assign a random number 0-100 to each agent to determine attrition
    backlog$attrition <- runif(nrow(backlog))
    
    # Use attrition probability to determine which agents get removed based on randomly generated number
    backlog <- backlog[backlog$attrition <= attrition_prob,]
    
    # Delete gender_diff, race_diff, demographic_score, score, and attrition from the backlog for proper looping
    backlog <- backlog[, -which(names(backlog) %in% c("gender_diff", "race_diff", "demographic_score", "score", "attrition"))]
    
    # Adjust weights after each batch, up to a certain limit
    if(incidence_weight > incidence_weight_min) {
      incidence_weight <- incidence_weight - weight_change_per_batch
      demographic_weight <- demographic_weight + weight_change_per_batch
      blinded_n <- reestimate_n(algorithm_output, method = ssmethod)
      weight_change_per_batch <- (incidence_weight-incidence_weight_min)/(blinded_n/5)
      
      if(high_demo_error_adjustment == TRUE){
        error_adjustment <- 1-min(c(1-abs((table(algorithm_output$Gender)/nrow(algorithm_output)-target_gender_comp)/target_gender_comp),1-abs((table(algorithm_output$Race)/nrow(algorithm_output)-target_race_comp)/target_race_comp)))
        weight_change_per_batch <- weight_change_per_batch + error_adjustment
      }
      
      print(paste("Current Weights:", incidence_weight, demographic_weight, "Weight Change Next Batch:", weight_change_per_batch))
    }
    
    # Calculate new estimated sample size requirement at the prespecified reestimation point
    if(agent_count >= reestimation_point && reestimate_completion == 0) {
      req_sample_size <- reestimate_n(algorithm_output, method = ssmethod)
      print(paste("Sample size reestimated to be:", req_sample_size))
      reestimate_completion <- 1
    }
    
    # Print current agent count
    print(paste("Current number of agents screened:", agent_count))
    
    # Print warning statement in console if work constraint is destined to be violated
    print_warning(work_constraint, agent_count, recruitment_per_batch, recruited_per_batch, algorithm_output)
  }
  
  # Print elapsed time for algorithm completion
  print(proc.time()-start_time)
  
  # Construct list to be returned
  output_list <- list("algorithm_output" = algorithm_output, "agent_count" = agent_count)
  return(output_list)
}


# Testing functionality of algorithm function
algorithm_list <- algorithm(recruitment_dataset = recruitment_pool, model_used = rsf_model, recruitment_per_batch = 50, recruited_per_batch = 5, trial_followup_years = 1.5, req_sample_size = 800, work_constraint = 8000)

algorithm_output <- algorithm_list$algorithm_output

# After running algorithm, calculate actual incidence based on simulation data and predicted incidence
table(algorithm_output$infected_by_trialend)[2]/nrow(algorithm_output)
mean(algorithm_output$infected_probability)

#########################################################################
#########################################################################
#########################################################################















#######################################################################
############################### ANALYSIS ##############################
#######################################################################
############ USED TO RUN MULTIPLE ALGORITHMS FOR ANALYSIS #############

# Initiation of output data frames:
expectation_vs_reality <- data.frame()
race_analysis <- data.frame()
gender_analysis <- data.frame()

number_of_runs <- 100

while(nrow(expectation_vs_reality) < number_of_runs) {
  
  algorithm_list <- algorithm(recruitment_dataset = recruitment_pool, model_used = cox_model, recruitment_per_batch = 50, recruited_per_batch = 5, trial_followup_years = 1.5, req_sample_size = 800, work_constraint = 8000)
  
  algorithm_output <- algorithm_list$algorithm_output
  
  newrow <- matrix(c(table(algorithm_output$infected_by_trialend)[2]/sum(table(algorithm_output$infected_by_trialend)), mean(algorithm_output$infected_probability), table(algorithm_output$chronic_by_trialend)[2]/sum(table(algorithm_output$chronic_by_trialend)), algorithm_list$agent_count), nrow = 1)
  colnames(newrow) <- c("actual", "expected", "chronic", "agent_count")
  expectation_vs_reality <- rbind(expectation_vs_reality, newrow)
  
  racerow <- matrix(table(algorithm_output$Race)/nrow(algorithm_output), nrow = 1)
  colnames(racerow) <- c("Hispanic", "NHBlack", "NHWhite", "Other")
  genderrow <- matrix(table(algorithm_output$Gender)/nrow(algorithm_output), nrow = 1)
  colnames(genderrow) <- c("Female", "Male")
  
  race_analysis <- rbind(race_analysis, racerow)
  gender_analysis <- rbind(gender_analysis, genderrow)
  
  print(paste("Completed runs:", nrow(expectation_vs_reality)))
}

# Incidence Calculations:
mean(expectation_vs_reality$actual)
mean(expectation_vs_reality$expected)

# Sample Size Calculations:
p1 <- mean(expectation_vs_reality$actual)
p2 <- 0.4 * p1
h <- 2*asin(sqrt(p1))-2*asin(sqrt(p2))
2*ceiling(pwr.2p.test(h = h, sig.level = 0.05, power = .80, alternative="greater")$n)

# Representativeness Max Error Calculation:
1-abs((colMeans(gender_analysis)-target_gender_comp)/target_gender_comp)
1-abs((colMeans(race_analysis)-target_race_comp)/target_race_comp)

# Mean Number of Agents Screened Before Reaching Req. SS
mean(expectation_vs_reality$agent_count)

## GRAPHS - Gender, Race
# Gender
gender_input <- colMeans(gender_analysis)
gender <- rbind(gender_input, target_gender_comp)
rownames(gender) = c("Algorithm", "Susceptible")
barplot(gender, xlab = "Gender", ylab = "Percent", beside=TRUE)
legend(x = "topleft", legend = c("Algorithm", "Susceptible"), fill = c("#4D4D4D", "#E6E6E6"), cex = 0.60)

# Race
race_input <- colMeans(race_analysis)
race <- rbind(race_input, target_race_comp)
rownames(race) = c("Algorithm", "Susceptible")
barplot(race,xlab = "Race", ylab = "Percent", beside=TRUE)
legend(x = "topleft", legend = c("Algorithm", "Susceptible"), fill = c("#4D4D4D", "#E6E6E6"), cex = 0.60)

