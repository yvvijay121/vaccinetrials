library(BART)
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

# Create BART Model
bart_model <- surv.bart()

