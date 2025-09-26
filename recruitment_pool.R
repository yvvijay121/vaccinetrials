library(magrittr)

# Set the working directory to the directory with the CSV of simulation data. The simulation data is derived from a 10-year run of the HepCEP algorithm
setwd("~/Documents/Programming_Projects/vaccinetrials")

# Read the CSV file as a data frame
simulation_data <- read.csv("data/events_pool.csv")

# Label data with the time each agent has their first infection, if applicable.
agents_infected_total <- simulation_data[simulation_data$Event == "infected",]
agents_infected_total <- agents_infected_total[!duplicated(agents_infected_total$Agent), c("Time", "Agent")]
colnames(agents_infected_total) <- c("infected_time", "Agent")

# Create recruitment_pool data frame with time of first infection and time of chronicity
recruitment_pool <- merge(simulation_data, agents_infected_total, by = "Agent", all.x = TRUE)

# Set recruitment pool only to unique agents when they are first activated
recruitment_pool <- recruitment_pool[recruitment_pool$Event == "activated",]

# Label if the agent was infected within 1.5 years of activation
recruitment_pool$infected_time_after_start <- recruitment_pool$infected_time - recruitment_pool$Time
recruitment_pool[is.na(recruitment_pool$infected_time_after_start),]$infected_time_after_start <- 999
recruitment_pool$infected_by_trialend <- 0
recruitment_pool[recruitment_pool$infected_time_after_start <= 1.5, ]$infected_by_trialend <- 1

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

# Label if each agent is susceptible or not for use in the algorithm. In real life use, PREDICTEE will not know if they are susceptible or not until they undergo a blood test, after which they will be removed and replaced with another agent who is not susceptible if they are determined to not be susceptible.
recruitment_pool$susceptible <- 0
recruitment_pool$susceptible[recruitment_pool$HCV == "susceptible"] <- 1