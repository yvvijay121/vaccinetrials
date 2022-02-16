library(magrittr)

# NOTE: the recruitment pool is derived from the test data used in the validation of the Cox Model; see cox_model.R

# Set the working directory to the directory with the CSV of data
setwd("C:/Users/richa/OneDrive - University of Illinois at Chicago/Stats/HCV")

# Read the CSV file as a data frame and create the training set using 20% of the data.
simulation_data <- read.csv("events_pool.csv")

# Create the same test data as seen in the Cox Model validation
test <- simulation_data[mod(simulation_data$DBLabel, 5) != 0,]

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







# Read CNEP data for use in demographic comparisons and as target population
cnep <- read.csv("cnep_plus_all_2018.02.13.csv")
cnep_susceptible <- cnep[cnep$HCV == "susceptible",]

# Stratify age categories in CNEP data for demographic matching
cnep_susceptible$Age_Category <- NA
cnep_susceptible[cnep_susceptible$Age < 20,]$Age_Category <- "<20"
cnep_susceptible[cnep_susceptible$Age >= 20 & cnep_susceptible$Age < 30,]$Age_Category <- "20-29"
cnep_susceptible[cnep_susceptible$Age >= 30 & cnep_susceptible$Age < 40,]$Age_Category <- "30-39"
cnep_susceptible[cnep_susceptible$Age >= 40 & cnep_susceptible$Age < 50,]$Age_Category <- "40-49"
cnep_susceptible[cnep_susceptible$Age >= 50,]$Age_Category <- "49+"