# Debug script to understand the Cox model issue
library(survival)

# Load the Cox model
cox_model <- readRDS("universal_cox.rds")

# Check the model structure
print("=== COX MODEL DIAGNOSTICS ===")
print("Class of cox_model:")
print(class(cox_model))

print("\nFormula:")
print(formula(cox_model))

print("\nVariable names in coefficients:")
print(names(cox_model$coefficients))

print("\nModel terms:")
print(terms(cox_model))

print("\nModel call:")
print(cox_model$call)

# Check if there are any issues with the model object
print("\nModel environment:")
print(environment(cox_model$terms))

# Try to see what's in the model's data
if("x" %in% names(cox_model)) {
  print("\nModel has 'x' component (design matrix)")
  print("Colnames of x:")
  print(colnames(cox_model$x))
} else {
  print("\nModel does not have 'x' component")
}

# Load CNEP data to create a small test dataset
cnep <- read.csv("data/cnep_plus_all_2018.02.13.csv")
print(paste("\nCNEP data loaded with", nrow(cnep), "rows"))
print("CNEP columns:")
print(colnames(cnep))

# Create a small test dataset with the required variables
test_data <- cnep[1:5, ]
test_data$Agent <- 1:5

# Add missing variables that the Cox model needs
required_vars <- c("Age", "Gender", "Syringe_source", "Drug_in_degree", "Drug_out_degree", 
                   "current_total_network_size", "Daily_injection_intensity", 
                   "Fraction_recept_sharing", "chicago_community_name")

print("\nChecking for required variables in CNEP data:")
for(var in required_vars) {
  if(var %in% colnames(test_data)) {
    print(paste(var, ": PRESENT"))
  } else {
    print(paste(var, ": MISSING - will add dummy data"))
  }
}

# Add missing variables
if(!"Drug_in_degree" %in% colnames(test_data)) {
  test_data$Drug_in_degree <- c(1,2,0,3,1)
}
if(!"Drug_out_degree" %in% colnames(test_data)) {
  test_data$Drug_out_degree <- c(2,1,1,2,0) 
}
if(!"current_total_network_size" %in% colnames(test_data)) {
  test_data$current_total_network_size <- c(5,8,3,12,4)
}
if(!"Daily_injection_intensity" %in% colnames(test_data)) {
  test_data$Daily_injection_intensity <- c(1.2, 2.1, 0.8, 3.4, 1.7)
}
if(!"Fraction_recept_sharing" %in% colnames(test_data)) {
  test_data$Fraction_recept_sharing <- c(0.1, 0.3, 0.0, 0.5, 0.2)
}

# Make sure factors are properly set
test_data$Gender <- factor(test_data$Gender)
test_data$Race <- factor(test_data$Race)

if(!"Syringe_source" %in% colnames(test_data)) {
  test_data$Syringe_source <- factor(c("pharmacy", "exchange", "other", "pharmacy", "exchange"))
} else {
  test_data$Syringe_source <- factor(test_data$Syringe_source)
}

if(!"chicago_community_name" %in% colnames(test_data)) {
  test_data$chicago_community_name <- factor(c("North", "South", "West", "North", "South"))
} else {
  test_data$chicago_community_name <- factor(test_data$chicago_community_name)
}

print("\nTest data created:")
print("Dimensions:")
print(dim(test_data))
print("Column names:")
print(colnames(test_data))

# Now try to predict with the Cox model
print("\n=== TESTING COX PREDICTION ===")
print("Adding survival variables for prediction...")
test_data$status <- NA
test_data$survival_time <- 1.5

print("Attempting prediction...")
tryCatch({
  probabilities <- predict(cox_model, test_data, type = "survival")
  print("SUCCESS! Prediction worked.")
  print("Predicted survival probabilities:")
  print(probabilities)
}, error = function(e) {
  print("ERROR in prediction:")
  print(e$message)
  print("Full error:")
  print(e)
})

# Try different prediction approaches
print("\n=== TRYING DIFFERENT PREDICTION APPROACHES ===")

# Approach 1: Using newdata parameter explicitly
tryCatch({
  probabilities <- predict(cox_model, newdata = test_data, type = "survival")
  print("SUCCESS with newdata parameter!")
  print(head(probabilities))
}, error = function(e) {
  print("ERROR with newdata parameter:")
  print(e$message)
})

# Approach 2: Check model frame
tryCatch({
  mf <- model.frame(cox_model, data = test_data)
  print("SUCCESS creating model frame!")
  print("Model frame dimensions:")
  print(dim(mf))
}, error = function(e) {
  print("ERROR creating model frame:")
  print(e$message)
})

print("\n=== DIAGNOSTIC COMPLETE ===")