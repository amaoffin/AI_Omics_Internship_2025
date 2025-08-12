#### Try It Yourself ####
# Practice Exercises 

# ----------------------------------------------------------------------------------------------------------------

# 1. Check Cholesterol level (using if) 
# Write an If statement to check cholesterol level is greater than 240, 
# if true, it will prints “High Cholesterol”

cholesterol <- 230
if (cholesterol > 240) {
  print("High Cholesterol")
}

# ----------------------------------------------------------------------------------------------------------------

# 2. Blood Pressure Status (using if...else)
# Write an if…else statement to check if blood pressure is normal.
# If it’s less than 120, print: “Blood Pressure is normal”
# If false then print: “Blood Pressure is high”

Systolic_bp <- 130
if (Systolic_bp < 120) {
  print("Blood Pressure is normal")
} else {
  print("Blood Pressure is high")
}

# ----------------------------------------------------------------------------------------------------------------

# 3. Automating Data Type Conversion with for loop

# Use patient_info.csv data and its metadata.
# Perform the following steps separately on each dataset (data and metadata).
# Create a copy of the dataset to work on.
# Identify all columns that should be converted to factor type.
# Store their names in a variable (factor_cols).

# Example: factor_cols <- c("gender", "smoking_status")

# Use a for loop to convert all the columns in factor_cols to factor type.
# Pass factor_cols to the loop as a vector.

# Hint:
# for (col in factor_cols) {
#   data[[col]] <- as.factor(data[[col]])  # Replace 'data' with the name of your dataset
# }

#getting datasets
raw_patient_data <- read.csv(file.choose())
patient_data <- raw_patient_data
raw_metadata <- read.csv(file.choose())
metadata <- raw_metadata

#identifying in patient data
factor_cols_patient <- c()
for (i in 1:ncol(patient_data)) {
  if (is.character(patient_data[[i]])) {
    factor_cols_patient <- c(factor_cols_patient, names(patient_data)[i])
  }
}
#identifying i metadata
factor_cols_metadata <- c()
for (i in 1:ncol(metadata)) {
  if (is.character(metadata[[i]])) {
    factor_cols_metadata <- c(factor_cols_metadata, names(metadata)[i])
  }
}

#converting patient data
for (col in factor_cols_patient) {
  patient_data[[col]] <- as.factor(patient_data[[col]])
}
#converting metadata
for (col in factor_cols_metadata) {
  metadata[[col]] <- as.factor(metadata[[col]])
}

str(patient_data)
str(metadata)

# ----------------------------------------------------------------------------------------------------------------

# 4. Converting Factors to Numeric Codes

# Choose one or more factor columns (e.g., smoking_status).
# Convert "Yes" to 1 and "No" to 0 using a for loop.

# Hint:
# binary_cols <- c("smoking_status")   # store column names in a vector
# use ifelse() condition inside the loop to replace Yes with 1 and No with 0
# for (col in binary_cols) {
#   data[[col]] <- # insert your ifelse() code here
# }

binary_cols <- c("smoker", "diagnosis")

for (col in binary_cols) {
  patient_data[[col]] <- ifelse(patient_data[[col]] == "Yes" | patient_data[[col]] == "Cancer", 1, 0)
}
# ----------------------------------------------------------------------------------------------------------------

# 3. Verification:
#    Compare the original and modified datasets to confirm changes.
str(raw_patient_data)
str(patient_data)
str(raw_metadata)
str(metadata)

# ----------------------------------------------------------------------------------------------------------------
