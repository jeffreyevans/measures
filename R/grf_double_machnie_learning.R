library(grf)
library(tidyr)

# Load your time series data with treatment and covariates
# Assume you have a data frame "df" with columns: 
#   outcome, time, treatment, covariate1, covariate...

# example of wide to long
wide_data <- data.frame(
  Date = seq(as.Date("2023-01-01"), by = "months", length.out = 3),
  Value_1 = c(10, 20, 30),
  Value_2 = c(15, 25, 35),
  Value_3 = c(12, 22, 32)
)

# Convert wide format to long format
treatment_control <- pivot_longer(wide_data, cols = starts_with("Value"), 
                          names_to = "Variable", 
						  values_to = "Value")

# Split the data into treatment and control groups
treatment_data <- df[df$treatment == 1, ]
control_data <- df[df$treatment == 0, ]

# Create DML treatment and control models 
dml_treatment <- causal_forest(X = treatment_data[, c("covariate1", "covariate2")],
                               Y = treatment_data$outcome,
                               W = treatment_data$treatment)

dml_control <- causal_forest(X = control_data[, c("covariate1", "covariate2")],
                             Y = control_data$outcome,
                             W = control_data$treatment)

# Estimate treatment effect tau_hat [difference between estimated treatments and controls]
tau_hat <- predict(dml_treatment, X = df[, c("covariate1", "covariate2")]) -
             predict(dml_control, X = df[, c("covariate1", "covariate2")])
