
# R Script: Model Comparison Tables for Five Datasets Using EZ and Competing Distributions

# Required packages
library(fitdistrplus)
library(goftest)
library(stats4)

# Placeholder function definitions for fitting models
fit_model <- function(data, model_type) {
  # model_type: "EZ", "EPI", "EE", "Z"
  # This is a placeholder. Actual implementation depends on model equations.
  return(list(
    W_star = runif(1, 1.1, 1.4),
    A_star = runif(1, 0.75, 0.95),
    AIC = runif(1, 115, 135),
    CAIC = runif(1, 116, 137),
    BIC = runif(1, 117, 138),
    HQIC = runif(1, 115, 137),
    KS = runif(1, 0.09, 0.14),
    p_value = runif(1, 0.05, 0.12)
  ))
}

# List of datasets (replace with actual data loading)
datasets <- list(
  "Glass Strength" = rnorm(50, mean = 5, sd = 1),
  "Air Conditioning Failures" = rexp(50, rate = 0.1),
  "Vinyl Chloride" = rgamma(50, shape = 2, rate = 1),
  "Electronic Components" = rexp(15, rate = 0.2),
  "Carbon Fibres" = rnorm(100, mean = 3, sd = 0.5)
)

models <- c("EPI", "EE", "Z", "EZ")

# Function to generate summary for one dataset
evaluate_dataset <- function(data) {
  results <- lapply(models, function(m) unlist(fit_model(data, m)))
  results_df <- as.data.frame(results)
  rownames(results_df) <- names(results[[1]])
  colnames(results_df) <- models
  return(results_df)
}

# Evaluate all datasets
results_list <- lapply(datasets, evaluate_dataset)

# Print results
for (name in names(results_list)) {
  cat("\nModel comparison for:", name, "\n")
  print(round(results_list[[name]], 3))
}
