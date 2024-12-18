#!/usr/bin/env Rscript

# DESCRIPTION:
# This R script performs LASSO classification on selected differentially expressed genes (DEGs).
# It reads upregulated and downregulated gene lists, normalizes gene expression data,
# applies LASSO regression with Leave-One-Out Cross-Validation (LOOCV), and calculates various evaluation metrics.
# The script also generates a ROC curve and exports evaluation results, including selected features, to a PDF.

# USAGE:
# Rscript 05_runLASSO_DEG.R <upregulated_genes_file> <downregulated_genes_file> <normalized_expression_file>

# Load necessary libraries
library(glmnet)
library(readxl)
library(dplyr)
library(pROC)
library(mltools)
library(gridExtra)

# Function to read gene lists
read_gene_lists <- function(upregulated_file, downregulated_file) {
  DEG_up <- read_excel(upregulated_file)      # Load upregulated genes
  DEG_down <- read_excel(downregulated_file)  # Load downregulated genes
  return(rbind(DEG_up, DEG_down))
}

# Function to preprocess expression data
preprocess_expression_data <- function(normalized_expression_file, DEG) {
  info_expr <- read.delim(normalized_expression_file)
  i_DEG <- which(rownames(info_expr) %in% DEG$ID)  # Identify rows in the expression matrix
  miRNAs_selected <- info_expr[i_DEG, ]             # Subset the expression data
  return(t(miRNAs_selected))                          # Transpose for further analysis
}

# Function to prepare target variable
prepare_target_variable <- function(snames_2) {
  split_names <- strsplit(snames_2, "_")            # Split names by "_"
  sample_group <- sapply(split_names, function(x) x[3]) # Extract the group information
  sample_group[is.na(sample_group)] <- "C"           # Replace NA with control group
  sample_group <- lapply(sample_group, function(x) ifelse(x == "B" | x == "BS", 1, 0)) # B = 1, V = 0
  return(unlist(sample_group))                        # Flatten the list to vector
}

# Function to perform LASSO classification
perform_lasso_classification <- function(X, y) {
  set.seed(123)  # Set seed for reproducibility
  lasso_model <- cv.glmnet(X, y, family = "binomial", type.measure = "class", nfolds = nrow(X))
  return(lasso_model)
}

# Function to compute evaluation metrics
compute_metrics <- function(y, y_pred, y_pred_prob) {
  conf_matrix <- table(Actual = y, Predicted = y_pred)
  
  # Calculate accuracy
  accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
  
  # Compute ROC curve and AUC
  roc_data <- roc(y, y_pred_prob)
  auc_score <- auc(roc_data)
  
  # Compute True Positives (TP), False Positives (FP), False Negatives (FN), True Negatives (TN)
  TP <- sum(y_pred == 1 & y == 1)
  FP <- sum(y_pred == 1 & y == 0)
  FN <- sum(y_pred == 0 & y == 1)
  TN <- sum(y_pred == 0 & y == 0)
  
  # Calculate Positive Predictive Value (PPV) and other metrics
  PPV <- TP / (TP + FP)
  
  # Calculate Matthews correlation coefficient (MCC)
  MCC <- mcc(TP = TP, TN = TN, FN = FN, FP = FP)
  
  # Calculate precision, recall, and F1 score
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  F1_score <- 2 * (precision * recall) / (precision + recall)
  
  # Create a data frame for evaluation metrics
  evaluation_metrics <- data.frame(
    Metric = c("Accuracy", "AUC", "PPV", "Precision", "Recall", "F1 Score", "MCC"),
    Value = c(round(accuracy, 3), round(auc_score, 3), round(PPV, 3),
              round(precision, 3), round(recall, 3), round(F1_score, 3), round(MCC, 3))
  )
  
  return(list(evaluation_metrics = evaluation_metrics, roc_data = roc_data))
}

# Function to plot ROC and export results to PDF
export_results_to_pdf <- function(roc_data, evaluation_metrics, selected_features_df) {
  pdf("evaluation_results.pdf", width = 10, height = 10)
  
  # Plot the ROC curve
  plot(roc_data, main = "ROC Curve", col = "blue", legacy.axes = TRUE, xlab = '1 - Specificity', ylab='Sensitivity')
  text(x = 0.4, y = 0.2, paste('AUC score:', round(auc(roc_data), 3)))  # Print AUC score on the plot
  
  # Add evaluation metrics table to the PDF
  gridExtra::grid.table(evaluation_metrics, rows = NULL, theme = ttheme_minimal())
  
  # Add selected features table to the PDF
  gridExtra::grid.table(selected_features_df, rows = NULL, theme = ttheme_minimal())
  
  # Close the PDF device
  dev.off()
}

# Main script execution
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Please provide three arguments: <upregulated_genes_file> <downregulated_genes_file> <normalized_expression_file>")
}

# Read command-line arguments
upregulated_genes_file <- args[1]
downregulated_genes_file <- args[2]
normalized_expression_file <- args[3]

# Read and process data
DEG <- read_gene_lists(upregulated_genes_file, downregulated_genes_file)
miRNAs_selected <- preprocess_expression_data(normalized_expression_file, DEG)

# Prepare the target variable
miRNAs_selected <- t(miRNAs_selected)  # Transpose for further analysis
snames_2 <- rownames(miRNAs_selected)
sample_group <- prepare_target_variable(snames_2)

# Prepare data matrix (X) and target vector (y)
X <- as.matrix(miRNAs_selected)  # Features matrix
y <- as.factor(sample_group)      # Target variable as factor

# Perform LASSO classification
lasso_model <- perform_lasso_classification(X, y)

# Get the optimal lambda value
best_lambda <- lasso_model$lambda.min

# Extract coefficients for the optimal lambda
lasso_coefficients <- coef(lasso_model, s = best_lambda)

# Extract coefficients for the positive class
class_coef <- as.matrix(lasso_coefficients)[-1, ]  # Remove intercept

# Combine coefficients and calculate their absolute values
original_feature_names <- colnames(X)
coefficients_df <- data.frame(feature = original_feature_names, abs_coef = abs(as.numeric(class_coef)))

# Sort features based on their absolute coefficient values
sorted_features <- coefficients_df[order(coefficients_df$abs_coef, decreasing = TRUE), ]

# Print the top features that were selected by the LASSO model
selected_features <- sorted_features[sorted_features$abs_coef != 0, "feature"]

# Filter data matrix to include only selected features
selected_features_df <- subset(coefficients_df, feature %in% selected_features)
X_selected <- X[, selected_features]  # Subset X for selected features

# Initialize prediction vectors
y_pred <- numeric(length(y))
y_pred_prob <- numeric(length(y))

# Implement Leave-One-Out Cross-Validation for predictions
for (i in 1:length(y)) {
  # Leave out the ith sample for testing
  X_train <- X_selected[-i, ]
  y_train <- y[-i]
  
  # Fit LASSO model
  lasso_model <- glmnet(X_train, y_train, family = "binomial")
  
  # Predict probabilities for the left-out sample
  y_pred_prob[i] <- predict(lasso_model, newx = X_selected[i, ], s = best_lambda, type = "response")
  
  # Assign class label based on probability threshold (0.5)
  y_pred[i] <- ifelse(y_pred_prob[i] > 0.5, 1, 0)
}

# Convert predicted labels to factor
y_pred <- as.factor(y_pred)

# Compute evaluation metrics
metrics_results <- compute_metrics(y, y_pred, y_pred_prob)

# Export results to PDF
export_results_to_pdf(metrics_results$roc_data, metrics_results$evaluation_metrics, selected_features_df)
