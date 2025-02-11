# CHECKLIST BEFORE RUNNING SCRIPT 
# 1. Check #Writing to file section 
# 2. Check #covscan and #virscan sections -- comment out the one not being used
  # 2-1. Use get_cov_organism_by_id, cov_lib for covscan data and get_organism_by_id, vir_lib for virscan data
# 3. If using subset data, check #TEST ON A SUBSET section


library(dplyr)
library(ggplot2)
library(cowplot)
library(randomForest)
library(ROCR)
library(reshape2)
library(here)
library(data.table)
library(openxlsx)
library(randomForestExplainer)
library(viridis)

# # if library doesn't exist, install it
# if (!requireNamespace("rfPermute", quietly = TRUE)) install.packages("rfPermute", repos = "http://cran.rstudio.com", dependencies = TRUE)
# library(rfPermute) # p-value computation for random forest





# library(gt)
# library(Matrix)
# library(stringr)
library(purrr)
# library(broom)
# library(parallel) 
# library(gridExtra)  # For arranging multiple plots/tables
# library(ggrepel)    # For enhanced plot labeling
# library(grid) 
# library(pbapply)    # For parallel processing with progress bar
library(RColorBrewer)
# library(patchwork)

# anchor project root directory 
here::i_am("scripts/04_eda_randomforest.R")


# helper function
# and loading original data
# TODO: put in a separate script; is repeated in 03_misc_eda_tasks.R
# also put loading data stuff into a separate script 

# this helper function is for covscan; for now comment out when doing stuff with virscan data
get_cov_organism_by_id <- function(id_chr, cov_lib) {
  id_numeric <- as.numeric(id_chr)
  
  # Filter the row corresponding to the given id
  row_data <- cov_lib %>% 
    filter(id == id_numeric)
  
  # Check if row_data is empty
  if (nrow(row_data) == 0) {
    return(list(
      Organism = NA,
      Protein = NA,
      Sequence = NA,
      Start = NA,
      End = NA
    ))  # Return a list with NAs
  }
  
  # Extract attributes
  organism <- row_data %>% pull(Organism) %>% first()
  protein <- row_data %>% pull(Protein) %>% first()
  sequence <- row_data %>% pull(Sequence) %>% first()
  start_pos <- row_data %>% pull(start) %>% first()
  end_pos <- row_data %>% pull(end) %>% first()
  
  # Return attributes as a named list
  return(list(
    Organism = organism,
    Protein = protein,
    Sequence = sequence,
    Start = start_pos,
    End = end_pos
  ))
}

# this helper function is for virscan; for now comment out when doing stuff with covscan data
get_organism_by_id <- function(id_chr, vir_lib) {
  id_numeric <- as.numeric(id_chr)
  
  # Filter the row corresponding to the given id
  row_data <- vir_lib %>% 
    filter(id == id_numeric)
  
  # Check if row_data is empty
  if (nrow(row_data) == 0) {
    return(NA)  # Return NA if no matching id is found
  }
  
  # Try to get the organism
  label <- row_data %>%
    pull(Organism) %>%
    first()
  
  # Fallbacks if Organism is NA
  if (is.na(label)) {
    protein_name <- row_data %>%
      pull(`Protein names`) %>%
      first()
    if (!is.na(protein_name)) {
      label <- paste("protein:", protein_name)
    } else {
      species <- row_data %>%
        pull(Species) %>%
        first()
      if (!is.na(species)) {
        label <- paste("species:", species)
      } else {
        sequence <- row_data %>%
          pull(Sequence) %>%
          first()
        label <- paste("sequence:", sequence)
      }
    }
  }
  
  return(label)
}

#============================================================================

# load data

# covscan data
cov_z_pat <- readRDS(here::here("data", "processed", "cov_z_pat.rds"))
cov_lib <- readRDS(here::here("data", "processed", "cov_lib.rds"))


# virscan data
vir_z_pat <- readRDS(here::here("data", "processed", "vir_z_pat.rds"))
vir_lib <- readRDS(here::here("data", "rds", "vir_lib.rds"))


# #Writing to file ===========================================================
#TEST ON A SUBSET
# Define number of peptides for testing - default empty string
NUM_TEST_PEPTIDES <- ""

# CHECK BEFORE RUNNING - CHANGE TO APPROPRIATE NAMES ####################################################################
RF_RESULTS_FNAME <- "rf-result"
DATASET_NAME <- paste0("covscan", "_", NUM_TEST_PEPTIDES)
BASE_FNAME <- paste0(RF_RESULTS_FNAME, "_", DATASET_NAME)
EXCEL_FNAME <- paste0(BASE_FNAME, "_", DATASET_NAME, ".xlsx")
# 25-01-28-TODO make a directory for each run 
# CHECK BEFORE RUNNING - CHANGE TO APPROPRIATE NAMES ####################################################################
# Writing to file ===========================================================







#============================================================================
# TODO put data prep stuff in a separate script 
### data prep

# #covscan ####################################################################
META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "Hospitalized") # no COVID19_status

X_y_cov_z <- cov_z_pat %>%
  select(-all_of(META_COL_NO_STATUS)) %>%
  mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
  # put COVID19_status as the first column
  select(COVID19_status, everything()) %>% 
  mutate(COVID19_status = as.factor(COVID19_status)) # convert to factor

cat("Step: Data Preparation\n")
cat("Dimensions of X_y_cov_z:", dim(X_y_cov_z), "\n")
print(head(X_y_cov_z))  # Inspect a small portion

# #TEST ON A SUBSET 
# Randomly select peptides
# TEST_COL_FNAME <- paste0(BASE_FNAME, "_", "test_columns.rds")
# test_peptide_columns <- sample(colnames(X_y_cov_z[,-1]), NUM_TEST_PEPTIDES)
# saveRDS(test_peptide_columns, here::here("results", "eda", "rf", TEST_COL_FNAME))
# X_y_cov_z <- X_y_cov_z %>% select(all_of(c("COVID19_status", test_peptide_columns))) 

X_cov_z <- X_y_cov_z %>%
    select(-COVID19_status)

y_cov_z <- X_y_cov_z$COVID19_status

# Confirm dimensions
cat("Peptide Data Dimensions:", dim(X_cov_z), "\n")
cat("COVID19_status Dimensions:", length(y_cov_z), "\n")
# #covscan ####################################################################



# # #virscan ####################################################################
# META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "IP_reagent", 
#                    "Hospitalized", "Pull_down_antibody", 
#                    "patient", "original_id_column") # excludes COVID19_status which will be a part of the final dataframe


# # X_y_vir_z includes z-scores and COVID19_status vs. X_vir_z includes z-scores only
# X_y_vir_z <- vir_z_pat %>%
#     select(-all_of(META_COL_NO_STATUS)) %>%
#     mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
#     # put COVID19_status as the first column
#     select(COVID19_status, everything()) %>% 
#     mutate(COVID19_status = as.factor(COVID19_status))

# cat("Step: Data Preparation\n")
# cat("Dimensions of X_y_vir_z:", dim(X_y_vir_z), "\n")
# print(head(X_y_vir_z))  # Inspect a small portion

# # Randomly select peptides #TEST ON A SUBSET - comment if using full dataset 
# TEST_COL_FNAME <- paste0(BASE_FNAME, "_", "test_columns.rds")
# test_peptide_columns <- sample(colnames(X_y_vir_z[,-1]), NUM_TEST_PEPTIDES)
# # 25-01-27-TODO: save test_peptide_columns to file -- accidentally overwrote wrong dataframe from another script
# saveRDS(test_peptide_columns, here::here("results", "eda", "rf", TEST_COL_FNAME))
# # # Create the subset dataframe 
# X_y_vir_z <- X_y_vir_z # %>% select(all_of(c("COVID19_status", test_peptide_columns))) # uncomment for testing on small subset

# X_vir_z <- X_y_vir_z %>%
#     select(-COVID19_status)


# y_vir_z <- X_y_vir_z$COVID19_status

# # Confirm dimensions
# cat("Peptide Data Dimensions:", dim(X_vir_z), "\n")
# cat("COVID19_status Dimensions:", length(y_vir_z), "\n")
# #virscan ####################################################################





# Set dataset
# X <- 
# y <- 


# ### random forest implementation 


set.seed(618)

# #covscan train test split ####################################################################
print("splitting covscan data")
train_idx <- sample(seq_len(nrow(X)), size = 0.8 * nrow(X))
X_train <- X[train_idx, ]
y_train <- y[train_idx]
X_test  <- X[-train_idx, ]
y_test  <- y[-train_idx]
cat("Step: Train-Test Split\n")
cat("Train Set Dimensions:", dim(X_train), "\n")
cat("Test Set Dimensions:", dim(X_test), "\n")
cat("Number of positive cases in training set:", sum(y_train == "1"), "\n")
cat("Number of positive cases in testing set:", sum(y_test == "1"), "\n")
# #covscan ####################################################################



# # #virscan train test split ####################################################################
# # Convert y_vir_z to a factor if it's not already
# y_vir_z <- as.factor(y_vir_z)
# print("splitting virscan data")
# train_idx <- sample(seq_len(nrow(X_vir_z)), size = 0.8 * nrow(X_vir_z))
# X_train <- X_vir_z[train_idx, ]
# y_train <- y_vir_z[train_idx]
# X_test  <- X_vir_z[-train_idx, ]
# y_test  <- y_vir_z[-train_idx]
# cat("Step: Train-Test Split\n")
# cat("Train Set Dimensions:", dim(X_train), "\n")
# cat("Test Set Dimensions:", dim(X_test), "\n")
# cat("Number of positive cases in training set:", sum(y_train == "1"), "\n")
# cat("Number of positive cases in testing set:", sum(y_test == "1"), "\n")
# # #virscan ####################################################################





TR_TE_SPLIT_FNAME <- paste0(BASE_FNAME, "train_test_data.rds")
saveRDS(list(X_train = X_train, y_train = y_train, X_test = X_test, y_test = y_test), here::here("results", "eda", "rf", TR_TE_SPLIT_FNAME))
print("fitting rf")
# Train the random forest model
rf_model <- randomForest(X_train, y_train, importance = TRUE)
cat("Step: Random Forest Training\n")
print(rf_model)  # Print summary of the model

# Get predictions
# get class label predictions
class_predictions <- predict(rf_model, X_test, type = "response")
# get predicted probabilities for COVID19_status == positive
prob_predictions <- predict(rf_model, X_test, type = "prob")[,2]
accuracy <- mean(class_predictions == y_test)

# save rf_model
RF_MODEL_FNAME <- paste0(BASE_FNAME, "_", "model-fit.rds")
print("saving rf model")
saveRDS(rf_model, here::here("results", "eda", "rf", RF_MODEL_FNAME))
print("saved rf model -- done")

# save predictions
RF_CLASS_PRED_FNAME <- paste0(BASE_FNAME, "_", "class-predictions.rds")
RF_PROB_PRED_FNAME <- paste0(BASE_FNAME, "_", "prob-predictions.rds")
RF_ACC_FNAME <- paste0(BASE_FNAME, "_", "accuracy.rds")
print("saving rf predictions")
saveRDS(class_predictions, here::here("results", "eda", "rf", RF_CLASS_PRED_FNAME))
saveRDS(prob_predictions, here::here("results", "eda", "rf", RF_PROB_PRED_FNAME))
saveRDS(accuracy, here::here("results", "eda", "rf", RF_ACC_FNAME))
print("saved rf predictions -- done")



# ============================================================================
# Post-model fitting 
# ============================================================================
test_peptide_columns <- readRDS(here::here("results", "eda", "rf", "test_peptide_columns.rds"))
# Load saved model
rf_model <- readRDS(here::here("results", "eda", "rf", RF_MODEL_FNAME))

# Explain the model and compute p-values
importance_frame <- measure_importance(rf_model)
importance_frame <- importance_frame %>%
    mutate(peptide = as.character(variable))

# Extract feature importance
importance_data <- as.data.frame(importance(rf_model))
importance_data$peptide <- rownames(importance_data)

# add p-values
importance_data <- importance_data %>%
    left_join(importance_frame %>% select(peptide, p_value), by = "peptide")

TOP_NUM <- min(20, nrow(importance_data))


# Function to get the first N words of a string
get_first_n_words <- function(text, n = 10) {
  words <- unlist(strsplit(text, "\\s+"))  # Split the text into words
  paste(head(words, n), collapse = " ")  # Combine the first N words
}

importance_data_topN <- importance_data %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice_head(n = TOP_NUM) %>%
  mutate(
    attributes = map(peptide, get_cov_organism_by_id, cov_lib),
    organism = map_chr(attributes, "Organism"),
    Protein = map_chr(attributes, "Protein"),
    Sequence = map_chr(attributes, "Sequence"),
    Start = map_dbl(attributes, "Start"),
    End = map_dbl(attributes, "End")
  ) %>%
  select(-attributes) %>% # Remove the intermediate list column
  filter(!is.na(organism)) %>%
  distinct(organism, .keep_all = TRUE) %>% 
  mutate(
  first_n_words = sapply(organism, get_first_n_words, n = 10),  # Get the first N words of `organism`
  first_char_protein = substr(Protein, 1, 1),                  # Get the first character of `Protein`
  start_end = paste0(Start, "-", End),                        # Concatenate Start and End with "-"
  PEP_FULL_LABEL = paste(first_n_words, first_char_protein, start_end, sep = ", ")  # Combine all components
) %>%
select(-first_n_words, -first_char_protein, -start_end)  # Optionally drop intermediate columns

# save as cvs
IMP_FNAME <- paste0(BASE_FNAME, "_", "importance_data_top", TOP_NUM, ".csv")
write.csv(importance_data_topN, here::here("results", "eda", "rf", IMP_FNAME))

# Debugging: Check resulting dataframe
print(head(importance_data_topN))
cat("Step: Feature Importance\n")
cat("Dimensions of importance_data:", dim(importance_data), "\n")
print(head(importance_data))  # Inspect top rows

cat("Step: Top Features Selection\n")
cat("Dimensions of importance_data_topN:", dim(importance_data_topN), "\n")
print(head(importance_data_topN))  # Inspect selected top features

cat("Step: Organism Mapping\n")
cat("Unique organisms in importance_data_topN:\n")
print(unique(importance_data_topN$organism))


# Assuming merged_data contains 'peptide', 'p_value', and 'MeanDecreaseGini'
# Filter and rank features by p-values
ranked_p_values <- importance_data_topN %>%
  filter(!is.na(p_value)) %>%
  filter(p_value < 0.05)%>%
  arrange(p_value) %>%
  distinct(organism, .keep_all = TRUE) %>%
  head(TOP_NUM)  # Top 20 features by p-value

RANKED_P_FNAME <- paste0(BASE_FNAME, "_", "ranked_p_values.csv")
# save as cvs
write.csv(ranked_p_values, here::here("results", "eda", "rf", RANKED_P_FNAME))
cat("Unique organisms in ranked_p_values:\n")
print(unique(ranked_p_values$organism))
cat("Number of distinct organisms:", n_distinct(ranked_p_values$organism), "\n")


# Generate the color palette based on the number of distinct organisms
importance_palette <- viridis(n_distinct(importance_data_topN$organism))
ranked_p_palette <- viridis(n_distinct(ranked_p_values$organism))

cat("Step: Ranked Features for Plotting\n")
cat("Dimensions of ranked_p_values:", dim(ranked_p_values), "\n")
print(ranked_p_values)  # Inspect ranked features

cat("Step: Palette Configuration\n")
cat("Number of distinct organisms:", n_distinct(ranked_p_values$organism), "\n")
cat("Number of colors in importance palette:", length(importance_palette), "\n")
cat("Number of colors in ranked p palette:", length(ranked_p_palette), "\n")



# Plot feature importance with color scheme
ggplot(importance_data_topN, 
  aes(x = reorder(paste(peptide), MeanDecreaseGini), 
      y = MeanDecreaseGini, 
      fill = PEP_FULL_LABEL)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = importance_palette, name = "Organism") +
  theme_minimal() +
  labs(
    title = "Random Forest - Feature Importance",
    x = "Features (Organism)",
    y = "Mean Decrease Gini"
  ) +
  theme(
    legend.position = "right",  # Place the legend on the right
    axis.text.y = element_text(size = 8),  # Adjust size for long feature names
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)  # Center and bold title
  ) 
IMPORTANCE_FNAME <- paste0(BASE_FNAME, "_", "feature-importance.png")
ggsave(here::here("results", "eda", "rf", IMPORTANCE_FNAME), width = 10, height = 6, dpi = 300)

# Create the plot with the same color scheme
ggplot(ranked_p_values, 
       aes(x = reorder(peptide, -p_value), 
           y = -log10(p_value), 
           fill = PEP_FULL_LABEL)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = ranked_p_palette, name = "Organism") +
  labs(
    title = "Random Forest: Top Features by -log10(p-value)",
    x = "Features (Peptides)",
    y = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Place the legend on the right
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Centered and bold title
    plot.margin = unit(c(1.5, 1, 1, 1), "cm")  # Adjust margins
  )
PVAL_FNAME <- paste0(BASE_FNAME, "_", "top-peptides-by-pval.png")
ggsave(here::here("results", "eda", "rf", PVAL_FNAME), width = 10, height = 6, dpi = 300)

# load predictions from file and calculate AUC
test_predictions <- readRDS(here::here("results", "eda", "rf", RF_PROB_PRED_FNAME))

pred_test <- prediction(test_predictions, as.numeric(y_test) - 1)
perf_test <- performance(pred_test, "tpr", "fpr")
auc_test  <- performance(pred_test, "auc")@y.values[[1]]
# Save pred, perf and auc to a file
saveRDS(list(pred = pred_test, perf = perf_test, auc = auc_test), here::here("results", "eda", "rf", "rf_pred_perf_auc_subset.rds"))
# Plot ROC curve
# Create a data frame for the ROC curve
roc_data <- data.frame(
  fpr = unlist(perf_test@x.values),
  tpr = unlist(perf_test@y.values)
)

# Plot ROC curve using ggplot2
ggplot(roc_data, aes(x = fpr, y = tpr)) +
  geom_line(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = "ROC Curve",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  annotate("text", x = 0.75, y = 0.25, label = sprintf("AUC = %.2f", auc_test), color = "blue") +
  theme_minimal()

# Save the ROC plot
ROC_FNAME <- paste0(BASE_FNAME, "_", "roc_curve.png")
ggsave(here::here("results", "eda", "rf", ROC_FNAME), width = 10, height = 6, dpi = 300)
print("saved roc curve -- done!")















#============================================================================ 
# permutation-based importance testing RF
#============================================================================

# Fit the random forest model with permutation-based importance testing
# rf_perm_model <- rfPermute(X_vir_z, y_vir_z, importance = TRUE, ntree = 500)
# # save rf_perm_model
# saveRDS(rf_perm_model, here::here("results", "eda", "rf", "rf_perm_model_subset.rds"))
# print("saved rf_perm_model -- done")


# # Load saved model
# rf_perm_model <- readRDS(here::here("results", "eda", "rf", "rf_perm_model_subset.rds"))
# # Obtain performance metrics
# predictions <- predict(rf_perm_model, X_vir_z, type = "prob")[, 2]  # Probabilities for positive class
# # Save predictions
# saveRDS(predictions, here::here("results", "eda", "rf", "rf_perm_predictions.rds"))
# print("saved predictions -- done")
# # rf_model <- readRDS(here::here("results", "eda", "rf", "rf_model_subset.rds"))
# # Retrieve importance results
# importance_data <- importance(rf_perm_model, scale = FALSE)

# # Retrieve p-values for feature importance
# importance_pvalues <- rf_perm_model$pval

# # Combine importance and p-values
# importance_results <- data.frame(
#   Feature = rownames(importance_data),
#   Importance = importance_data[, "MeanDecreaseGini"],
#   PValue = importance_pvalues[, "MeanDecreaseGini"]
# )
# # View top features by importance
# head(importance_results)
# saveRDS(importance_results, here::here("results", "eda", "rf", "rf_perm_importance_results.rds"))
# print("saved importance_results -- done")


# # Load predictions
# predictions <- readRDS(file = "rf_perm_predictions.rds")

# # Create a prediction object
# pred <- prediction(predictions, as.numeric(y_vir_z) - 1)
# # Compute ROC performance
# perf <- performance(pred, "tpr", "fpr")
# # Compute AUC
# auc <- performance(pred, "auc")@y.values[[1]]

# # plot ROC curve
# plot(perf, col = "blue", main = "ROC Curve")
# abline(a = 0, b = 1, lty = 2, col = "gray")
# legend("bottomright", legend = sprintf("AUC = %.2f", auc), col = "blue", lwd = 2)





