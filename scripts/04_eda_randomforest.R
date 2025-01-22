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
vir_z_pat <- readRDS(here::here("data", "processed", "vir_z_pat.rds"))
vir_lib <- readRDS(here::here("data", "rds", "vir_lib.rds"))
# vir_z <- readRDS(here::here("data", "rds", "vir_z.rds"))

#============================================================================
# TODO put data prep stuff in a separate script 
### data prep
META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "IP_reagent", 
                   "Hospitalized", "Pull_down_antibody", 
                   "patient", "original_id_column") # excludes COVID19_status which will be a part of the final dataframe


# X_y_vir_z includes z-scores and COVID19_status vs. X_vir_z includes z-scores only
X_y_vir_z <- vir_z_pat %>%
    select(-all_of(META_COL_NO_STATUS)) %>%
    mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
    # put COVID19_status as the first column
    select(COVID19_status, everything()) %>% #factorize covid status
    mutate(COVID19_status = as.factor(COVID19_status))

# TEST ON A SUBSET
# Define number of peptides for testing
num_test_peptides <- 1000

# Randomly select peptides
test_peptide_columns <- sample(colnames(X_y_vir_z[,-1]), num_test_peptides)

# Create the subset dataframe
X_y_vir_z <- X_y_vir_z %>% select(all_of(c("COVID19_status", test_peptide_columns))) # uncomment for testing on small subset

X_vir_z <- X_y_vir_z %>%
    select(-COVID19_status)


y_vir_z <- X_y_vir_z$COVID19_status

# Confirm dimensions
cat("Peptide Data Dimensions:", dim(X_vir_z), "\n")

### random forest implementation 

set.seed(123)
# Convert y_vir_z to a factor if it's not already
y_vir_z <- as.factor(y_vir_z)

# print("fitting rf")
# # Train the random forest model
# rf_model <- randomForest(X_vir_z, y_vir_z, importance = TRUE)
# # save rf_model
# print("saving rf model")
# saveRDS(rf_model, here::here("results", "eda", "rf", "rf_model_subset.rds"))
# print("saved rf model -- done")

# Load saved model
rf_model <- readRDS(here::here("results", "eda", "rf", "rf_model_subset.rds"))

# Explain the model and compute p-values
importance_frame <- measure_importance(rf_model)
importance_frame <- importance_frame %>%
    mutate(peptide = as.character(variable))
# print(importance_frame)
# saveRDS(importance_frame, here::here("results", "eda", "rf", "rf_importance_frame_subset.rds"))

# Extract feature importance
importance_data <- as.data.frame(importance(rf_model))
importance_data$peptide <- rownames(importance_data)

# add p-values
importance_data <- importance_data %>%
    left_join(importance_frame %>% select(peptide, p_value), by = "peptide")

importance_data_top20 <- importance_data %>%
    arrange(desc(MeanDecreaseGini)) %>%
    slice_head(n = 20) %>%
    mutate(organism = map_chr(peptide, get_organism_by_id, vir_lib)) %>%
    # if any of the organism values is NA, drop the rows
    filter(!is.na(organism))

# Generate the color palette based on the number of distinct organisms
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(n_distinct(ranked_p_values$organism))
# Plot feature importance with color scheme
ggplot(importance_data_top20, 
       aes(x = reorder(paste(peptide), MeanDecreaseGini), 
           y = MeanDecreaseGini, 
           fill = organism)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = color_palette, name = "Organism") +
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
ggsave(here::here("results", "eda", "rf", "rf_feature_importance_subset_test.png"), width = 10, height = 6, dpi = 300)

# status: working
# ggplot(importance_data_top20, 
#     aes(x = reorder(paste(peptide, organism, sep = " - "), MeanDecreaseGini), 
#         y = MeanDecreaseGini)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +
#   theme_minimal() +
#   labs(title = "Feature Importance",
#     x = "Features (Organism)",
#     y = "Mean Decrease Gini")
# save plot
ggsave(here::here("results", "eda", "rf", "rf_feature_importance_subset_test.png"), width = 10, height = 6, dpi = 300)


# Assuming merged_data contains 'peptide', 'p_value', and 'MeanDecreaseGini'
# Filter and rank features by p-values
ranked_p_values <- importance_data_top20 %>%
  filter(!is.na(p_value)) %>%
  filter(p_value < 0.05)%>%
  arrange(p_value) %>%
  head(20)  # Top 20 features by p-value


# Create the plot with the same color scheme
ggplot(ranked_p_values, 
       aes(x = reorder(peptide, -p_value), 
           y = -log10(p_value), 
           fill = organism)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = color_palette, name = "Organism") +
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
ggsave(here::here("results", "eda", "rf", "rf_top_peptides_by_pval.png"), width = 10, height = 6, dpi = 300)
# # Plot ranked p-values on -log scale - status: works 
# ggplot(ranked_p_values, aes(x = reorder(peptide, -p_value), y = -log10(p_value))) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   coord_flip() +
#   theme_minimal() +
#   labs(
#     title = "Random Forest - Top Features by -log10(p-value)",
#     x = "Features (Peptides)",
#     y = "-log10(p-value)"
#   ) +
#   theme(axis.text.y = element_text(size = 8))

ggsave(here::here("results", "eda", "rf", "rf_feature_importance.png"), width = 10, height = 6, dpi = 300)


# # Obtain predictions and calculate AUC
# predictions <- predict(rf_model, X_vir_z, type = "prob")[,2]
# pred <- prediction(predictions, as.numeric(y_vir_z) - 1)
# perf <- performance(pred, "tpr", "fpr")
# auc <- performance(pred, "auc")@y.values[[1]]
# # Save pred, perf and auc to a file
# saveRDS(list(pred = pred, perf = perf, auc = auc), here::here("results", "eda", "rf", "rf_pred_perf_auc_subset.rds"))

# # Plot ROC curve
# plot(perf, col = "blue", main = "ROC Curve")
# abline(a = 0, b = 1, lty = 2, col = "gray")
# legend("bottomright", legend = sprintf("AUC = %.2f", auc), col = "blue", lwd = 2)

#============================================================================ 
# permutation-based importance testing RF
#============================================================================

# Fit the random forest model with permutation-based importance testing
rf_perm_model <- rfPermute(X_vir_z, y_vir_z, importance = TRUE, ntree = 500)
# save rf_perm_model
saveRDS(rf_perm_model, here::here("results", "eda", "rf", "rf_perm_model_subset.rds"))
print("saved rf_perm_model -- done")

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





