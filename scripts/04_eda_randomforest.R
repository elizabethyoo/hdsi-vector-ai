library(dplyr)
library(ggplot2)
library(cowplot)
library(randomForest)
library(ROCR)
library(reshape2)
library(here)
library(data.table)
library(openxlsx)
# library(gt)
# library(Matrix)
# library(stringr)
# library(purrr)
# library(broom)
# library(parallel) 
# library(gridExtra)  # For arranging multiple plots/tables
# library(ggrepel)    # For enhanced plot labeling
# library(grid)
# library(pbapply)    # For parallel processing with progress bar
# library(RColorBrewer)
# library(patchwork)

# anchor project root directory 
here::i_am("scripts/04_eda_randomforest.R")

#============================================================================

# load data
vir_z_pat <- readRDS(here::here("data", "processed", "vir_z_pat.rds"))
# vir_lib <- readRDS(here::here("data", "rds", "vir_lib.rds"))
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
X_y_vir_z <- X_y_vir_z # %>% select(all_of(c("COVID19_status", test_peptide_columns))) # uncomment for testing on small subset

X_vir_z <- X_y_vir_z %>%
    select(-COVID19_status)


y_vir_z <- X_y_vir_z$COVID19_status

# Confirm dimensions
cat("Peptide Data Dimensions:", dim(X_vir_z), "\n")

### random forest implementation 

set.seed(123)
# Convert y_vir_z to a factor if it's not already
y_vir_z <- as.factor(y_vir_z)

print("fitting rf")
# Train the random forest model
rf_model <- randomForest(X_vir_z, y_vir_z, importance = TRUE)
# save rf_model
print("saving rf model")
saveRDS(rf_model, here::here("results", "eda", "rf", "rf_model_subset.rds"))
print("saved rf model -- done")

# Load saved model
# rf_model <- readRDS(here::here("results", "eda", "rf", "rf_model_subset.rds"))

# # Extract feature importance
# importance_data <- as.data.frame(importance(rf_model))
# importance_data$Feature <- rownames(importance_data)

# importance_data_top20 <- importance_data %>%
#     arrange(desc(MeanDecreaseGini)) %>%
#     slice_head(n = 20)

# # Plot feature importance (Mean Decrease Gini)
# ggplot(importance_data_top20, aes(x = reorder(Feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
#   geom_bar(stat = 'identity') +
#   coord_flip() +
#   theme_minimal() +
#   labs(title = "Feature Importance", x = "Features", y = "Mean Decrease Gini")


