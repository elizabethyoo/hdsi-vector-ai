if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr,
  ggplot2,
  cowplot,
  randomForest,
  ROCR,
  reshape2,
  here,
  data.table,
  openxlsx,
  purrr
)

# # if library doesn't exist, install it
# if (!requireNamespace("rfPermute", quietly = TRUE)) install.packages("rfPermute", repos = "http://cran.rstudio.com", dependencies = TRUE)
# library(rfPermute) # p-value computation for random forest

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

lasso_results <- readRDS("/n/home01/egraff/hdsi-vector-ai/results/lasso_full_result.rds")

head(lasso_results)
head(lasso_results$coefficients)

dim(lasso_results$coefficients)

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


# anchor project root directory 
here::i_am("scripts/04_eda_randomforest.R")

# #Writing to file ===========================================================
#TEST ON A SUBSET
# Define number of peptides for testing - default empty string
NUM_TEST_PEPTIDES <- ""

# CHECK BEFORE RUNNING - CHANGE TO APPROPRIATE NAMES ####################################################################
RF_RESULTS_FNAME <- "rf-result"
DATASET_NAME <- paste0("covscan", "_", NUM_TEST_PEPTIDES)
BASE_FNAME <- paste0(RF_RESULTS_FNAME, "_", DATASET_NAME)
EXCEL_FNAME <- paste0(BASE_FNAME, "_", DATASET_NAME, ".xlsx")

# Create directory for this run
RUN_DIR <- here::here("results", "eda", "rf", paste0(BASE_FNAME, "_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))
dir.create(RUN_DIR, recursive = TRUE)

# covscan data============================================================================
# load data
# covscan data
cov_z_pat <- readRDS(here::here("data", "processed", "cov_z_pat.rds"))
cov_lib <- readRDS(here::here("data", "processed", "cov_lib.rds"))

# Prepare data
META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "Hospitalized")

X_y_cov_z <- cov_z_pat %>%
  select(-all_of(META_COL_NO_STATUS)) %>%
  mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
  select(COVID19_status, everything()) %>% 
  mutate(COVID19_status = as.factor(COVID19_status))

# TEST ON A SUBSET ==========================================================
if (NUM_TEST_PEPTIDES != "") {
  TEST_COL_FNAME <- paste0(BASE_FNAME, "_", "test_columns.rds")
  test_peptide_columns <- sample(colnames(X_y_cov_z[,-1]), as.numeric(NUM_TEST_PEPTIDES))
  saveRDS(test_peptide_columns, here::here("results", "eda", "rf", TEST_COL_FNAME))
  X_y_cov_z <- X_y_cov_z %>% select(all_of(c("COVID19_status", test_peptide_columns)))
}

X_cov_z <- X_y_cov_z %>% select(-COVID19_status)
y_cov_z <- X_y_cov_z$COVID19_status

# Perform 5-fold cross-validation
set.seed(618)
#============================================================================

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

# Extract feature importance
importance_data <- as.data.frame(importance(rf_model))
importance_data$peptide <- rownames(importance_data)
# add organism names

importance_data$Organism <- sapply(importance_data$Feature, get_organism_by_id, vir_lib = vir_lib)

importance_data_top20 <- importance_data %>%
    arrange(desc(MeanDecreaseGini)) %>%
    slice_head(n = 20) %>%
    mutate(organism = map_chr(peptide, get_organism_by_id, vir_lib)) %>%
    # if any of the organism values is NA, drop the rows
    filter(!is.na(organism))

ggplot(importance_data_top20, 
    aes(x = reorder(paste(peptide, organism, sep = " - "), MeanDecreaseGini), 
        y = MeanDecreaseGini)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Feature Importance",
    x = "Features (Organism)",
    y = "Mean Decrease Gini")
# save plot
ggsave(here::here("results", "eda", "rf", "rf_feature_importance_subset_test.png"), width = 10, height = 6, dpi = 300)

# Obtain predictions and calculate AUC
predictions <- predict(rf_model, X_vir_z, type = "prob")[,2]
pred <- prediction(predictions, as.numeric(y_vir_z) - 1)
perf <- performance(pred, "tpr", "fpr")
auc <- performance(pred, "auc")@y.values[[1]]
# Save pred, perf and auc to a file
saveRDS(list(pred = pred, perf = perf, auc = auc), here::here("results", "eda", "rf", "rf_pred_perf_auc_subset.rds"))

# Plot ROC curve
plot(perf, col = "blue", main = "ROC Curve")
abline(a = 0, b = 1, lty = 2, col = "gray")
legend("bottomright", legend = sprintf("AUC = %.2f", auc), col = "blue", lwd = 2)

#============================================================================ 
# permutation-based importance testing RF
#============================================================================

# Fit the random forest model with permutation-based importance testing
rf_perm_model <- rfPermute(X_vir_z, y_vir_z, importance = TRUE, ntree = 500)
# save rf_perm_model
saveRDS(rf_perm_model, here::here("results", "eda", "rf", "rf_perm_model_subset.rds"))
print("saved rf_perm_model -- done")

# # Load saved model
rf_perm_model <- readRDS(here::here("results", "eda", "rf", "rf_perm_model_subset.rds"))
# Obtain performance metrics
predictions <- predict(rf_perm_model, X_vir_z, type = "prob")[, 2]  # Probabilities for positive class
# Save predictions
saveRDS(predictions, here::here("results", "eda", "rf", "rf_perm_predictions.rds"))
print("saved predictions -- done")
# rf_model <- readRDS(here::here("results", "eda", "rf", "rf_model_subset.rds"))
# Retrieve importance results
importance_data <- importance(rf_perm_model, scale = FALSE)

# Retrieve p-values for feature importance
importance_pvalues <- rf_perm_model$pval

# Combine importance and p-values
importance_results <- data.frame(
  Feature = rownames(importance_data),
  Importance = importance_data[, "MeanDecreaseGini"],
  PValue = importance_pvalues[, "MeanDecreaseGini"]
)
# View top features by importance
head(importance_results)
saveRDS(importance_results, here::here("results", "eda", "rf", "rf_perm_importance_results.rds"))
print("saved importance_results -- done")


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


# 25-01-28 rf on covidscan data scratchwork
model_fit_test <- readRDS("/n/home01/egraff/hdsi-vector-ai/results/eda/rf/rf-result_covscan__model-fit.rds")
prob_pred_test <- readRDS("/n/home01/egraff/hdsi-vector-ai/results/eda/rf/rf-result_covscan__prob-predictions.rds")

# 25-02-11 combine virscan and covscan data 
# CHECKLIST BEFORE RUNNING SCRIPT 
# 1. Check #Writing to file section 
# 2. Check #covscan and #virscan sections -- comment out the one not being used
  # 2-1. Use get_cov_organism_by_id, cov_lib for covscan data and get_organism_by_id, vir_lib for virscan data
# 3. If using subset data, check #TEST ON A SUBSET section


if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
dplyr,
ggplot2,
reshape2,
here,
gt,
stringr,
data.table,
openxlsx,
pbapply,
grplasso,
parallel
)

# anchor project root directory 
here::i_am("scripts/05_eda_grplasso.R")


# #Writing to file ===========================================================
#TEST ON A SUBSET
# Define number of peptides for testing - default empty string
NUM_TEST_PEPTIDES <- ""

# CHECK BEFORE RUNNING - CHANGE TO APPROPRIATE NAMES ####################################################################
GRPLASSO_RESULTS_FNAME <- "grplasso-result"
DATASET_NAME <- paste0("covscan", "_", NUM_TEST_PEPTIDES)
BASE_FNAME <- paste0(GRPLASSO_RESULTS_FNAME, "_", DATASET_NAME)
EXCEL_FNAME <- paste0(BASE_FNAME, "_", DATASET_NAME, ".xlsx")

# Create directory for this run
RUN_DIR <- here::here("results", "eda", "grplasso", paste0(BASE_FNAME, "_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))
dir.create(RUN_DIR, recursive = TRUE)
# CHECK BEFORE RUNNING - CHANGE TO APPROPRIATE NAMES ####################################################################
# Writing to file ===========================================================

# Get number of cores from SLURM, fallback to detectCores() - 1 if not set
num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(num_cores) || num_cores <= 0) {
    num_cores <- parallel::detectCores() - 1
}

cat("Using", num_cores, "cores for parallel processing\n")

# TODO - put in a separate script; is repeated in 03_misc_eda_tasks.R
# also put loading data stuff into a separate script
# # import helper functions script
# source(here("scripts", "00_helper_funcs.R"))

DEBUG_LOG <- "grplasso_debug_log.txt"

# load data
vir_z_pat <- readRDS(here::here("data", "processed", "vir_z_pat.rds"))
vir_lib <- readRDS(here::here("data", "rds", "vir_lib.rds"))
vir_z <- readRDS(here::here("data", "rds", "vir_z.rds"))


# covscan data============================================================================
# load data
# covscan data
cov_z_pat <- readRDS(here::here("data", "processed", "cov_z_pat.rds"))
cov_lib <- readRDS(here::here("data", "processed", "cov_lib.rds"))

# Prepare data
META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "Hospitalized")

X_y_cov_z <- cov_z_pat %>%
  select(-all_of(META_COL_NO_STATUS)) %>%
  mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
  select(COVID19_status, everything()) %>% 
  mutate(COVID19_status = as.factor(COVID19_status))

# TEST ON A SUBSET ==========================================================
seed.set(123)
if (NUM_TEST_PEPTIDES != "") {
  TEST_COL_FNAME <- paste0(BASE_FNAME, "_", "test_columns.rds")
  test_peptide_columns <- sample(colnames(X_y_cov_z[,-1]), as.numeric(NUM_TEST_PEPTIDES))
  saveRDS(test_peptide_columns, here::here("results", "eda", "rf", TEST_COL_FNAME))
  X_y_cov_z <- X_y_cov_z %>% select(all_of(c("COVID19_status", test_peptide_columns)))
}

X_cov_z <- X_y_cov_z %>% select(-COVID19_status)
y_cov_z <- X_y_cov_z$COVID19_status

X <- X_cov_z %>%
  as.matrix()

y <- y_cov_z %>%
  as.numeric()

cat("Peptide Data Dimensions:", dim(X_cov_z), "\n")
cat("Response Vector Length:", length(y_cov_z), "\n")

# modify get_organism_from_id to account for covid_lib
# TODO change get_organism_from_id to get_vir_organism_by_id
# Need to extract: Protein, start, end, (peptide) sequence
# also cannot query on the headnode for rhinovirus A, Human Herpesvirus4 and HIV-1
get_cov_organism_by_id <- function(id_chr, cov_lib) {
  id_numeric <- as.numeric(id_chr)
  
  # Filter the row corresponding to the given id
  row_data <- cov_lib %>% 
    filter(id == id_numeric)
  
  # Check if row_data is empty
  if (nrow(row_data) == 0) {
    return(NA)  # Return NA if no matching id is found
  }
  # Try to get the Organism, Protein, Sequence, start, and end
  label <- row_data %>%
    pull(Organism) %>%
    first()
  
  protein <- row_data %>%
    pull(Protein) %>%
    first()
  
  sequence <- row_data %>%
    pull(Sequence) %>%
    first()
  
  start_pos <- row_data %>%
    pull(start) %>%
    first()
  
  end_pos <- row_data %>%
    pull(end) %>%
    first()
  
  # Combine the extracted information into a single string
  label <- paste("Organism:", label, 
                 "Protein:", protein, 
                 "Sequence:", sequence, 
                 "Start:", start_pos, 
                 "End:", end_pos)
  
  return(label)
}


#' get_organism_by_id
#'
#' This helper function returns descriptors for one or more peptide IDs based on data 
#' in the provided vir_lib data frame. The priority for the label is:
#' 1) Organism
#' 2) Protein names (prefixed "protein:")
#' 3) Species (prefixed "species:")
#' 4) Sequence (prefixed "sequence:")
#' If no row is found for an ID, NA is returned.
#'
#' @param id_chr A character vector (or single value) of peptide IDs.
#' @param vir_lib A data frame with fields such as "id", "Organism", "Protein names", 
#'   "Species", and "Sequence".
#' @return A character vector of labels or NA where not found.
#'
#' @examples
#' # Assuming vir_lib is a data frame with the necessary columns
#' get_organism_by_id(c("1234", "9999"), vir_lib)
#'
#' @export
get_organism_by_id <- function(id_chr, vir_lib) {
    ids_numeric <- suppressWarnings(as.numeric(id_chr))
    
    sapply(ids_numeric, function(single_id) {
        row_data <- vir_lib %>%
            filter(id == single_id)
        
        if (nrow(row_data) == 0) {
            return(NA_character_)
        }
        
        label <- row_data %>%
            pull(Organism) %>%
            first()
        
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
        label
    })
}

#============================================================================
# virscan data============================================================================
### data prep
META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "IP_reagent", 
                   "Hospitalized", "Pull_down_antibody", 
                   "patient", "original_id_column") # excludes COVID19_status which will be a part of the final dataframe


# X_y_vir_z includes z-scores and COVID19_status vs. X_vir_z includes z-scores only
X_y_vir_z <- vir_z_pat %>%
    select(-all_of(META_COL_NO_STATUS)) %>%
    mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
    # put COVID19_status as the first column
    select(COVID19_status, everything())

# TEST ON A SUBSET ==========================================================
seed.set(123)
if (NUM_TEST_PEPTIDES != "") {
  TEST_COL_FNAME <- paste0(BASE_FNAME, "_", "test_columns.rds")
  test_peptide_columns <- sample(colnames(X_y_vir_z[,-1]), as.numeric(NUM_TEST_PEPTIDES))
  saveRDS(test_peptide_columns, here::here("results", "eda", "rf", TEST_COL_FNAME))
  X_y_vir_z <- X_y_vir_z %>% select(all_of(c("COVID19_status", test_peptide_columns)))
}

X_vir_z <- X_y_vir_z %>% select(-COVID19_status)
y_vir_z <- X_y_vir_z$COVID19_status

X <- X_vir_z %>%
  as.matrix()

y <- y_vir_z %>%
  as.numeric()

# Confirm dimensions
sink(DEBUG_LOG, append = TRUE)
cat("Peptide Data Dimensions:", dim(X_vir_z), "\n")
cat("Response Vector Length:", length(y_vir_z), "\n")
