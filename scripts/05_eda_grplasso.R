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

# parallel processing ver

DEBUG_LOG <- "grplasso_debug_log.txt"

run_grplasso_reg_parallel <- function(X, Y, groups = NULL, split_type = c("train_validate_test", "train_test"), 
                                      train_ratio = 0.6, valid_ratio = 0.2, visualize = FALSE, 
                                      desired_number_of_groups = 5, num_cores) {
    split_type <- match.arg(split_type)
    n <- nrow(X)
    p <- ncol(X)

    sink(DEBUG_LOG, append = TRUE)
    cat("Shape of X:", dim(X), "\n")
    cat("Length of Y:", length(Y), "\n")
    cat("First few rows of X:\n")
    print(head(X[, 1:min(5, ncol(X))]))
    cat("First few values of Y:\n")
    print(head(Y))
    cat("Number of missing values in X:", sum(is.na(X)), "\n")
    cat("Number of missing values in Y:", sum(is.na(Y)), "\n")

    zero_var_cols <- which(apply(X, 2, var) == 0)
    if (length(zero_var_cols) > 0) {
        cat("Columns with zero variance:", zero_var_cols, "\n")
    }
    sink()

    if (is.null(groups)) {
        stop("Error: Groups not provided. Please specify a grouping option.")
    } else if (groups == "corr") {
        corr_matrix <- cor(X)
        dist_matrix <- as.dist(1 - abs(corr_matrix))
        hc <- hclust(dist_matrix, method = "average")
        groups <- cutree(hc, k = desired_number_of_groups)
    } else if (groups == "pca") {
        pca_result <- prcomp(X, scale. = TRUE)
        loading_matrix <- abs(pca_result$rotation[, 1:min(ncol(X), 10)])
        groups <- apply(loading_matrix, 1, which.max)
    } else if (groups == "org") {

    }
    else {
        groups <- rep(1:ncol(X), each = 1)
    }
    index_final <- c(NA, groups)

    sink(DEBUG_LOG, append = TRUE)
    cat("Number of groups provided:", length(groups), "\n")
    cat("Number of features in X:", ncol(X), "\n")
    sink()

    if (length(groups) != ncol(X)) {
        stop("Error: Mismatch between group length and number of features in X.")
    }

    set.seed(123)
    train_idx <- sample(1:n, size = train_ratio * n)
    if (split_type == "train_validate_test") {
        valid_idx <- setdiff(sample(1:n, size = (train_ratio + valid_ratio) * n), train_idx)
        test_idx <- setdiff(1:n, c(train_idx, valid_idx))
    } else {
        test_idx <- setdiff(1:n, train_idx)
    }

    X_train <- cbind(Intercept = 1, X[train_idx, ])
    Y_train <- Y[train_idx]
    X_test <- cbind(Intercept = 1, X[test_idx, ])
    Y_test <- Y[test_idx]

    sink(DEBUG_LOG, append = TRUE)
    cat("Train set size:", length(train_idx), "\n")
    if (split_type == "train_validate_test") {
        cat("Validation set size:", length(valid_idx), "\n")
    }
    cat("Test set size:", length(test_idx), "\n")
    cat("Dimensions of training data:", dim(X_train), "\n")
    cat("Dimensions of test data:", dim(X_test), "\n")
    sink()

    lambdas <- exp(seq(log(1e-4), log(10), length.out = 100))

    if (split_type == "train_validate_test") {
        X_valid <- cbind(Intercept = 1, X[valid_idx, ])
        Y_valid <- Y[valid_idx]

        cl <- makeCluster(num_cores)
        clusterExport(cl, varlist = c("X_train", "Y_train", "index_final", "X_valid", "Y_valid", "grplasso", "DEBUG_LOG"))
        clusterEvalQ(cl, library(grplasso))

        cv_results <- parSapply(cl, lambdas, function(lambda) {
            sink(DEBUG_LOG, append = TRUE)
            cat("Running cross-validation for lambda =", lambda, "\n")
            sink()

            fit <- grplasso(X_train, Y_train, index = index_final, lambda = lambda, model = LinReg())
            mse <- mean((Y_valid - predict(fit, X_valid))^2)

            sink(DEBUG_LOG, append = TRUE)
            cat("MSE for lambda =", lambda, ":", mse, "\n")
            sink()
            
            return(mse)
        })

        stopCluster(cl)
        best_lambda <- lambdas[which.min(cv_results)]
    } else {
        best_lambda <- lambdas[1]
    }

    final_model <- grplasso(X_train, Y_train, index = index_final, lambda = best_lambda, model = LinReg())
    Y_test_pred <- predict(final_model, X_test)
    test_mse <- mean((Y_test - Y_test_pred)^2)

    coefs <- as.numeric(coef(final_model))

    sink(DEBUG_LOG, append = TRUE)
    if (any(is.na(coefs))) {
        stop("NA values found in coefficients")
    }
    cat("First few coefficients:\n")
    print(head(coefs))

    coefs_no_intercept <- data.frame(
        feats = colnames(X_train),
        coef = coefs,
        group = index_final,
        is_intercept = is.na(index_final)
    ) %>% arrange(desc(abs(coef)))

    cat("Coefficient extraction successful. First few sorted coefficients:\n")
    print(head(coefs_no_intercept))

    cat("Best lambda selected:", best_lambda, "\n")
    cat("Final test MSE:", test_mse, "\n")
    sink()

    list(
        model = final_model,
        test_mse = test_mse,
        coefs = coefs_no_intercept,
        best_lambda = best_lambda
    )
}



# non parallel processing ver - is very slow

# run_grplasso_reg <- function(X, Y, groups = NULL, split_type = c("train_validate_test", "train_test"), 
#                              train_ratio = 0.6, valid_ratio = 0.2, visualize = FALSE, desired_number_of_groups = 5) {
#     split_type <- match.arg(split_type)
#     n <- nrow(X)
#     p <- ncol(X)
#     # Error messages
#     cat("Shape of X:", dim(X), "\n")
#     cat("Length of Y:", length(Y), "\n")
#     cat("First few rows of X:\n")
#     print(head(X[, 1:5]))  # Print first 5 columns
#     cat("First few values of Y:\n")
#     print(head(Y))

#     # Check for NAs
#     cat("Number of missing values in X:", sum(is.na(X)), "\n")
#     cat("Number of missing values in Y:", sum(is.na(Y)), "\n")
    
#     # Verify column variance
#     zero_var_cols <- which(apply(X, 2, var) == 0)
#     if (length(zero_var_cols) > 0) {
#     cat("Columns with zero variance:", zero_var_cols, "\n")
#     }
#     # Automatically determine groups if not provided
#     if (is.null(groups)) {
#         # throw an error asking for user to specify group option -- "none" if default
#         stop("Error: Groups not provided. Please specify a grouping option.")
#     } else if (groups == "corr") {
#         # Correlation-based grouping
#         corr_matrix <- cor(X)
#         dist_matrix <- as.dist(1 - abs(corr_matrix))
#         hc <- hclust(dist_matrix, method = "average")
#         if (is.null(desired_number_of_groups)) {
#              desired_number_of_groups <- p    
#         }
#         groups <- cutree(hc, k = desired_number_of_groups)
#         # debugging mesages
#         cat("Grouping length:", length(groups), "\n")
#         print(table(groups))  # Show distribution of groups
#     } else if (groups == "pca") {
#         # PCA-based grouping
#         # Perform PCA
#         pca_result <- prcomp(X, scale. = TRUE)
        
#         # Limit to the number of available features to avoid mismatch
#         num_pcs <- min(ncol(X), 10)  
#         loading_matrix <- abs(pca_result$rotation[, 1:num_pcs])
        
#         groups <- apply(loading_matrix, 1, which.max)
        
#         # Ensure correct length
#         if (length(groups) != ncol(X)) {
#             groups <- rep(1:ncol(X))  # Assign each feature to its own group as a fallback
#             warning("Mismatch between groups and features. Assigning default groups.")
#         }
#         # debugging messages
#         cat("Grouping length:", length(groups), "\n")
#         cat("distribution of grouops: ", table(groups), "\n")  # Show distribution of groups
#     } else {
#        groups <- rep(1:ncol(X), each = 1)  # Assign each feature to its own group, default option
#     }
#     index_final <- c(NA, groups) # account for intercept column 

#     # Ensure directories exist
#     if (visualize) {
#         plot_dir <- here::here("results", "eda", "grplasso")
#         if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
#     }
#     sink(DEBUG_LOG, append = TRUE)
#     cat("Checking groups alignment...\n")
#     cat("Number of groups provided:", length(groups), "\n")
#     cat("Number of features in X:", ncol(X), "\n")

#     if (length(groups) != ncol(X)) {
#         stop("Error: Mismatch between group length and number of features in X.")
#     }
#     sink()

#     # Train/validate/test split
#     set.seed(123)  # For reproducibility
#     train_idx <- sample(1:n, size = train_ratio * n)
#     if (split_type == "train_validate_test") {
#         valid_idx <- setdiff(sample(1:n, size = (train_ratio + valid_ratio) * n), train_idx)
#         test_idx <- setdiff(1:n, c(train_idx, valid_idx))
#     } else {
#         test_idx <- setdiff(1:n, train_idx)
#     }

#     X_train <- X[train_idx, ]
#     X_train <- cbind(Intercept = 1, X_train) # intercept
#     Y_train <- Y[train_idx]
#     X_test <- X[test_idx, ]
#     X_test <- cbind(Intercept = 1, X_test) # intercept
#     Y_test <- Y[test_idx]

#     # debug statements 
#     sink(DEBUG_LOG, append = TRUE)
#     cat("Train set size:", length(train_idx), "\n")
#     if (split_type == "train_validate_test") {
#         cat("Validation set size:", length(valid_idx), "\n")
#     }
#     cat("Test set size:", length(test_idx), "\n")

#     cat("Dimensions of training data:", dim(X_train), "\n")
#     cat("Dimensions of test data:", dim(X_test), "\n")
#     sink() 

#     # Fit grouped LASSO and select best lambda
#     lambdas <- exp(seq(log(1e-4), log(10), length.out = 100))


#     # Debugging print statements
#     sink(DEBUG_LOG, append = TRUE)
#     cat("Starting grouped LASSO analysis\n")
#     cat("Data dimensions: ", dim(X), "\n")
#     cat("Response vector length: ", length(Y), "\n")
#     cat("First few rows of X:\n")
#     print(head(X[, 1:min(5, ncol(X))]))  # Print first 5 columns safely
#     cat("Summary of Y:\n")
#     print(summary(Y))

#     cat("Number of groups provided:", length(groups), "\n")
#     cat("Unique groups:", unique(groups), "\n")
#     sink()

#     if (split_type == "train_validate_test") {
#         X_valid <- X[valid_idx, ]
#         X_valid <- cbind(Intercept = 1, X_valid) # intercept
#         Y_valid <- Y[valid_idx]
        
#         cv_results <- sapply(lambdas, function(lambda) {
#             sink(DEBUG_LOG, append = TRUE)
#             cat("Running cross-validation for lambda =", lambda, "\n")

#             fit <- grplasso(X_train, Y_train, index = index_final, lambda = lambda, model = LinReg())
#             mse <- mean((Y_valid - predict(fit, X_valid))^2)

#             cat("MSE for lambda =", lambda, ":", mse, "\n")
#             sink()  # Close sink after logging

#             return(mse)
#         })

#         best_lambda <- lambdas[which.min(cv_results)]
#     } else {
#         best_lambda <- lambdas[1]
#     }

#     # Final Model Fitting and Logging
#     final_model <- grplasso(X_train, Y_train, index = index_final, lambda = best_lambda, model = LinReg())
#     Y_test_pred <- predict(final_model, X_test)
#     test_mse <- mean((Y_test - Y_test_pred)^2)

#     sink(DEBUG_LOG, append = TRUE)
#     # Save coefficients
#     coefs <- as.numeric(coef(final_model))
#     if (any(is.na(coefs))) {
#         stop("NA values found in coefficients")
#     }
#     if (is.list(coefs)) {
#     cat("Warning: Coefficients are stored in a list. Converting to numeric...\n")
#     coefs <- unlist(coefs)
#     }

#     cat("First few coefficients:\n")
#     print(head(coefs))

#     # Ensure coefficients are numeric
#     if (!is.numeric(coefs)) {
#         stop("Error: Coefficients contain non-numeric values.")
#     }

#     coefs_no_intercept <- data.frame(
#         feats = colnames(X_train),
#         coef = coefs,
#         group = index_final,
#         is_intercept = is.na(index_final) # TRUE for intercept row
#     ) %>% arrange(desc(abs(coef)))

#     cat("Coefficient extraction successful. First few sorted coefficients:\n")
#     print(head(coefs_no_intercept))
#     sink()



#     sink(DEBUG_LOG, append = TRUE)

#     cat(str(coef(final_model)))  # Inspect the structure
#     cat(head(coef(final_model))) # Print some values
   
#     cat("Best lambda selected:", best_lambda, "\n")
#     cat("Final test MSE:", test_mse, "\n")
#     cat("Coefficients:\n")
#     print(head(coef(final_model)))
#     sink()

#     # Return results
#     list(
#         model = final_model,
#         test_mse = test_mse,
#         coefs = coefs_no_intercept, 
#         best_lambda = best_lambda
#     )

# }

#### TEST ####
# load data 
# vir_z_pat <- readRDS(here::here("data", "processed", "vir_z_pat.rds"))
# vir_lib <- readRDS(here::here("data", "rds", "vir_lib.rds"))
# vir_z <- readRDS(here::here("data", "rds", "vir_z.rds"))

#============================================================================
# parallel call to function 
#============================================================================
group_options <- c("corr", "pca") # ,"none")
results <- list()

# Detect number of available cores, leaving 1 core free for system stability
num_cores <- detectCores() - 1

# Define a function to process each group option in parallel
process_group_option <- function(group_option) {
    sink(DEBUG_LOG, append = TRUE)
    cat("Processing group option:", group_option, "\n")
    sink()
    
    result <- tryCatch({
        run_grplasso_reg_parallel(X_vir_z, y_vir_z, groups = group_option, split_type = "train_validate_test", visualize = TRUE)
    }, error = function(e) {
        sink(DEBUG_LOG, append = TRUE)
        cat("Error encountered for group option:", group_option, "\n")
        cat("Error message:", conditionMessage(e), "\n")
        sink()
        return(NULL)  # Return NULL for failed cases
    })

    if (!is.null(result)) {
        sink(DEBUG_LOG, append = TRUE)
        cat("Group option:", group_option, "\n")
        cat("Test MSE:", result$test_mse, "\n")
        cat("Best Lambda:", result$best_lambda, "\n\n")
        sink()
        return(list(group_option = group_option, result = result))
    } else {
        sink(DEBUG_LOG, append = TRUE)
        cat("Skipping group option due to an error:", group_option, "\n")
        sink()
        return(NULL)
    }
}

# Run in parallel across group options
parallel_results <- mclapply(group_options, process_group_option, mc.cores = num_cores)

# Convert list to named format for easier access
for (res in parallel_results) {
    if (!is.null(res)) {
        results[[res$group_option]] <- res$result
    }
}

# Save results as an RDS file
saveRDS(results, here::here("results", "eda", "grplasso", "grplasso_result_parallel.rds"))

# Load results and add organism names
results <- readRDS(here::here("results", "eda", "grplasso", "grplasso_result_parallel.rds"))
for (opt in names(results)) {
    organisms <- get_organism_by_id(results[[opt]]$coefs$feats[-1], vir_lib)
    organisms <- c("NA - intercept", organisms)
    results[[opt]]$coefs$organism <- organisms
}

saveRDS(results, here::here("results", "eda", "grplasso", "grplasso_result_parallel.rds"))

# Save results to a spreadsheet with separate sheets per grouping option
wb <- createWorkbook()
for (setting in names(results)) {
    addWorksheet(wb, setting)
    writeData(wb, sheet = setting, results[[setting]]$coefs)
}

results_dir <- here::here("results", "eda", "grplasso")
if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}


saveWorkbook(wb, file = file.path(RUN_DIR, EXCEL_FNAME), overwrite = TRUE)  


#============================================================================
# non parallel call to function, comment when testing parallel version
#============================================================================


# group_options <- c("corr", "pca") # , "none")
# results <- list()

# # Run the function with different group options and compare printed results
# for (group_option in group_options) {
#     sink(DEBUG_LOG, append = TRUE)

#     cat("Processing group option:", group_option, "\n")
#     result <- tryCatch({
#         run_grplasso_reg(X_vir_z, y_vir_z, groups = group_option, split_type = "train_validate_test", visualize = TRUE)
#     }, error = function(e) {
#         cat("Error encountered for group option:", group_option, "\n")
#         cat("Error message:", conditionMessage(e), "\n")
#         return(NULL)  # Return NULL for failed cases
#     })
    
#     if (!is.null(result)) {
#         results[[as.character(group_option)]] <- result
#         cat("Group option:", group_option, "\n")
#         cat("Test MSE:", result$test_mse, "\n")
#         cat("Best Lambda:", result$best_lambda, "\n\n")
#     } else {
#         cat("Skipping group option due to an error:", group_option, "\n")
#     }
#     sink()
# }

# # TODO logged 25-01-22-wed: still experiencing mismatch between number of groups andnumber of columns; try: remove intercept column when grouping 
# # TODO logged 25-01-23-thurs: forgot to include best_lambda when packaging up results into dataframe
# # save results as a RDS
# saveRDS(results, here::here("results", "eda", "grplasso", "grplasso_result.rds"))




# # Load results
# results <- readRDS(here::here("results", "eda", "grplasso", "grplasso_result.rds"))
# for (opt in names(results)) {
#     organisms <- get_organism_by_id(results[[opt]]$coefs$feats[-1], vir_lib)
#     organisms <- c("NA - intercept", organisms)
#     results[[opt]]$coefs$organism <- organisms
# }

# saveRDS(results, here::here("results", "eda", "grplasso", "grplasso_result.rds"))

# # save results to a spreadsheet, with separate sheet for each grouping option 
# wb <- createWorkbook()
# for (setting in names(results)) {
#     addWorksheet(wb, setting)
#     writeData(wb, sheet = setting, results[[setting]]$coef)
# }

# results_dir <- here::here("results", "eda", "grplasso")
# if (!dir.exists(results_dir)) {
#     dir.create(results_dir, recursive = TRUE)
# }

# grplasso_results_file <- "grplasso_result"
# dataset_name <- "virscan_all"
# final_file_name <- paste0(grplasso_results_file, "_", dataset_name, ".xlsx")
# saveWorkbook(wb, file = file.path(results_dir, final_file_name), overwrite = TRUE)
