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
parallel,
purrr
)

# anchor project root directory 
here::i_am("scripts/05_eda_grplasso.R")

# helper function -- works on any of covscan, virscan, or combination of both so long as
# peptide id is in "c_" or "v_" format
# see formatted data in /n/home01/egraff/hdsi-vector-ai/data/processed/Shrock_vir_cov_combined
get_organism_by_id <- function(id_chr, cov_lib, vir_lib) {
  # Ensure the id columns are character for comparison if they are factors
  if (is.factor(cov_lib$id)) cov_lib$id <- as.character(cov_lib$id)
  if (is.factor(vir_lib$id)) vir_lib$id <- as.character(vir_lib$id)
  
  if (startsWith(id_chr, "c_")) {
    # If cov_lib$id is numeric, extract a numeric value; otherwise, use the full id string.
    if (is.numeric(cov_lib$id)) {
      id_val <- as.numeric(sub("^c_", "", id_chr))
    } else {
      id_val <- id_chr
    }
    row_data <- cov_lib %>% filter(id == id_val)
    if (nrow(row_data) == 0) {
      return(list(
        Organism = NA,
        Protein  = NA,
        Sequence = NA,
        Start    = NA,
        End      = NA
      ))
    }
    return(list(
      Organism = row_data %>% pull(Organism) %>% first(),
      Protein  = row_data %>% pull(Protein) %>% first(),
      Sequence = row_data %>% pull(Sequence) %>% first(),
      Start    = row_data %>% pull(start)    %>% first(),
      End      = row_data %>% pull(end)      %>% first()
    ))
    
  } else if (startsWith(id_chr, "v_")) {
    if (is.numeric(vir_lib$id)) {
      id_val <- as.numeric(sub("^v_", "", id_chr))
    } else {
      id_val <- id_chr
    }
    row_data <- vir_lib %>% filter(id == id_val)
    if (nrow(row_data) == 0) {
      return(list(
        Organism = NA,
        Protein  = NA,
        Sequence = NA,
        Start    = NA,
        End      = NA
      ))
    }
    return(list(
      Organism = row_data %>% pull(Organism) %>% first(),
      Protein  = row_data %>% pull(`Protein names`) %>% first(),
      Sequence = row_data %>% pull(Sequence) %>% first(),
      Start    = row_data %>% pull(start)    %>% first(),
      End      = row_data %>% pull(end)      %>% first()
    ))
  } else {
    stop("Unrecognized peptide id prefix (must be 'c_' or 'v_').")
  }
}

# Helper function to get the first N words of a string
get_first_n_words <- function(text, n = 10) {
  words <- unlist(strsplit(text, "\\s+"))  # Split the text into words
  paste(head(words, n), collapse = " ")  # Combine the first N words
}


# #Writing to file ===========================================================
#TEST ON A SUBSET
# Define number of peptides for testing - default empty string
NUM_TEST_PEPTIDES <- "100"

# CHECK BEFORE RUNNING - CHANGE TO APPROPRIATE NAMES ####################################################################
GRPLASSO_RESULTS_FNAME <- "grplasso-result"
DATASET_NAME <- paste0("combined", "_", NUM_TEST_PEPTIDES)
BASE_FNAME <- paste0(GRPLASSO_RESULTS_FNAME, "_", DATASET_NAME)
EXCEL_FNAME <- paste0(BASE_FNAME, "_", DATASET_NAME, ".xlsx")

# Create directory for this run
RUN_DIR <- here::here("results", "eda", "grplasso", paste0(BASE_FNAME, "_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))
if (!dir.exists(RUN_DIR)) {
    dir.create(RUN_DIR, recursive = TRUE)
}
# CHECK BEFORE RUNNING - CHANGE TO APPROPRIATE NAMES ####################################################################


# Writing to file ===========================================================

# Get number of cores from SLURM, fallback to detectCores() - 1 if not set
num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores() - 1))
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
# pull formatted dfs from /n/home01/egraff/hdsi-vector-ai/data/processed/Shrock_vir_cov_combined
DATA_DIR <- here::here("data", "processed", "Shrock_vir_cov_combined")
cat("Step: Loading and preparing data...\n")
com_z_pat <- readRDS(file.path(DATA_DIR, "combined_z_pat.rds")) # combined
vir_z_pat <- readRDS(file.path(DATA_DIR, "vir_z_pat_renamed_ids.rds")) # virscan
cov_z_pat <- readRDS(file.path(DATA_DIR, "cov_z_pat_renamed_ids.rds")) # covscan

cov_lib   <- readRDS(file.path(DATA_DIR, "cov_lib_renamed_ids.rds"))  # if needed
vir_lib   <- readRDS(file.path(DATA_DIR, "vir_lib_renamed_ids.rds"))  # if needed

# CHECK BEFORE RUNNING - CHANGE TO APPROPRIATE NAMES ####################################################################
# Prepare data  
# TODO isolate data processing step in a separate script s.t. nothing in analysis script is specific to dataset being used;
# abstract to X and y

META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "Hospitalized")

X_y_cov_z <- cov_z_pat %>%
  select(-all_of(META_COL_NO_STATUS)) %>%
  mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
  select(COVID19_status, everything())
  
# TEST ON A SUBSET ==========================================================
if (NUM_TEST_PEPTIDES != "") {
  set.seed(12345)  
  TEST_COL_FNAME <- paste0(BASE_FNAME, "_", "test_columns.rds")
  test_peptide_columns <- sample(colnames(X_y_cov_z[,-1]), as.numeric(NUM_TEST_PEPTIDES))
  # save in RUN_DIR
  saveRDS(test_peptide_columns, file.path(RUN_DIR, TEST_COL_FNAME))
  X_y_cov_z <- X_y_cov_z %>% select(all_of(c("COVID19_status", test_peptide_columns)))
}

X_cov_z <- X_y_cov_z %>% select(-COVID19_status)
y_cov_z <- X_y_cov_z$COVID19_status

cat("Peptide Data Dimensions:", dim(X_cov_z), "\n")
cat("Response Vector Length:", length(y_cov_z), "\n")

X <- X_cov_z %>%
  as.matrix()

y <- y_cov_z %>%
  as.numeric()
# CHECK BEFORE RUNNING - CHANGE TO APPROPRIATE NAMES ####################################################################




# #============================================================================
# # virscan data============================================================================
# ### data prep
# META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "IP_reagent", 
#                    "Hospitalized", "Pull_down_antibody", 
#                    "patient", "original_id_column") # excludes COVID19_status which will be a part of the final dataframe


# # X_y_vir_z includes z-scores and COVID19_status vs. X_vir_z includes z-scores only
# X_y_vir_z <- vir_z_pat %>%
#     select(-all_of(META_COL_NO_STATUS)) %>%
#     mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
#     # put COVID19_status as the first column
#     select(COVID19_status, everything())

# # TEST ON A SUBSET ==========================================================
# set.seed(123)
# if (NUM_TEST_PEPTIDES != "") {
#   TEST_COL_FNAME <- paste0(BASE_FNAME, "_", "test_columns.rds")
#   test_peptide_columns <- sample(colnames(X_y_vir_z[,-1]), as.numeric(NUM_TEST_PEPTIDES))
#   saveRDS(test_peptide_columns, file.path(RUN_DIR, TEST_COL_FNAME))
#   X_y_vir_z <- X_y_vir_z %>% select(all_of(c("COVID19_status", test_peptide_columns)))
# }

# X_vir_z <- X_y_vir_z %>% select(-COVID19_status)
# y_vir_z <- X_y_vir_z$COVID19_status

# X <- X_vir_z %>%
#   as.matrix()

# y <- y_vir_z %>%
#   as.numeric()

# # Confirm dimensions
# sink(DEBUG_LOG, append = TRUE)
# cat("Peptide Data Dimensions:", dim(X), "\n")
# cat("Response Vector Length:", length(y), "\n")

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

    if (groups == "Organism") {
    # Extract peptide IDs from the column names of X
    peptide_ids <- colnames(X)
    if (is.null(peptide_ids)) {
        stop("Error: X must have column names representing peptide IDs.")
    }
    
    # For each peptide ID, apply get_organism_by_id() to retrieve the organism value
    organisms <- sapply(peptide_ids, function(pid) {
        result <- get_organism_by_id(pid, cov_lib, vir_lib)
        return(result$Organism)
    })
    
    # Identify indices for valid (non-NA) organism entries
    valid_idx <- which(!is.na(organisms))
    if (length(valid_idx) < length(organisms)) {
        warning("Excluding ", length(organisms) - length(valid_idx), " peptides with NA organism.")
    }
    
    # Subset X (and peptide_ids and organisms) to exclude peptides with NA organism
    X <<- X[, valid_idx, drop = FALSE]  # Using <<- to update X in the parent environment if needed
    peptide_ids <- peptide_ids[valid_idx]
    organisms <- organisms[valid_idx]
    
    # Convert organism labels to a numeric group index and build index_final
    group_indices <- as.integer(as.factor(organisms))
    index_final <- c(NA, group_indices)
    }

    set.seed(12345)
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
        clusterExport(cl, varlist = c("X_train", "Y_train", "index_final", "X_valid", "Y_valid", "grplasso", "DEBUG_LOG"), envir = environment())
        clusterEvalQ(cl, library(grplasso))
        clusterEvalQ(cl, library(dplyr))

        cv_results <- parSapply(cl, lambdas, function(lambda) {
            sink(DEBUG_LOG, append = TRUE)
            cat("Running cross-validation for lambda =", lambda, "\n")
            sink()

            fit <- grplasso(X_train, Y_train, index = index_final, lambda = lambda, model = LogReg())
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

    final_model <- grplasso(X_train, Y_train, index = index_final, lambda = best_lambda, model = LogReg())
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

#============================================================================
# parallel call to function 
#============================================================================

# Detect number of available cores, leaving 1 core free for system stability
num_cores <- detectCores() - 1
PAR_RESULTS_FNAME <- paste0(BASE_FNAME, "_", "fin.rds")
one_grp <- c("Organism")
# call grplasso function for single group option "Organism"
results <- run_grplasso_reg_parallel(X, y, groups = one_grp, split_type = "train_validate_test", visualize = TRUE, num_cores = num_cores)
saveRDS(results, file.path(RUN_DIR, PAR_RESULTS_FNAME))
cat("save results as RDS in", RUN_DIR, "\n")

# # Define a function to process each group option in parallel
# process_group_option <- function(group_option) {
#   sink(DEBUG_LOG, append = TRUE)
#   cat("Processing group option:", group_option, "\n")
#   sink()
  
#   result <- tryCatch({
#     run_grplasso_reg_parallel(X, y, groups = group_option, split_type = "train_validate_test", visualize = TRUE, num_cores = num_cores)
#   }, error = function(e) {
#     sink(DEBUG_LOG, append = TRUE)
#     cat("Error encountered for group option:", group_option, "\n")
#     cat("Error message:", conditionMessage(e), "\n")
#     sink()
#     return(NULL)  # Return NULL for failed cases
#   })

#   if (!is.null(result)) {
#     sink(DEBUG_LOG, append = TRUE)
#     cat("Group option:", group_option, "\n")
#     cat("Test MSE:", result$test_mse, "\n")
#     cat("Best Lambda:", result$best_lambda, "\n\n")
#     sink()
#     return(list(group_option = group_option, result = result))
#   } else {
#     sink(DEBUG_LOG, append = TRUE)
#     cat("Skipping group option due to an error:", group_option, "\n")
#     sink()
#     return(NULL)
#   }
# }


# results <- list()
# # Convert list to named format for easier access
# for (res in parallel_results) {
#     if (!is.null(res)) {
#         results[[res$group_option]] <- res$result
#     }
# }


#============================================================================
# post-processing
#============================================================================
RUN_DIR <- "/n/home01/egraff/hdsi-vector-ai/results/eda/grplasso/grplasso-result_covscan_100_2025-02-25_13-16-48"
PAR_RESULTS_FNAME <- "grplasso-result_covscan_100_fin.rds"
# Load results and add organism names
results <- readRDS(file.path(RUN_DIR, PAR_RESULTS_FNAME))

# Summarize key model results
cat("Best lambda selected:", results$best_lambda, "\n")
cat("Test MSE:", results$test_mse, "\n")

# Create a summary table of coefficients excluding the intercept
coefs_df <- results$coefs
coefs_no_int <- subset(coefs_df, is_intercept == FALSE)
coefs_no_int <- coefs_no_int %>%
  arrange(desc(coef)) %>% # Sort by coefficient value
  mutate(
    attributes = map(feats, get_organism_by_id, cov_lib, vir_lib),
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
  select(-first_n_words, -first_char_protein, -start_end)

# # sanity check -- verified there are ten unique (group, organism) pairs
# unique_pairs <- coefs_no_int %>%
#   distinct(group, organism)

# Save the coefficients to a CSV file
COEF_FNAME_CSV <- paste0(BASE_FNAME, "_", "coefs.csv")
write.csv(coefs_no_int, here::here(RUN_DIR, COEF_FNAME_CSV), row.names = FALSE)

# Visualization =================================================================
# Define a common color scale using the default hue scale
common_colors <- scale_fill_hue()

# Plot 1: Horizontal bar plot (ensure same color scale)

p1 <- ggplot(coefs_no_int, aes(x = reorder(feats, coef), y = coef, fill = factor(group))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Peptide ID", y = "Coefficient",
       fill = paste0("Group (Organism) ", organism),
       title = "Group Lasso Coefficients (Excluding Intercept)") +
  common_colors +
  theme_minimal()

print(p1)

print(p1)

# Plot 2: Boxplot using the exact same color scale as p1
p2 <- ggplot(coefs_no_int, aes(x = factor(group), y = coef, fill = factor(group))) +
  geom_boxplot() +
  labs(x = "Group", y = "Coefficient",
       title = "Coefficient Distribution by Group") +
  common_colors +
  theme_minimal()

print(p2)

for (opt in names(results)) {
  organisms <- sapply(results[[opt]]$coefs$feats[-1], function(x) {
      get_organism_by_id(x, cov_lib, vir_lib)$Organism
  })
  organisms <- c("NA - intercept", organisms)
  results[[opt]]$coefs$organism <- organisms
}

# # save results as RDS in RUN_DIR
# saveRDS(results, file.path(RUN_DIR, PAR_RESULTS_FNAME))

# # Save results to a spreadsheet with separate sheets per grouping option
# wb <- createWorkbook()
# for (setting in names(results)) {
#     addWorksheet(wb, setting)
#     writeData(wb, sheet = setting, results[[setting]]$coefs)
# }

# saveWorkbook(wb, file = file.path(RUN_DIR, EXCEL_FNAME), overwrite = TRUE)  

# # end of script
# cat("Script complete! Check outputs in:", RUN_DIR, "\n")







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

# Detect number of available cores, leaving 1 core free for system stability
# num_cores <- detectCores() - 1

# # Define a function to process each group option in parallel
# process_group_option <- function(group_option) {
#     sink(DEBUG_LOG, append = TRUE)
#     cat("Processing group option:", group_option, "\n")
#     sink()
    
#     result <- tryCatch({
#         run_grplasso_reg_parallel(X_vir_z, y_vir_z, groups = group_option, split_type = "train_validate_test", visualize = TRUE)
#     }, error = function(e) {
#         sink(DEBUG_LOG, append = TRUE)
#         cat("Error encountered for group option:", group_option, "\n")
#         cat("Error message:", conditionMessage(e), "\n")
#         sink()
#         return(NULL)  # Return NULL for failed cases
#     })

#     if (!is.null(result)) {
#         sink(DEBUG_LOG, append = TRUE)
#         cat("Group option:", group_option, "\n")
#         cat("Test MSE:", result$test_mse, "\n")
#         cat("Best Lambda:", result$best_lambda, "\n\n")
#         sink()
#         return(list(group_option = group_option, result = result))
#     } else {
#         sink(DEBUG_LOG, append = TRUE)
#         cat("Skipping group option due to an error:", group_option, "\n")
#         sink()
#         return(NULL)
#     }
# }

# # Run in parallel across group options
# parallel_results <- mclapply(group_options, process_group_option, mc.cores = num_cores)

# # Convert list to named format for easier access
# for (res in parallel_results) {
#     if (!is.null(res)) {
#         results[[res$group_option]] <- res$result
#     }
# }

# # Load results and add organism names
# results <- readRDS(file.path(RUN_DIR, PAR_RESULTS_FNAME))


# for (opt in names(results)) {
#     organisms <- get_organism_by_id(results[[opt]]$coefs$feats[-1], cov_lib, vir_lib)
#     organisms <- c("NA - intercept", organisms)
#     results[[opt]]$coefs$organism <- organisms
# }

# # save results as RDS in RUN_DIR
# saveRDS(results, file.path(RUN_DIR, PAR_RESULTS_FNAME))

# # Save results to a spreadsheet with separate sheets per grouping option
# wb <- createWorkbook()
# for (setting in names(results)) {
#     addWorksheet(wb, setting)
#     writeData(wb, sheet = setting, results[[setting]]$coefs)
# }

# saveWorkbook(wb, file = file.path(RUN_DIR, EXCEL_FNAME), overwrite = TRUE)  


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
