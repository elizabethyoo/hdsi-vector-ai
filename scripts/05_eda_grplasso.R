library(dplyr)
library(ggplot2)
library(reshape2)
library(here)
library(gt)
library(stringr)
library(data.table)
library(openxlsx)
library(pbapply)    # For parallel processing with progress bar
library(grplasso)

# anchor project root directory 
here::i_am("scripts/05_eda_grplasso.R")

#============================================================================


# Randomly select peptides
test_peptide_columns <- sample(colnames(X_y_vir_z[,-1]), num_test_peptides)

# Create the subset dataframe
X_y_vir_z <- X_y_vir_z %>% select(all_of(c("COVID19_status", test_peptide_columns))) # uncomment for testing on small subset

X_vir_z <- X_y_vir_z %>%
    select(-COVID19_status)


y_vir_z <- X_y_vir_z$COVID19_status

# Confirm dimensions
cat("Peptide Data Dimensions:", dim(X_vir_z), "\n")

set.seed(123)

run_grplasso_reg <- function(X, Y, groups = NULL, split_type = c("train_validate_test", "train_test"), 
                             train_ratio = 0.6, valid_ratio = 0.2, visualize = FALSE, desired_number_of_groups = 5) {
    split_type <- match.arg(split_type)
    n <- nrow(X)
    p <- ncol(X)

    # Automatically determine groups if not provided
    if (is.null(groups)) {
        groups <- 1:p  # Default: Treat each feature as its own group
    } else if (groups == "corr") {
        # Correlation-based grouping
        corr_matrix <- cor(X)
        dist_matrix <- as.dist(1 - abs(corr_matrix))
        hc <- hclust(dist_matrix, method = "average")
        groups <- cutree(hc, k = desired_number_of_groups)
    } else if (groups == "pca") {
        # PCA-based grouping
        pca_result <- prcomp(X, scale. = TRUE)
        num_pcs <- min(10, ncol(pca_result$rotation))  # Limit PCs based on data dimensions
        loading_matrix <- abs(pca_result$rotation[, 1:num_pcs])
        groups <- apply(loading_matrix, 1, which.max)
    }

    # Ensure directories exist
    if (visualize) {
        plot_dir <- here::here("results", "eda", "grplasso")
        if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
    }

    # Train/validate/test split
    set.seed(123)  # For reproducibility
    train_idx <- sample(1:n, size = train_ratio * n)
    if (split_type == "train_validate_test") {
        valid_idx <- setdiff(sample(1:n, size = (train_ratio + valid_ratio) * n), train_idx)
        test_idx <- setdiff(1:n, c(train_idx, valid_idx))
    } else {
        test_idx <- setdiff(1:n, train_idx)
    }
    X_train <- X[train_idx, ]
    Y_train <- Y[train_idx]
    X_test <- X[test_idx, ]
    Y_test <- Y[test_idx]

    # Fit grouped LASSO and select best lambda
    lambdas <- exp(seq(log(1e-4), log(1), length.out = 100))
    if (split_type == "train_validate_test") {
        X_valid <- X[valid_idx, ]
        Y_valid <- Y[valid_idx]
        cv_results <- sapply(lambdas, function(lambda) {
            fit <- grplasso(X_train, Y_train, index = groups, lambda = lambda, model = LinReg())
            mse <- mean((Y_valid - predict(fit, X_valid))^2)
            return(mse)
        })
        best_lambda <- lambdas[which.min(cv_results)]
    } else {
        best_lambda <- lambdas[1]
    }

    # Fit final model and evaluate
    final_model <- grplasso(X_train, Y_train, index = groups, lambda = best_lambda, model = LinReg())
    Y_test_pred <- predict(final_model, X_test)
    test_mse <- mean((Y_test - Y_test_pred)^2)

    # Save coefficients
    coefficients <- data.frame(coef = coef(final_model), group = groups, index = 1:p) %>%
        arrange(desc(abs(coef)))

    # Return results
    list(model = final_model, test_mse = test_mse, coefficients = coefficients)
}

#### TEST ####
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
num_test_peptides <- 100

# Run the function with a subset of the data 
result <- run_grplasso_reg(X_vir_z, y_vir_z, groups = "corr", split_type = "train_validate_test", visualize = TRUE)

# save results as a RDS
saveRDS(result, here::here("results", "eda", "grplasso", "grplasso_result.rds"))

# Print the results
print(result$test_mse)
print(result$coefficients)
# simulate data to test the function
# Run the function with different group options and compare results
group_options <- c("corr", "pca", NULL)
results <- list()

for (group_option in group_options) {
    result <- run_grplasso_reg(X, Y, groups = group_option, split_type = "train_validate_test", visualize = TRUE)
    results[[as.character(group_option)]] <- result
    cat("Group option:", group_option, "\n")
    cat("Test MSE:", result$test_mse, "\n")
    cat("Best Lambda:", result$best_lambda, "\n\n")
}
