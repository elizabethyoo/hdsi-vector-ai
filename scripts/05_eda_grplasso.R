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
DEBUG_LOG <- "grplasso_debug_log.txt"

# load data
vir_z_pat <- readRDS(here::here("data", "processed", "vir_z_pat.rds"))
vir_lib <- readRDS(here::here("data", "rds", "vir_lib.rds"))
vir_z <- readRDS(here::here("data", "rds", "vir_z.rds"))
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
    select(COVID19_status, everything())

# TEST ON A SUBSET
# Define number of peptides for testing
num_test_peptides <- 60

# Randomly select peptides
test_peptide_columns <- sample(colnames(X_y_vir_z[,-1]), num_test_peptides)

# Create the subset dataframe
X_y_vir_z <- X_y_vir_z %>% select(all_of(c("COVID19_status", test_peptide_columns))) # uncomment for testing on small subset

X_vir_z <- X_y_vir_z %>%
    select(-COVID19_status) %>%
    as.matrix()


y_vir_z <- X_y_vir_z$COVID19_status %>%
    as.numeric()

# Confirm dimensions
sink(DEBUG_LOG, append = TRUE)
cat("Peptide Data Dimensions:", dim(X_vir_z), "\n")

### random forest implementation 

set.seed(123)

run_grplasso_reg <- function(X, Y, groups = NULL, split_type = c("train_validate_test", "train_test"), 
                             train_ratio = 0.6, valid_ratio = 0.2, visualize = FALSE, desired_number_of_groups = 5) {
    split_type <- match.arg(split_type)
    n <- nrow(X)
    p <- ncol(X)
    # Error messages
    cat("Shape of X:", dim(X), "\n")
    cat("Length of Y:", length(Y), "\n")
    cat("First few rows of X:\n")
    print(head(X[, 1:5]))  # Print first 5 columns
    cat("First few values of Y:\n")
    print(head(Y))

    # Check for NAs
    cat("Number of missing values in X:", sum(is.na(X)), "\n")
    cat("Number of missing values in Y:", sum(is.na(Y)), "\n")
    
    # Verify column variance
    zero_var_cols <- which(apply(X, 2, var) == 0)
    if (length(zero_var_cols) > 0) {
    cat("Columns with zero variance:", zero_var_cols, "\n")
    }
    # Automatically determine groups if not provided
    if (is.null(groups)) {
        groups <- 1:p  # Default: Treat each feature as its own group
        # debugging messages
        # Check grouping
        cat("Grouping length:", length(groups), "\n")
        print(table(groups))  # Show distribution of groups
    } else if (groups == "corr") {
        # Correlation-based grouping
        corr_matrix <- cor(X)
        dist_matrix <- as.dist(1 - abs(corr_matrix))
        hc <- hclust(dist_matrix, method = "average")
        if (is.null(desired_number_of_groups)) {
             desired_number_of_groups <- p    
        }
        groups <- cutree(hc, k = desired_number_of_groups)
        # debugging mesages
        cat("Grouping length:", length(groups), "\n")
        print(table(groups))  # Show distribution of groups
    } else if (groups == "pca") {
        # PCA-based grouping
        # Perform PCA
        pca_result <- prcomp(X, scale. = TRUE)
        
        # Limit to the number of available features to avoid mismatch
        num_pcs <- min(ncol(X), 10)  
        loading_matrix <- abs(pca_result$rotation[, 1:num_pcs])
        
        groups <- apply(loading_matrix, 1, which.max)
        
        # Ensure correct length
        if (length(groups) != ncol(X)) {
            groups <- rep(1:ncol(X))  # Assign each feature to its own group as a fallback
            warning("Mismatch between groups and features. Assigning default groups.")
        }
        # debugging messages
        cat("Grouping length:", length(groups), "\n")
        cat("distribution of grouops: ", table(groups), "\n")  # Show distribution of groups
    }

    # Ensure directories exist
    if (visualize) {
        plot_dir <- here::here("results", "eda", "grplasso")
        if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
    }
    sink(DEBUG_LOG, append = TRUE)
    cat("Checking groups alignment...\n")
    cat("Number of groups provided:", length(groups), "\n")
    cat("Number of features in X:", ncol(X), "\n")

    if (length(groups) != ncol(X)) {
        stop("Error: Mismatch between group length and number of features in X.")
    }
    sink()

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
    X_train <- cbind(Intercept = 1, X_train) # intercept
    Y_train <- Y[train_idx]
    X_test <- X[test_idx, ]
    X_test <- cbind(Intercept = 1, X_test) # intercept
    Y_test <- Y[test_idx]

    # debug statements 
    sink(DEBUG_LOG, append = TRUE)
    cat("Train set size:", length(train_idx), "\n")
    if (split_type == "train_validate_test") {
        cat("Validation set size:", length(valid_idx), "\n")
    }
    cat("Test set size:", length(test_idx), "\n")

    cat("Dimensions of training data:", dim(X_train), "\n")
    cat("Dimensions of test data:", dim(X_test), "\n")
    sink() 

    # Fit grouped LASSO and select best lambda
    lambdas <- exp(seq(log(1e-4), log(10), length.out = 100))


    # Debugging print statements
    sink(DEBUG_LOG, append = TRUE)
    cat("Starting grouped LASSO analysis\n")
    cat("Data dimensions: ", dim(X), "\n")
    cat("Response vector length: ", length(Y), "\n")
    cat("First few rows of X:\n")
    print(head(X[, 1:min(5, ncol(X))]))  # Print first 5 columns safely
    cat("Summary of Y:\n")
    print(summary(Y))

    cat("Number of groups provided:", length(groups), "\n")
    cat("Unique groups:", unique(groups), "\n")
    sink()

    if (split_type == "train_validate_test") {
        X_valid <- X[valid_idx, ]
        X_valid <- cbind(Intercept = 1, X_valid) # intercept
        Y_valid <- Y[valid_idx]
        
        cv_results <- sapply(lambdas, function(lambda) {
            sink(DEBUG_LOG, append = TRUE)
            cat("Running cross-validation for lambda =", lambda, "\n")

            fit <- grplasso(X_train, Y_train, index = groups, lambda = lambda, model = LinReg())
            mse <- mean((Y_valid - predict(fit, X_valid))^2)

            cat("MSE for lambda =", lambda, ":", mse, "\n")
            sink()  # Close sink after logging

            return(mse)
        })

        best_lambda <- lambdas[which.min(cv_results)]
    } else {
        best_lambda <- lambdas[1]
    }

    # Final Model Fitting and Logging
    final_model <- grplasso(X_train, Y_train, index = groups, lambda = best_lambda, model = LinReg())
    Y_test_pred <- predict(final_model, X_test)
    test_mse <- mean((Y_test - Y_test_pred)^2)

    sink(DEBUG_LOG, append = TRUE)
    # Save coefficients
    coefs <- as.numeric(coef(final_model))
    if (any(is.na(coefs))) {
        stop("NA values found in coefficients")
    }
    if (is.list(coefs)) {
    cat("Warning: Coefficients are stored in a list. Converting to numeric...\n")
    coefs <- unlist(coefs)
    }


    cat("First few coefficients:\n")
    print(head(coefs))

    # Ensure coefficients are numeric
    if (!is.numeric(coefs)) {
        stop("Error: Coefficients contain non-numeric values.")
    }


    coefficients <- data.frame(coef = coefs, group = groups, index = 1:p) %>%
        arrange(desc(abs(coef)))


    cat("Coefficient extraction successful. First few sorted coefficients:\n")
    print(head(coefficients))
    sink()



    sink(DEBUG_LOG, append = TRUE)

    cat(str(coef(final_model)))  # Inspect the structure
    cat(head(coef(final_model))) # Print some values
   
    cat("Best lambda selected:", best_lambda, "\n")
    cat("Final test MSE:", test_mse, "\n")
    cat("Coefficients:\n")
    print(head(coef(final_model)))
    sink()

    # Return results
    list(model = final_model, test_mse = test_mse, coefficients = coefficients)

}

#### TEST ####
# load data 
# vir_z_pat <- readRDS(here::here("data", "processed", "vir_z_pat.rds"))
# vir_lib <- readRDS(here::here("data", "rds", "vir_lib.rds"))
# vir_z <- readRDS(here::here("data", "rds", "vir_z.rds"))

#============================================================================


group_options <- c("corr", "pca", NULL)
results <- list()

# Run the function with different group options and compare printed results
for (group_option in group_options) {
    sink(DEBUG_LOG, append = TRUE)

    cat("Processing group option:", group_option, "\n")
    result <- tryCatch({
        run_grplasso_reg(X_vir_z, y_vir_z, groups = group_option, split_type = "train_validate_test", visualize = TRUE)
    }, error = function(e) {
        cat("Error encountered for group option:", group_option, "\n")
        cat("Error message:", conditionMessage(e), "\n")
        return(NULL)  # Return NULL for failed cases
    })
    
    if (!is.null(result)) {
        results[[as.character(group_option)]] <- result
        cat("Group option:", group_option, "\n")
        cat("Test MSE:", result$test_mse, "\n")
        cat("Best Lambda:", result$best_lambda, "\n\n")
    } else {
        cat("Skipping group option due to an error:", group_option, "\n")
    }
    sink()
}

# TODO logged 25-01-22-wed: still experiencing mismatch between number of groups andnumber of columns; try: remove intercept column when grouping 

# save results as a RDS
saveRDS(results, here::here("results", "eda", "grplasso", "grplasso_result.rds"))
