library(grplasso)
library(here)
library(dplyr)

run_grplasso_reg <- function(X, Y, groups = NULL, split_type = c("train_validate_test", "train_test"), train_ratio = 0.6, valid_ratio = 0.2, visualize = FALSE) {
    split_type <- match.arg(split_type)
    
    n <- nrow(X)
    p <- ncol(X)

    # Automatically determine groups if not provided
    if (is.null(groups)) {
        groups <- 1:p  # Default: Treat each feature as its own group
        # Validate groups (only applicable for default grouping)
        if (length(groups) != p) {
            stop("Length of 'groups' must match the number of columns in X when groups are treated individually.")
        }
    } else if (groups == "corr") {
        # Correlation-based grouping
        corr_matrix <- cor(X)
        dist_matrix <- as.dist(1 - abs(corr_matrix))  # Convert to a distance matrix
        hc <- hclust(dist_matrix, method = "average")  # Hierarchical clustering
        desired_number_of_groups <- 5  # Define the desired number of groups
        groups <- cutree(hc, k = desired_number_of_groups)  # Define group assignments
    } else if (groups == "pca") {
        # PCA-based grouping
        pca_result <- prcomp(X, scale. = TRUE)
        num_pcs <- 10  # Number of principal components to use for grouping
        loading_matrix <- abs(pca_result$rotation[, 1:num_pcs])
        groups <- apply(loading_matrix, 1, which.max)  # Assign each feature to its strongest PC
    }

    
    # Initialize result list
    result <- list()
    
    # Train/validate/test split
    if (split_type == "train_validate_test") {
        train_idx <- sample(1:n, size = train_ratio * n)
        valid_idx <- setdiff(sample(1:n, size = (train_ratio + valid_ratio) * n), train_idx)
        test_idx <- setdiff(1:n, c(train_idx, valid_idx))
        
        X_train <- X[train_idx, ]
        Y_train <- Y[train_idx]
        X_valid <- X[valid_idx, ]
        Y_valid <- Y[valid_idx]
        X_test <- X[test_idx, ]
        Y_test <- Y[test_idx]
    } else if (split_type == "train_test") {
        train_idx <- sample(1:n, size = train_ratio * n)
        test_idx <- setdiff(1:n, train_idx)
        
        X_train <- X[train_idx, ]
        Y_train <- Y[train_idx]
        X_test <- X[test_idx, ]
        Y_test <- Y[test_idx]
    }
    
    # Fit grouped LASSO with cross-validation
    lambdas <- exp(seq(log(1e-4), log(1), length.out = 100))  # Define a range of lambdas
    if (split_type == "train_validate_test") {
        cv_results <- sapply(lambdas, function(lambda) {
            fit <- grplasso(X_train, Y_train, index = groups, lambda = lambda, model = LinReg())
            mse <- mean((Y_valid - predict(fit, X_valid))^2)
            return(mse)
        })
    } else {
        # For train_test, no validation set is available, so no cross-validation
        cv_results <- NULL
    }
    
    # Select best lambda and fit final model
    if (!is.null(cv_results)) {
        best_lambda <- lambdas[which.min(cv_results)]
    } else {
        best_lambda <- lambdas[1]  # Default to first lambda if no CV
    }
    
    final_model <- grplasso(X_train, Y_train, index = groups, lambda = best_lambda, model = LinReg())
    
    # Predict and evaluate performance on test set
    Y_test_pred <- predict(final_model, X_test)
    test_mse <- mean((Y_test - Y_test_pred)^2)
    
    # Save coefficients
    coefficients <- as.numeric(coef(final_model))
    # Dataframe with coefficients in descending order of magnitude, their corresponding group, and index in the original data matrix
    coefficients <- data.frame(coef = coefficients, group = groups, index = 1:p) %>%
        arrange(desc(abs(coef)))


    # Save results
    result$model <- final_model
    result$best_lambda <- best_lambda
    result$test_mse <- test_mse
    result$coefficients <- coefficients
    result$cv_results <- if (!is.null(cv_results)) data.frame(lambda = lambdas, mse = cv_results) else NULL
    result$train_idx <- train_idx
    result$valid_idx <- if (split_type == "train_validate_test") valid_idx else NULL
    result$test_idx <- test_idx
    
    if (visualize && !is.null(cv_results)) {
        # Cross-validation MSE vs lambda plot
        plot_path <- here::here("results", "02_run_grplasso_cv_results.png")
        png(plot_path)
        plot(log(lambdas), cv_results, type = "b", col = "blue", xlab = "Log(Lambda)", ylab = "Validation MSE",
             main = "Cross-Validation MSE vs Lambda")
        abline(v = log(best_lambda), col = "red", lty = 2)  # Optimal lambda
        dev.off()
        result$plot_path <- plot_path
    }
    
    return(result)
}

#### TEST ####
# Simulate data
set.seed(123)
n <- 100  # Number of samples
p <- 50   # Number of features
X <- matrix(rnorm(n * p), nrow = n, ncol = p)

# Simulate sparse data
beta <- rep(0, p)
sparse_indices <- sample(1:p, size = 10)  # Only 10 non-zero coefficients
beta[sparse_indices] <- rnorm(10)

# Generate response variable with some noise
Y <- X %*% beta + rnorm(n)

# Run the function with the simulated data
result <- run_grplasso_reg(X, Y, groups = "corr", split_type = "train_validate_test", visualize = TRUE)


# # Print the results
# print(result$test_mse)
# print(result$coefficients)
# # simulate data to test the function
# # Run the function with different group options and compare results
# group_options <- c("corr", "pca", NULL)
# results <- list()

# for (group_option in group_options) {
#     result <- run_grplasso_reg(X, Y, groups = group_option, split_type = "train_validate_test", visualize = TRUE)
#     results[[as.character(group_option)]] <- result
#     cat("Group option:", group_option, "\n")
#     cat("Test MSE:", result$test_mse, "\n")
#     cat("Best Lambda:", result$best_lambda, "\n\n")
# }
