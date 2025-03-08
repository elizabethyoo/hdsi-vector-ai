# Given a m by n matrix dataframe or tibble, where m is the number of rows and n is the number of columns, 
# and a vector of length m, this script runs a LASSO regression and returns the coefficients of the model.

# load packages
library(glmnet)
library(here)
library(dplyr)

# anchor project root directory 
here::i_am("scripts/01_run_regression.R")


# Function to perform LASSO regression with optional train/validate/test split and visualization
run_lasso_regression <- function(X, y, split_type = c("train_validate_test", "train_test"), train_ratio = 0.6, valid_ratio = 0.2, visualize = FALSE) {
    split_type <- match.arg(split_type)
    
    n <- nrow(X)
    row_names <- rownames(X)
    
    # Initialize result list
    result <- list()
    
    if (split_type == "train_validate_test") {
        train_idx <- sample(1:n, size = train_ratio * n)
        valid_idx <- setdiff(sample(1:n, size = (train_ratio + valid_ratio) * n), train_idx)
        test_idx <- setdiff(1:n, c(train_idx, valid_idx))
        
        X_train <- X[train_idx, ]
        y_train <- y[train_idx]
        X_valid <- X[valid_idx, ]
        y_valid <- y[valid_idx]
        X_test <- X[test_idx, ]
        y_test <- y[test_idx]
        
        # Cross-validated LASSO
        cv_fit <- cv.glmnet(X_train, y_train, alpha = 1)
        best_lambda <- cv_fit$lambda.min
        
        # Predict on validation and test sets
        valid_pred <- as.numeric(predict(cv_fit, s = best_lambda, newx = X_valid))
        test_pred <- as.numeric(predict(cv_fit, s = best_lambda, newx = X_test))
        
        # Evaluate performance
        mse_valid <- mean((y_valid - valid_pred)^2)
        mse_test <- mean((y_test - test_pred)^2)
        
        # Coefficients at the best lambda
        lasso_coefficients <- coef(cv_fit, s = best_lambda)
        
        # Save relevant information to result
        result$model <- cv_fit
        result$best_lambda <- best_lambda
        result$mse_valid <- mse_valid
        result$mse_test <- mse_test
        result$coefficients <- as.matrix(lasso_coefficients)  # Convert to matrix for easier handling
        result$lambda_sequence <- cv_fit$lambda
        result$cvm <- cv_fit$cvm  # Mean cross-validation errors
        result$train_idx <- train_idx
        result$valid_idx <- valid_idx
        result$test_idx <- test_idx
        
        if (visualize) {
            # Visualization of the process of optimizing lambda
            plot_path_lambda <- here::here("results", paste0("01_lasso_fit_lambda_", deparse(substitute(y)), "_", deparse(substitute(X)), ".png"))
            png(filename = plot_path_lambda)
            plot(cv_fit)
            abline(v = log(cv_fit$lambda.min), col = "blue", lty = 2) # Mark the optimal lambda
            title("Cross-Validation Curve for LASSO")
            dev.off()
            
            # Visualization of how the MSEs evolve with the process
            plot_path_mse <- here::here("results", paste0("01_lasso_fit_mse_lambda_", deparse(substitute(y)), "_", deparse(substitute(X)), ".png"))
            png(filename = plot_path_mse)
            plot(log(cv_fit$lambda), cv_fit$cvm, type = "b", pch = 19, col = "red", xlab = "Log(Lambda)", ylab = "Mean Squared Error", main = "MSE vs Log(Lambda)")
            points(log(cv_fit$lambda.min), min(cv_fit$cvm), col = "blue", pch = 19)
            abline(v = log(cv_fit$lambda.min), col = "blue", lty = 2)
            legend("topright", legend = c("MSE", "Optimal Lambda"), col = c("red", "blue"), pch = 19, lty = 1:2)
            dev.off()
            
            result$plot_path_lambda <- plot_path_lambda
            result$plot_path_mse <- plot_path_mse
        }
        
    } else if (split_type == "train_test") {
        train_idx <- sample(1:n, size = train_ratio * n)
        test_idx <- setdiff(1:n, train_idx)
        
        X_train <- X[train_idx, ]
        y_train <- y[train_idx]
        X_test <- X[test_idx, ]
        y_test <- y[test_idx]
        
        # Cross-validated LASSO
        cv_fit <- cv.glmnet(X_train, y_train, alpha = 1)
        best_lambda <- cv_fit$lambda.min
        
        # Predict on test set
        test_pred <- as.numeric(predict(cv_fit, s = best_lambda, newx = X_test))
        
        # Evaluate performance
        mse_test <- mean((y_test - test_pred)^2)
        
        # Coefficients at the best lambda
        lasso_coefficients <- coef(cv_fit, s = best_lambda)
        # Dataframe with coefficients in descending order of magnitude, their corresponding group, and index in the original data matrix
        lasso_coefficients <- data.frame(coef = as.numeric(lasso_coefficients), index = 1:ncol(X)) %>%
            arrange(desc(abs(coef)))
            
        # Save relevant information to result
        result$model <- cv_fit
        result$best_lambda <- best_lambda
        result$mse_test <- mse_test
        result$coefficients <- as.matrix(lasso_coefficients)
        result$lambda_sequence <- cv_fit$lambda
        result$cvm <- cv_fit$cvm
        result$train_idx <- train_idx
        result$test_idx <- test_idx
    }
    
    return(result)
}

vir_z_pat <- readRDS(here::here("data", "processed", "vir_z_pat.rds"))
vir_z_pat <- as_tibble(vir_z_pat)
PAT_META_COLS <- c("rep_id", "Sample", "Library", "IP_reagent", "COVID19_status", "Hospitalized", "Pull_down_antibody", "patient", "original_id_column")
X_vir_z <- vir_z_pat %>% select(-all_of(PAT_META_COLS))
# convert y to numeric 1 if "positive" and 0 if "negative"
y_vir_z <- vir_z_pat %>% mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>% pull(COVID19_status)

# Assuming X_syn is a tibble with column names and row names
# and y_syn is a vector with row names matching X_syn
# Run the LASSO regression
result <- run_lasso_regression(X = as.matrix(X_vir_z), y = y_vir_z, visualize = TRUE)
# Save lasso results for future inspection 
saveRDS(result, here::here("results", "lasso_full_result.rds"))
result <- readRDS(here::here("results", "lasso_full_result.rds"))
print(result)

# # Map the coefficients back to their corresponding column names
# coefficients <- result$coefficients
# coefficients_df <- as.data.frame(as.matrix(coefficients))
# coefficients_df$feature <- rownames(coefficients_df)
# coefficients_df <- coefficients_df[coefficients_df$feature != "(Intercept)", ]

# # Print the coefficients with their corresponding feature names
# print(coefficients_df)