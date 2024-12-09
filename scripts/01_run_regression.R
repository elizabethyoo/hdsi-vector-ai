# Given a m by n matrix dataframe or tibble, where m is the number of rows and n is the number of columns, 
# and a vector of length m, this script runs a LASSO regression and returns the coefficients of the model.

# load packages
library(glmnet)
library(here)

# anchor project root directory 
here::i_am("scripts/01_run_regression.R")

# Function to perform LASSO regression with optional train/validate/test split and visualization
run_lasso_regression <- function(X, y, split_type = c("train_validate_test", "train_test"), train_ratio = 0.6, valid_ratio = 0.2, visualize = FALSE) {
    split_type <- match.arg(split_type)
    
    n <- nrow(X)
    
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
        
        if (visualize) {
            # Visualization of the process of optimizing lambda
            png(filename = here::here("results", paste0("01_lasso_fit_lambda_", deparse(substitute(y)), "_", deparse(substitute(X)), ".png")))
            plot(cv_fit)
            abline(v = log(cv_fit$lambda.min), col = "blue", lty = 2) # Mark the optimal lambda
            title("Cross-Validation Curve for LASSO")
            dev.off()
            
            # Visualization of how the MSEs evolve with the process
            png(filename = here::here("results", paste0("01_lasso_fit_mse_lambda_", deparse(substitute(y)), "_", deparse(substitute(X)), ".png")))
            plot(log(cv_fit$lambda), cv_fit$cvm, type = "b", pch = 19, col = "red", xlab = "Log(Lambda)", ylab = "Mean Squared Error", main = "MSE vs Log(Lambda)")
            points(log(cv_fit$lambda.min), min(cv_fit$cvm), col = "blue", pch = 19)
            abline(v = log(cv_fit$lambda.min), col = "blue", lty = 2)
            legend("topright", legend = c("MSE", "Optimal Lambda"), col = c("red", "blue"), pch = 19, lty = 1:2)
            dev.off()
        }
        
        return(list(valid_mse = mse_valid, test_mse = mse_test, coefficients = lasso_coefficients))
        
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
        # Visualization of the process of optimizing lambda
        png(filename = here::here("results", paste0("01_lasso_fit_lambda_", deparse(substitute(y)), "_", deparse(substitute(X)), ".png")))
        plot(cv_fit)
        abline(v = log(cv_fit$lambda.min), col = "blue", lty = 2) # Mark the optimal lambda
        title("Cross-Validation Curve for LASSO")
        dev.off()
        
        # Visualization of how the MSEs evolve with the process
        png(filename = here::here("results", paste0("01_lasso_fit_mse_lambda_", deparse(substitute(y)), "_", deparse(substitute(X)), ".png")))
        plot(log(cv_fit$lambda), cv_fit$cvm, type = "b", pch = 19, col = "red", xlab = "Log(Lambda)", ylab = "Mean Squared Error", main = "MSE vs Log(Lambda)")
        points(log(cv_fit$lambda.min), min(cv_fit$cvm), col = "blue", pch = 19)
        abline(v = log(cv_fit$lambda.min), col = "blue", lty = 2)
        legend("topright", legend = c("MSE", "Optimal Lambda"), col = c("red", "blue"), pch = 19, lty = 1:2)
        dev.off()
        }
        
        return(list(test_mse = mse_test, coefficients = lasso_coefficients))
    }
}

# example usage for run_lasso_regression
result <- run_lasso_regression(X = X_syn, y = y_syn, visualize = TRUE)
