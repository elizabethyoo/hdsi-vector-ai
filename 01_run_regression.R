# Given a m by n matrix dataframe or tibble, where m is the number of rows and n is the number of columns, 
# and a vector of length m, this script runs a LASSO regression and returns the coefficients of the model.

# Install and load the glmnet package for LASSO
if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
library(glmnet)

# Function to split data into train/validation/test or train/test
split_data <- function(X, y, validation = TRUE, train_ratio = 0.6, valid_ratio = 0.2) {
    n <- nrow(X)
    train_idx <- sample(1:n, size = train_ratio * n)
    
    if (validation) {
        valid_idx <- setdiff(sample(1:n, size = (train_ratio + valid_ratio) * n), train_idx)
        test_idx <- setdiff(1:n, c(train_idx, valid_idx))
        
        X_train <- X[train_idx, ]
        y_train <- y[train_idx]
        X_valid <- X[valid_idx, ]
        y_valid <- y[valid_idx]
        X_test <- X[test_idx, ]
        y_test <- y[test_idx]
        
        return(list(X_train = X_train, y_train = y_train, X_valid = X_valid, y_valid = y_valid, X_test = X_test, y_test = y_test))
    } else {
        test_idx <- setdiff(1:n, train_idx)
        
        X_train <- X[train_idx, ]
        y_train <- y[train_idx]
        X_test <- X[test_idx, ]
        y_test <- y[test_idx]
        
        return(list(X_train = X_train, y_train = y_train, X_test = X_test, y_test = y_test))
    }
}

# Example usage
data_split <- split_data(X, y, validation = TRUE)

X_train <- data_split$X_train
y_train <- data_split$y_train
X_test <- data_split$X_test
y_test <- data_split$y_test

if (exists("X_valid")) {
    X_valid <- data_split$X_valid
    y_valid <- data_split$y_valid
}

# Cross-validated LASSO
cv_fit <- cv.glmnet(X_train, y_train, alpha = 1) # alpha=1 for LASSO
best_lambda <- cv_fit$lambda.min

# Predict on validation set if it exists
if (exists("X_valid")) {
    valid_pred <- predict(cv_fit, s = best_lambda, newx = X_valid)
    mse_valid <- mean((y_valid - valid_pred)^2)
    cat("Validation MSE:", mse_valid, "\n")
}

# Predict on test set
test_pred <- predict(cv_fit, s = best_lambda, newx = X_test)
mse_test <- mean((y_test - test_pred)^2)
cat("Test MSE:", mse_test, "\n")

# Coefficients at the best lambda
lasso_coefficients <- coef(cv_fit, s = best_lambda)
print(lasso_coefficients)

# Visualization
# Cross-validation plot for lambda
plot(cv_fit)
abline(v = log(cv_fit$lambda.min), col = "blue", lty = 2) # Mark the optimal lambda
title("Cross-Validation Curve for LASSO")

# Coefficients shrinkage as lambda varies
lasso_fit <- glmnet(X_train, y_train, alpha = 1)
plot(lasso_fit, xvar = "lambda", label = TRUE)
title("Coefficient Paths for LASSO")
