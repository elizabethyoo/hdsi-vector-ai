###############################################################################
## 04_eda_randomforest_combined.R
## - Combined script that merges:
##    (A) Single train/test approach with feature importance (from 04_eda_randomforest.R)
##    (B) 10-fold cross-validation logic with mean ROC + std dev ribbons (from 04_2_eda_cv_randomfoest.R)
###############################################################################

# 1) Load Libraries ---------------------------------------------------------
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
  randomForestExplainer,
  viridis,
  caret,
  purrr,       # for map(), map_chr() if needed
  RColorBrewer
)

# ----------------------------------------------------------------------------
# 2) Anchor the project root directory with {here}
here::i_am("scripts/04_eda_randomforest_combined.R")

# ----------------------------------------------------------------------------
# 3) Helper Functions
#    - get_organism_by_id() from your original script
#    - You can add or remove any others as needed
get_organism_by_id <- function(id_chr, cov_lib, vir_lib) {
  # If the id comes from covidscan data (prefix "c_")
  if (startsWith(id_chr, "c_")) {
    id_numeric <- as.numeric(sub("^c_", "", id_chr))
    row_data <- cov_lib %>% filter(id == id_numeric)
    
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
    # If the id comes from virscan data (prefix "v_")
    id_numeric <- as.numeric(sub("^v_", "", id_chr))
    row_data <- vir_lib %>% filter(id == id_numeric)
    
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

# ----------------------------------------------------------------------------
# 4) Parameters & Output Filenames
NUM_TEST_PEPTIDES <- ""  # e.g., "100" for testing subset or "" for full
RF_RESULTS_FNAME <- "rf-result"
DATASET_NAME <- paste0("covscan", "_", NUM_TEST_PEPTIDES)
BASE_FNAME <- paste0(RF_RESULTS_FNAME, "_", DATASET_NAME)
EXCEL_FNAME <- paste0(BASE_FNAME, "_", DATASET_NAME, ".xlsx")

# Create results directory, unique for this run
RUN_DIR <- here::here("results", "eda", "rf", paste0(BASE_FNAME, "_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))
dir.create(RUN_DIR, recursive = TRUE)

# ----------------------------------------------------------------------------
# 5) Data Loading & Preparation
cat("Step: Loading and preparing data...\n")

# Example: Covscan Data (Adjust as needed for your actual dataset selection)
cov_z_pat <- readRDS(here::here("data", "processed", "cov_z_pat.rds"))
com_z_pat <- readRDS(here::here("data", "processed", "Shrock_vir_cov_combined" "combined_z_pat.rds"))
cov_lib   <- readRDS(here::here("data", "processed", "cov_lib.rds"))
vir_lib   <- readRDS(here::here("data", "rds", "vir_lib.rds"))  # if needed

# Remove meta columns, turn COVID19_status into factor, etc.
META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "Hospitalized")
X_y_df <- cov_z_pat %>%
  select(-all_of(META_COL_NO_STATUS)) %>%
  mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
  select(COVID19_status, everything()) %>%
  mutate(COVID19_status = as.factor(COVID19_status))

# (Optional) Test on a Subset
if (NUM_TEST_PEPTIDES != "") {
  TEST_COL_FNAME <- paste0(BASE_FNAME, "_", "test_columns.rds")
  test_peptide_columns <- sample(colnames(X_y_df[,-1]), as.numeric(NUM_TEST_PEPTIDES))
  saveRDS(test_peptide_columns, here::here("results", "eda", "rf", TEST_COL_FNAME))
  X_y_df <- X_y_df %>% select(all_of(c("COVID19_status", test_peptide_columns)))
}

cat("Data dimensions:", dim(X_y_df), "\n")

# ----------------------------------------------------------------------------
# 6) Define Functions for Random Forest Approaches

# 6A) Single Train/Test Split ------------------------------------------------
rf_train_test <- function(X_y, base_fname, run_dir, cov_lib, vir_lib) {
  # X_y is a dataframe with first column = COVID19_status
  # base_fname, run_dir for file saving
  # cov_lib, vir_lib for get_organism_by_id (feature importance)
  
  # Separate features and target
  X <- X_y %>% select(-COVID19_status)
  y <- X_y$COVID19_status
  
  cat(">>> Train/Test Split...\n")
  set.seed(618)
  train_idx <- sample(seq_len(nrow(X)), size = 0.8 * nrow(X))
  X_train   <- X[train_idx, ]
  y_train   <- y[train_idx]
  X_test    <- X[-train_idx, ]
  y_test    <- y[-train_idx]
  
  cat("Train Set Dimensions:", dim(X_train), "\n")
  cat("Test Set Dimensions:",  dim(X_test), "\n")
  cat("Number of positive cases (train):", sum(y_train == "1"), "\n")
  cat("Number of positive cases (test): ", sum(y_test  == "1"), "\n")
  
  # Save the train/test split
  TR_TE_SPLIT_FNAME <- paste0(base_fname, "_train_test_data.rds")
  saveRDS(list(X_train = X_train, y_train = y_train, X_test = X_test, y_test = y_test),
          file.path(run_dir, TR_TE_SPLIT_FNAME))
  
  cat(">>> Fitting Random Forest (train/test)\n")
  rf_model <- randomForest(X_train, y_train, importance = TRUE)
  print(rf_model)
  
  # Predictions
  class_predictions <- predict(rf_model, X_test, type = "response")
  prob_predictions  <- predict(rf_model, X_test, type = "prob")[,2]
  accuracy          <- mean(class_predictions == y_test)
  
  # Save model and predictions
  RF_MODEL_FNAME      <- paste0(base_fname, "_model-fit.rds")
  RF_CLASS_PRED_FNAME <- paste0(base_fname, "_class-predictions.rds")
  RF_PROB_PRED_FNAME  <- paste0(base_fname, "_prob-predictions.rds")
  RF_ACC_FNAME        <- paste0(base_fname, "_accuracy.rds")
  
  saveRDS(rf_model,        file.path(run_dir, RF_MODEL_FNAME))
  saveRDS(class_predictions, file.path(run_dir, RF_CLASS_PRED_FNAME))
  saveRDS(prob_predictions,  file.path(run_dir, RF_PROB_PRED_FNAME))
  saveRDS(accuracy,          file.path(run_dir, RF_ACC_FNAME))
  
  # Compute ROC / AUC
  pred_test <- ROCR::prediction(prob_predictions, as.numeric(y_test) - 1)
  perf_test <- ROCR::performance(pred_test, "tpr", "fpr")
  auc_test  <- ROCR::performance(pred_test, "auc")@y.values[[1]]
  
  # Save ROC objects
  saveRDS(list(pred = pred_test, perf = perf_test, auc = auc_test),
          file.path(run_dir, paste0(base_fname, "_rf_pred_perf_auc.rds")))
  
  # Generate and save ROC curve
  roc_data <- data.frame(
    fpr = unlist(perf_test@x.values),
    tpr = unlist(perf_test@y.values)
  )
  roc_plot <- ggplot(roc_data, aes(x = fpr, y = tpr)) +
    geom_line(color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    labs(
      title = "ROC Curve (Train/Test)",
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    annotate("text", x = 0.75, y = 0.25, label = sprintf("AUC = %.2f", auc_test), color = "blue") +
    theme_minimal()
  
  ROC_FNAME <- paste0(base_fname, "_roc_curve.png")
  ggsave(file.path(run_dir, ROC_FNAME), plot = roc_plot, width = 10, height = 6, dpi = 300)
  
  # Return model & data for further usage (feature importance, etc.)
  list(
    model          = rf_model,
    X_test         = X_test,
    y_test         = y_test,
    prob_predictions = prob_predictions,
    accuracy       = accuracy,
    auc_test       = auc_test
  )
}

# 6B) Cross-Validation (N-fold) ----------------------------------------------
rf_cross_validation <- function(X_y, base_fname, run_dir, num_folds = 10, cov_lib = NULL, vir_lib = NULL) {
  # X_y is a dataframe with first column = COVID19_status
  # base_fname, run_dir for file saving
  # num_folds: number of cross-validation folds
  # cov_lib, vir_lib for annotation (optional)
  
  cat(">>> Cross-Validation (", num_folds, "-fold )...\n", sep = "")
  
  # Separate features and target
  X <- X_y %>% dplyr::select(-COVID19_status)
  y <- X_y$COVID19_status
  
  set.seed(618)
  folds <- caret::createFolds(y, k = num_folds, list = TRUE, returnTrain = TRUE)
  
  roc_data_list    <- list()         # Store ROC points
  auc_values       <- numeric()      # Store fold-wise AUC
  importance_list  <- list()         # NEW: Store fold-wise feature importance
  
  # Loop over folds
  for (i in seq_along(folds)) {
    cat("Processing Fold", i, "\n")
    train_idx <- folds[[i]]
    
    X_train <- X[train_idx, ]
    y_train <- y[train_idx]
    X_test  <- X[-train_idx, ]
    y_test  <- y[-train_idx]
    
    # Train Random Forest
    rf_model <- randomForest::randomForest(X_train, y_train, importance = TRUE)
    
    # 1) Store Feature Importance for this fold
    fold_importance <- as.data.frame(randomForest::importance(rf_model))
    fold_importance$Feature <- rownames(fold_importance)
    fold_importance$Fold    <- i
    importance_list[[i]]    <- fold_importance
    
    # 2) Predicted probabilities for ROC
    prob_predictions <- predict(rf_model, X_test, type = "prob")[,2]
    
    # 3) ROC/AUC
    pred <- ROCR::prediction(prob_predictions, as.numeric(y_test) - 1)
    perf <- ROCR::performance(pred, "tpr", "fpr")
    auc  <- ROCR::performance(pred, "auc")@y.values[[1]]
    auc_values <- c(auc_values, auc)
    
    roc_data_list[[i]] <- data.frame(
      fold = i,
      fpr  = unlist(perf@x.values),
      tpr  = unlist(perf@y.values)
    )
  }
  
  # ---------------------------
  # (A) Combine & Aggregate Feature Importance
  # ---------------------------
  
  # Merge all fold importances into one data frame
  all_importance <- dplyr::bind_rows(importance_list)
  
  # For example, if you'd like to average 'MeanDecreaseGini' across folds:
  # (Adjust to your chosen metric; e.g., 'IncNodePurity' or others)
  importance_summary <- all_importance %>%
    dplyr::group_by(Feature) %>%
    dplyr::summarise(
      MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE),
      sd_MDg           = sd(MeanDecreaseGini, na.rm = TRUE),
      # Add any other summary stats if desired
      folds_included   = n()  # how many folds had that feature (should be all folds)
    ) %>%
    dplyr::arrange(desc(MeanDecreaseGini))
  
  # (Optional) Annotate peptides with organism information
  # if your 'Feature' values match what get_organism_by_id expects:
  # 
  # importance_summary <- importance_summary %>%
  #   mutate(attributes = map(Feature, ~ get_organism_by_id(.x, cov_lib, vir_lib)),
  #          organism   = map_chr(attributes, "Organism")) %>%
  #   select(-attributes)
  
  # Save this aggregated importance data
  importance_csv <- paste0(base_fname, "_cv_importance_summary.csv")
  readr::write_csv(importance_summary, file.path(run_dir, importance_csv))
  
  # ---------------------------
  # (B) Combine & Aggregate ROC Curves
  # ---------------------------
  
  roc_data <- do.call(rbind, roc_data_list)
  
  # Compute mean +/- SD TPR at fixed FPR points
  fpr_values <- seq(0, 1, by = 0.01)
  tpr_matrix <- sapply(roc_data_list, function(df) {
    approx(df$fpr, df$tpr, xout = fpr_values)$y
  })
  mean_tpr <- rowMeans(tpr_matrix, na.rm = TRUE)
  sd_tpr   <- apply(tpr_matrix, 1, sd, na.rm = TRUE)
  
  roc_summary <- data.frame(
    fpr       = fpr_values,
    mean_tpr  = mean_tpr,
    lower_tpr = mean_tpr - sd_tpr,
    upper_tpr = mean_tpr + sd_tpr
  )
  
  # Plot CV ROC
  roc_plot <- ggplot() +
    geom_line(data = roc_data, aes(x = fpr, y = tpr, group = fold), 
              color = "lightblue", alpha = 0.7) +
    geom_line(data = roc_summary, aes(x = fpr, y = mean_tpr), 
              color = "blue", size = 1.2) +
    geom_ribbon(data = roc_summary, aes(x = fpr, ymin = lower_tpr, ymax = upper_tpr), 
                fill = "gray", alpha = 0.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    labs(
      title = paste0(num_folds, "-Fold Cross-Validation ROC"),
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    theme_minimal()
  
  # Save AUC values and the final plot
  AUC_FNAME <- paste0(base_fname, "_cv_auc_values.csv")
  readr::write_csv(
    data.frame(Fold = seq_along(auc_values), AUC = auc_values), 
    file.path(run_dir, AUC_FNAME)
  )
  
  ROC_PLOT_FNAME <- paste0(base_fname, "_cv_roc_plot.png")
  ggsave(file.path(run_dir, ROC_PLOT_FNAME), plot = roc_plot, width = 10, height = 6, dpi = 300)
  
  cat("Cross-validation complete.\n")
  cat("Feature importance saved to: ", file.path(run_dir, importance_csv), "\n")
  cat("AUC values saved to:         ", file.path(run_dir, AUC_FNAME), "\n")
  cat("ROC plot saved to:           ", file.path(run_dir, ROC_PLOT_FNAME), "\n")
  
  # Return summary
  list(
    roc_data_list      = roc_data_list,
    roc_summary        = roc_summary,
    auc_values         = auc_values,
    importance_summary = importance_summary  # aggregated feature importance
  )
}

# ----------------------------------------------------------------------------
# 7) Run the Workflow
#    (A) Single train/test approach
#    (B) Cross-validation approach
#    Then (C) Feature importance from the single train/test model
# ----------------------------------------------------------------------------

### 7A) Single Train/Test
train_test_res <- rf_train_test(
  X_y = X_y_df, 
  base_fname = BASE_FNAME,
  run_dir = RUN_DIR,
  cov_lib = cov_lib,
  vir_lib = vir_lib
)

### 7B) Cross-Validation
cv_res <- rf_cross_validation(
  X_y = X_y_df,
  base_fname = BASE_FNAME,
  run_dir = RUN_DIR,
  num_folds = 10  # e.g., 10-fold
)

# ----------------------------------------------------------------------------
# 8) Feature Importance & Plotting (Train/Test Model) ------------------------
cat("Step: Feature Importance (Train/Test Model)\n")

rf_model <- train_test_res$model

# Extract raw importance
importance_data <- as.data.frame(randomForest::importance(rf_model))
importance_data$peptide <- rownames(importance_data)

# measure_importance() from randomForestExplainer
importance_frame <- measure_importance(rf_model) %>%
  mutate(peptide = as.character(variable))

# Merge p-values
importance_data <- importance_data %>%
  left_join(importance_frame %>% select(peptide, p_value), by = "peptide")

# How many features to show?
TOP_NUM <- min(20, nrow(importance_data))

# Detailed annotation
importance_data_topN <- importance_data %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice_head(n = TOP_NUM) %>%
  mutate(
    attributes = map(peptide, ~ get_organism_by_id(.x, cov_lib, vir_lib)),
    organism   = map_chr(attributes, "Organism"),
    Protein    = map_chr(attributes, "Protein"),
    Sequence   = map_chr(attributes, "Sequence"),
    Start      = map_dbl(attributes, "Start"),
    End        = map_dbl(attributes, "End")
  ) %>%
  select(-attributes) %>%
  filter(!is.na(organism)) %>%
  distinct(organism, .keep_all = TRUE)

# Write top importance to CSV
IMP_FNAME <- paste0(BASE_FNAME, "_importance_data_top", TOP_NUM, ".csv")
write.csv(importance_data_topN, file.path(RUN_DIR, IMP_FNAME))

cat("Top features saved to:", IMP_FNAME, "\n")
cat("Feature Importance dimension:", dim(importance_data_topN), "\n")

# Filter & rank by p-value
ranked_p_values <- importance_data_topN %>%
  filter(!is.na(p_value)) %>%
  filter(p_value < 0.05) %>%
  arrange(p_value) %>%
  distinct(organism, .keep_all = TRUE) %>%
  head(TOP_NUM)

RANKED_P_FNAME <- paste0(BASE_FNAME, "_ranked_p_values.csv")
write.csv(ranked_p_values, file.path(RUN_DIR, RANKED_P_FNAME))

# ----------------------------------------------------------------------------
# 9) Plot Feature Importance (Top 20, Train/Test) ---------------------------
cat("Step: Plotting Feature Importance...\n")

# Construct a label for plotting if desired
importance_data_topN <- importance_data_topN %>%
  mutate(
    pep_label = paste(peptide, organism, sep = " | ")
  )

# Choose a palette
importance_palette <- viridis::viridis(n_distinct(importance_data_topN$organism))

# Bar plot of MeanDecreaseGini
p_imp <- ggplot(importance_data_topN, 
                aes(x = reorder(pep_label, MeanDecreaseGini), 
                    y = MeanDecreaseGini, 
                    fill = organism)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = importance_palette, name = "Organism") +
  theme_minimal() +
  labs(
    title = "Random Forest - Feature Importance",
    x = "Peptide (Organism)",
    y = "Mean Decrease Gini"
  ) +
  theme(
    legend.position = "right",
    axis.text.y     = element_text(size = 8),
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5)
  )

IMP_PLOT_FNAME <- paste0(BASE_FNAME, "_feature_importance.png")
ggsave(file.path(RUN_DIR, IMP_PLOT_FNAME), plot = p_imp, width = 10, height = 6, dpi = 300)

cat("Feature importance plot saved to:", IMP_PLOT_FNAME, "\n")

# ----------------------------------------------------------------------------
# 10) Done
cat("Script complete! Check outputs in:", RUN_DIR, "\n")
