# hodge podge of EDA tasks
# 2024-12-18-wed: 
# - examine percentage of zeros in each feature or sample

# load packages
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(here)
library(gt)
library(Matrix)


# anchor project root directory 
here::i_am("scripts/03_misc_eda_tasks.R")

# load data
vir_lib <- readRDS(here::here("data", "rds", "vir_lib.rds"))
vir_z_pat <- readRDS(here::here("data", "processed", "vir_z_pat.rds"))
vir_z <- readRDS(here::here("data", "rds", "vir_z.rds"))

PAT_META_COLS <- c("rep_id", "Sample", "Library", "IP_reagent", 
                   "COVID19_status", "Hospitalized", "Pull_down_antibody", 
                   "patient", "original_id_column")

# Extract X and y
X_vir_z <- vir_z_pat %>%
    select(-all_of(PAT_META_COLS)) 
X_vir_z_mx <- X_vir_z %>% as.matrix()

y_vir_z <- vir_z_pat %>%
    mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
    pull(COVID19_status)



# TODO: Percentage of zeros 

#==============================================================================
# Computing correlations 
#==============================================================================

# # Note: this part requires additional memory and takes a few min to run 

# # For each column of X_vir_z, calculate the correlation with y_vir_z
# # Suppose X is 844 x 115753, y is length 844
# n <- nrow(X_vir_z_mx)      # 844
# p <- ncol(X_vir_z_mx)      # 115753

# # Center y
# y_mean <- mean(y_vir_z)
# y_c    <- y_vir_z - y_mean
# # Precompute the denominator part for y
# den_y  <- sqrt(sum(y_c^2))  # a single number

# # Center each column of X
# X_means     <- colMeans(X_vir_z_mx)                       # length p
# X_centered  <- sweep(X_vir_z_mx, 2, X_means, FUN = "-")   # same dim as X

# # Numerators: crossprod of X_centered with y_c
# # crossprod(X_centered, y_c) is p x 1
# num <- crossprod(X_centered, y_c)  # dimension: 115753 x 1

# # Denominators: sqrt(colSums(X_centered^2)) * den_y
# # colSums(X_centered^2) is length p
# den_x <- sqrt(colSums(X_centered^2))  # length p
# den   <- den_x * den_y

# # Final correlations: length p
# corrs <- as.vector(num / den)


# # Save the correlation vector to disk in the results directory 
# saveRDS(corrs, file = here::here("results", "pep_corrs.rds"))
# # corrs is now a numeric vector of length 115,753 
# # with cor(X[,j], y) in position j
# print("done running the correlation calculations")

#==============================================================================

# load the correlation vector from disk
# 1) Read in your correlations (a numeric vector)
corrs <- readRDS(file = here::here("results", "pep_corrs.rds"))

# 2) Create corrs_df
corrs_df <- data.frame(
  id = seq_along(corrs),     # or 1:length(corrs)
  correlation = corrs
) %>%
  arrange(desc(correlation)) %>% 
  left_join(
    vir_lib %>% select(id, Organism, `Protein names`), 
    by = "id"
  ) %>%
  drop_na(id, Organism, `Protein names`)

# 3) Create corrs_abs_top20_gt (top 20 by absolute correlation)
corrs_abs_top20_gt <- corrs_df %>%
  arrange(desc(abs(correlation))) %>%
  head(20) %>%
  gt() %>%
  tab_header(
    title = "top 20 peptides by absolute correlation with covid-19 status"
  )

# 4) Or create corrs_top20_gt (top 20 by correlation itself)
corrs_top20_gt <- corrs_df %>%
  arrange(desc(correlation)) %>%
  head(20) %>%
  gt() %>%
  tab_header(
    title = "top 20 peptides by correlation with covid-19 status"
  )

# save corrs_df as csv
# write.csv(corrs_df, file = here::here("results", "pep_sgned_corrs_df.csv"))


#==============================================================================
# two-sample test of means for COVID-19 positive and negative samples
#==============================================================================
# 1) Create a data frame with X_vir_z_mx and y_vir_z
META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "IP_reagent", 
                   "COVID19_status", "Hospitalized", "Pull_down_antibody", 
                   "patient", "original_id_column")
X_y_vir_z <- cbind(as.data.frame(X_vir_z), y_vir_z)
X_vir_z <-vir_z_pat %>%
    select(-all_of(PAT_META_COLS)) 


y_vir_z <- vir_z_pat %>%
    mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
    pull(COVID19_status)

# checking for sparsity
# num_zero  <- sum(X_vir_z_mx == 0)
# total_elems <- length(X_vir_z_mx)
# sparsity_ratio <- num_zero / total_elems

# sparsity_ratio

# Number of nonzero entries:
X_sp <- as(X_vir_z_mx, "dgCMatrix")
nnz <- Matrix::nnzero(X_sp)  
cat("Number of nonzero elements:", nnz, "\n")
cat("Fraction of nonzero elements:", nnz / length(X_vir_z_mx), "\n")

#==============================================================================
