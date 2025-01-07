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
library(stringr)
library(data.table)

# anchor project root directory 
here::i_am("scripts/03_misc_eda_tasks.R")

#============================================================================

# load data
# vir_lib <- readRDS(here::here("data", "rds", "vir_lib.rds"))
# vir_z_pat <- readRDS(here::here("data", "processed", "vir_z_pat.rds"))
vir_z <- readRDS(here::here("data", "rds", "vir_z.rds"))
vir_z_vec <- unlist(vir_z[, 3:ncol(vir_z)], use.names = FALSE)# vir_z z-scores flattened into a single vector -- for all peptide histograms
summary(vir_z_vec)

# PAT_META_COLS <- c("rep_id", "Sample", "Library", "IP_reagent", 
#                    "COVID19_status", "Hospitalized", "Pull_down_antibody", 
#                    "patient", "original_id_column")

# # Extract X and y
# X_vir_z <- vir_z_pat %>%
#     select(-all_of(PAT_META_COLS)) 
# X_vir_z_mx <- X_vir_z %>% as.matrix()

# y_vir_z <- vir_z_pat %>%
#     mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
#     pull(COVID19_status)



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

# # load the correlation vector from disk
# # 1) Read in your correlations (a numeric vector)
# corrs <- readRDS(file = here::here("results", "pep_corrs.rds"))

# # 2) Create corrs_df
# corrs_df <- data.frame(
#   id = seq_along(corrs),     # or 1:length(corrs)
#   correlation = corrs
# ) %>%
#   arrange(desc(correlation)) %>% 
#   left_join(
#     vir_lib %>% select(id, Organism, `Protein names`), 
#     by = "id"
#   ) %>%
#   drop_na(id, Organism, `Protein names`)

# # 3) Create corrs_abs_top20_gt (top 20 by absolute correlation)
# corrs_abs_top20_gt <- corrs_df %>%
#   arrange(desc(abs(correlation))) %>%
#   head(20) %>%
#   gt() %>%
#   tab_header(
#     title = "top 20 peptides by absolute correlation with covid-19 status"
#   )

# # 4) Or create corrs_top20_gt (top 20 by correlation itself)
# corrs_top20_gt <- corrs_df %>%
#   arrange(desc(correlation)) %>%
#   head(20) %>%
#   gt() %>%
#   tab_header(
#     title = "top 20 peptides by correlation with covid-19 status"
#   )

# save corrs_df as csv
# write.csv(corrs_df, file = here::here("results", "pep_sgned_corrs_df.csv"))

# TODO two-sample test of means for COVID-19 positive and negative samples

#==============================================================================
# histogramming z-scores, stratified by COVID-19 status (in order to determine which two-sample test to use)
#==============================================================================

# META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "IP_reagent", 
#                    "Hospitalized", "Pull_down_antibody", 
#                    "patient", "original_id_column") # excludes COVID19_status which will be a part of the final dataframe
# # X_y_vir_z includes z-scores and COVID19_status vs. X_vir_z includes z-scores only
# X_y_vir_z <-vir_z_pat %>%
#     select(-all_of(META_COL_NO_STATUS)) %>%
#     mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
#     # put COVID19_status as the first column
#     select(COVID19_status, everything()) 

# #1. Plot histograms for a random subset of peptides
# # Suppose your full data frame is called df
# # df has 844 rows (patients), 1 column "status", and 115,753 peptide columns.

# # TODO for approaches #1, 2, 4: replace df with X_y_vir_z


# # Randomly pick, say, 5 peptide columns (excluding the 'status' column)
# set.seed(123)
# peptide_cols <- sample(names(X_y_vir_z)[-1], 20)  

# # TODO: put in a helper functions script 
# # Function to get organism by col_name
# get_organism_by_id <- function(id_chr, vir_lib) {
#   id_numeric <- as.numeric(id_chr)
#   organism <- vir_lib %>%
#     filter(id == id_numeric) %>%
#     pull(Organism)
#   return(organism)
# }

# PDF_NAME <- "sampled_zscores_histograms_perpeptide_stratified.pdf"
# pdf(here("results", "eda", "viz", PDF_NAME), width = 10, height = 6)  # opens a multi-page PDF

# # Now for each selected peptide, plot histograms by disease status
# for (col_name in peptide_cols) {
#   covid_z <- X_y_vir_z[X_y_vir_z$COVID19_status == 1, col_name][[1]]
#   control_z <- X_y_vir_z[X_y_vir_z$COVID19_status == 0, col_name][[1]]
  
#   # Plot side-by-side histograms in base R
#   par(mfrow = c(1, 2), mar = c(7, 4, 4, 2) + 1.0)  # 1 row, 2 columns in the plotting window, increase margins
  
#   hist(covid_z,
#      main = "covid+",
#      xlab = "z-score",
#      col = "orange",
#      breaks = 30)

#   # Now add a subtitle manually on a separate line
#   title(sub = paste("id:", col_name, "\norganism:", str_wrap(get_organism_by_id(col_name, vir_lib), width = 50)),
#         line = 6)  

#   hist(control_z,
#     main = "control",
#     xlab = "z-score",
#     col = "blue",
#     breaks = 30)
  
#   title(sub = paste("id:", col_name, "\norganism:", str_wrap(get_organism_by_id(col_name, vir_lib), width = 50)),
#       line = 6) 
#   # Reset plot layout
#   par(mfrow = c(1, 1))
# }
# dev.off() 

#2. Plot histograms by flattening all peptides into one "lump" -- lose info on individual peptides but ah well


# Reshape from wide -> long, excluding 'status'
# This could be large in memory if your dataset is truly massive, 
# because it will create ~ (844 * 115,753) rows!
# use data.table for speed 
# setDT(X_y_vir_z)

# # Melt the data.table
# X_y_vir_z_long <- melt(
#   data     = X_y_vir_z,
#   id.vars  = "COVID19_status",  # columns to keep fixed
#   variable.name = "peptide", 
#   value.name    = "z_score"
# )

# saveRDS(X_y_vir_z_long, file = here::here("data", "processed", "X_y_vir_z_long.rds"))
X_y_vir_z_long <- readRDS(here::here("data", "processed", "X_y_vir_z_long.rds"))


pdf(here("results", "eda", "viz", "all_zscores_allpeptides_histogram.pdf"), width = 10, height = 6) 

ggplot(X_y_vir_z_long, aes(x = z_score, fill = COVID19_status)) +
  geom_histogram(alpha = 0.5, bins = 50) +
  scale_fill_manual(values = c("control" = "blue", "covid+" = "orange")) +
  labs(
    x = "z-score",
    y = "Count",
    fill = "Status",
    title = "Overall Distribution of Z-scores by Disease Status",
    subtitle = "Full Dataset"
  ) +
  theme_minimal()


# sampling_percentages <- seq(0.1, 1.0, by = 0.1)  # 10%, 20%, ..., 100%

# # Then simply partition by status and plot a histogram
# for (p in sampling_percentages) {
#   # Calculate the percentage as a whole number for labeling
#   percentage_label <- paste0(p * 100, "% Sample")
  
#   ### ERASE LATER
#   message(paste("Processing:", percentage_label))

#   # Randomly sample the data
#   sampled_data <- X_y_vir_z_long %>%
#     slice_sample(prop = p)
  
#   ###ERASE  # Check if sampled_data has data
#   if (nrow(sampled_data) == 0) {
#     warning(paste("No data sampled for", percentage_label))
#     next  # Skip to the next iteration
#   }
  
#   # Generate the histogram
#   plot <- ggplot(sampled_data, aes(x = z_score, fill = factor(COVID19_status))) +
#     geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
#     scale_fill_manual(values = c("blue", "orange"), 
#                       labels = c("control", "covid+")) +
#     labs(
#       x = "z-score",
#       y = "Count",
#       fill = "Status",
#       title = paste("Distribution of Z-scores by Disease Status (", percentage_label, ")", sep = ""),
#       subtitle = "Randomly Sampled Data"
#     ) +
#     theme_minimal() +  # Optional: for a cleaner look
#     theme(
#       plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
#       plot.subtitle = element_text(hjust = 0.5, size = 12)
#     )
  
#   # Print the plot to the PDF
#   print("Generating plot...")
#   print(plot)
#   print("Plot generated.")
# }

dev.off()



# #3. Plot histograms for the top 20 peptides by correlation with COVID-19 status
# # TODO 

# #4. Plot summary statistics instead of full histograms 


# library(dplyr)

# # Summarize mean z-score per (status, peptide)
# summaries <- df %>%
#   group_by(status) %>%
#   summarize(across(-status, 
#                    list(mean = ~mean(.x, na.rm = TRUE)), 
#                    .names = "{.col}_{.fn}"))







# #==============================================================================
# # miscellanea 
# #==============================================================================
# # checking for sparsity -- this crashes if you run on the head node
# # num_zero  <- sum(X_vir_z_mx == 0)
# # total_elems <- length(X_vir_z_mx)
# # sparsity_ratio <- num_zero / total_elems
# # sparsity_ratio