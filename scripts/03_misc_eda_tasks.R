# hodge podge of EDA tasks

# CHECKLIST BEFORE RUNNING SCRIPT 
# 1. Check #Writing to file section 
# 2. Check #covscan and #virscan sections -- comment out the one not being used
  # 2-1. Use get_cov_organism_by_id, cov_lib for covscan data and get_organism_by_id, vir_lib for virscan data
# 3. If using subset data, check #TEST ON A SUBSET section



if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tibble,
  tidyr,
  dplyr,
  ggplot2,
  reshape2,
  here,
  gt,
  Matrix,
  stringr,
  data.table,
  purrr,
  broom,
  parallel,
  gridExtra,      # For arranging multiple plots/tables
  ggrepel,        # For enhanced plot labeling
  openxlsx,
  grid,
  pbapply,        # For parallel processing with progress bar
  RColorBrewer,
  patchwork,
  viridis
)
# anchor project root directory 
here::i_am("scripts/03_misc_eda_tasks.R")

#============================================================================

# Define number of peptides for testing - default empty string
NUM_TEST_PEPTIDES <- 500

#Writing to file ===========================================================
# CHECK BEFORE RUNNING - CHANGE TO APPROPRIATE NAMES ####################################################################
TT_RESULTS_FNAME <- "tt-result"
DATASET_NAME <- paste0("covscan", "_", NUM_TEST_PEPTIDES)
BASE_FNAME <- paste0(TT_RESULTS_FNAME, "_", DATASET_NAME)
EXCEL_FNAME <- paste0(BASE_FNAME, "_", DATASET_NAME, ".xlsx")

# Create directory for this run
RUN_DIR <- here::here("results", "eda", "ttest", paste0(BASE_FNAME, "_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))
dir.create(RUN_DIR, recursive = TRUE)


# load data
vir_z_pat <- readRDS(here::here("data", "processed", "vir_z_pat.rds"))
vir_lib <- readRDS(here::here("data", "rds", "vir_lib.rds"))
vir_z <- readRDS(here::here("data", "rds", "vir_z.rds"))

cov_z_pat <- readRDS(here::here("data", "processed", "cov_z_pat.rds"))
cov_lib <- readRDS(here::here("data", "processed", "cov_lib.rds"))
# BOOKMARK - sanity check file names, and checklist

# vir_z_vec <- unlist(vir_z[, 3:ncol(vir_z)], use.names = FALSE)# vir_z z-scores flattened into a single vector -- for all peptide histograms
# summary(vir_z_vec)

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

# this helper function is for covscan; for now comment out when doing stuff with virscan data
get_cov_organism_by_id <- function(id_chr, cov_lib) {
  id_numeric <- as.numeric(id_chr)
  
  # Filter the row corresponding to the given id
  row_data <- cov_lib %>% 
    filter(id == id_numeric)
  
  # Check if row_data is empty
  if (nrow(row_data) == 0) {
    return(list(
      Organism = NA,
      Protein = NA,
      Sequence = NA,
      Start = NA,
      End = NA
    ))  # Return a list with NAs
  }
  
  # Extract attributes
  organism <- row_data %>% pull(Organism) %>% first()
  protein <- row_data %>% pull(Protein) %>% first()
  sequence <- row_data %>% pull(Sequence) %>% first()
  start_pos <- row_data %>% pull(start) %>% first()
  end_pos <- row_data %>% pull(end) %>% first()
  
  # Return attributes as a named list
  return(list(
    Organism = organism,
    Protein = protein,
    Sequence = sequence,
    Start = start_pos,
    End = end_pos
  ))
}


# # TODO: put in a helper functions script 
# # Function to get organism by col_name
# TODO: if organism is missing, get protein name, and then species, and then sequence
get_organism_by_id <- function(id_chr, vir_lib) {
  id_numeric <- as.numeric(id_chr)
  
  # Filter the row corresponding to the given id
  row_data <- vir_lib %>% 
    filter(id == id_numeric)
  
  # Check if row_data is empty
  if (nrow(row_data) == 0) {
    return(NA)  # Return NA if no matching id is found
  }
  
  # Try to get the organism
  label <- row_data %>%
    pull(Organism) %>%
    first()
  
  # Fallbacks if Organism is NA
  if (is.na(label)) {
    protein_name <- row_data %>%
      pull(`Protein names`) %>%
      first()
    if (!is.na(protein_name)) {
      label <- paste("protein:", protein_name)
    } else {
      species <- row_data %>%
        pull(Species) %>%
        first()
      if (!is.na(species)) {
        label <- paste("species:", species)
      } else {
        sequence <- row_data %>%
          pull(Sequence) %>%
          first()
        label <- paste("sequence:", sequence)
      }
    }
  }
  
  return(label)
}


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
# X_y_vir_z_long <- readRDS(here::here("data", "processed", "X_y_vir_z_long.rds"))


# pdf(here("results", "eda", "viz", "all_zscores_allpeptides_histogram.pdf"), width = 10, height = 6) 

# ggplot(X_y_vir_z_long, aes(x = z_score, fill = COVID19_status)) +
#   geom_histogram(alpha = 0.5, bins = 50) +
#   scale_fill_manual(values = c("control" = "blue", "covid+" = "orange")) +
#   labs(
#     x = "z-score",
#     y = "Count",
#     fill = "Status",
#     title = "Overall Distribution of Z-scores by Disease Status",
#     subtitle = "Full Dataset"
#   ) +
#   theme_minimal()


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

# dev.off()



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

#==============================================================================
# TODO two-sample test of means for COVID-19 positive and negative samples
#==============================================================================

# # rows of virlib dataframe where any of the columns contain NA
# vir_lib_na_rows <- vir_lib %>% 
#   filter(if_any(Sequence, is.na))
# vir_lib_na_rows_lab <- vir_lib %>% 
#   filter(if_any(c(Organism, `Protein names`, Species), is.na)) # note: some of these labels may be missing but sequence data is available for all entries of virlib -- not too big of a deal if a label is missing

# #covscan ####################################################################
META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "Hospitalized") # no COVID19_status

X_y_cov_z <- cov_z_pat %>%
  select(-all_of(META_COL_NO_STATUS)) %>%
  mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
  # put COVID19_status as the first column
  select(COVID19_status, everything()) %>% 
  mutate(COVID19_status = as.factor(COVID19_status)) # convert to factor

cat("Step: Data Preparation\n")
cat("Dimensions of X_y_cov_z:", dim(X_y_cov_z), "\n")
print(head(X_y_cov_z))  # Inspect a small portion

X_y_cov_z_dt <- as.data.table(X_y_cov_z)
X_y_dt <- X_y_cov_z_dt 
X_dt <- X_y_cov_z_dt %>% select(-COVID19_status)

#TEST ON A SUBSET 
Randomly select peptides
TEST_COL_FNAME <- paste0(BASE_FNAME, "_", "test_columns.rds")
test_peptide_columns <- sample(colnames(X_y_cov_z[,-1]), NUM_TEST_PEPTIDES)
saveRDS(test_peptide_columns, here::here(RUN_DIR, TEST_COL_FNAME))
X_y_cov_z <- X_y_cov_z %>% select(all_of(c("COVID19_status", test_peptide_columns))) 

X_cov_z <- X_y_cov_z %>%
    select(-COVID19_status)

y_cov_z <- X_y_cov_z$COVID19_status

# Confirm dimensions
cat("Peptide Data Dimensions:", dim(X_cov_z), "\n")
cat("COVID19_status Dimensions:", length(y_cov_z), "\n")
# #covscan ####################################################################


# # # #virscan ####################################################################
# ### data prep
# META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "IP_reagent", 
#                    "Hospitalized", "Pull_down_antibody", 
#                    "patient", "original_id_column") # excludes COVID19_status which will be a part of the final dataframe


# # X_y_vir_z includes z-scores and COVID19_status vs. X_vir_z includes z-scores only
# X_y_vir_z <- vir_z_pat %>%
#     select(-all_of(META_COL_NO_STATUS)) %>%
#     mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
#     # put COVID19_status as the first column
#     select(COVID19_status, everything()) %>% #factorize covid status
#     mutate(COVID19_status = as.factor(COVID19_status))
# X_y_vir_z_dt <- as.data.table(X_y_vir_z)


# X_vir_z <- X_y_vir_z %>%
#     select(-COVID19_status)
# X_vir_z_dt <- as.data.table(X_vir_z)


# # Confirm dimensions
# cat("Peptide Data Dimensions:", dim(X_vir_z_dt), "\n")
# # # #virscan ####################################################################


# Define an output directory for test plots and tables
output_dir <- RUN_DIR

# Define the PDF and Excel file paths for test
pdf_file <- paste0(here(output_dir, BASE_FNAME),".pdf")
excel_file <- paste0(here(output_dir, BASE_FNAME),".xlsx")

# Set seed for reproducibility
set.seed(123)

# dataset-neutral input matrices


# Define number of peptides for testing
num_test_peptides <- NUM_TEST_PEPTIDES

# Randomly select peptides
test_peptide_columns <- sample(colnames(X_dt), num_test_peptides)

# Create the subset dataframe
X_y <- X_y_dt  # %>% select(all_of(c("COVID19_status", test_peptide_columns))) # uncomment for testing on small subset


# Extract peptide data by excluding disease status column
peptide_data_subset <- X_y %>%
  select(-COVID19_status)
disease_status_subset <- X_y$COVID19_status
# Confirm dimensions
cat("Peptide Data Subset Dimensions:", dim(peptide_data_subset), "\n")

# Define the number of top peptides to visualize (adjust as needed)
top_n <- 5 # Reduced for the subset

# Remove peptides with zero variance or any NA values to ensure valid t-tests
peptide_data_subset <- peptide_data_subset %>%
  select(where(~ var(., na.rm = TRUE) > 0)) %>%  # Corrected usage
  select(where(~ !any(is.na(.))))

# Define t-test function with error handling (serial processing for small subset)
t_test_function <- function(peptide_values) {
  test_result <- try(t.test(peptide_values ~ disease_status_subset, var.equal = FALSE), silent = TRUE)
  
  if(class(test_result) == "try-error") {
    return(tibble(statistic = NA, p.value = NA, parameter = NA, conf.low = NA, conf.high = NA, estimate1 = NA, estimate2 = NA))
  } else {
    return(tibble(statistic = test_result$statistic,
                  p.value = test_result$p.value,
                  parameter = test_result$parameter,
                  conf.low = test_result$conf.int[1],
                  conf.high = test_result$conf.int[2],
                  estimate1 = test_result$estimate[1],
                  estimate2 = test_result$estimate[2]))
  }
}

# Perform t-tests using apply (serially)
t_tests_subset <- pbapply(peptide_data_subset, 2, t_test_function)


# Function to get the first N words of a string
get_first_n_words <- function(text, n = 10) {
  words <- unlist(strsplit(text, "\\s+"))  # Split the text into words
  paste(head(words, n), collapse = " ")  # Combine the first N words
}

# Combine results into a dataframe
results_unequal_subset <- bind_rows(t_tests_subset, .id = "peptide") %>%  
# add a column Organism using get_organism_by_id function # CHANGE TO get_cov_organism_by_id for covscan data
  mutate(
  attributes = map(peptide, get_cov_organism_by_id, cov_lib),
  organism = map_chr(attributes, "Organism"),
  Protein = map_chr(attributes, "Protein"),
  Sequence = map_chr(attributes, "Sequence"),
  Start = map_dbl(attributes, "Start"),
  End = map_dbl(attributes, "End")
  )%>%
  select(-attributes) %>% # Remove the intermediate list column
  # if any of the organism values is NA, drop the rows
  filter(!is.na(organism)) %>%
  # rank by p-value
  arrange(p.value) %>%
  mutate(
  first_n_words = sapply(organism, get_first_n_words, n = 10),  # Get the first N words of `organism`
  first_char_protein = substr(Protein, 1, 1),                  # Get the first character of `Protein`
  start_end = paste0(Start, "-", End),                        # Concatenate Start and End with "-"
  PEP_FULL_LABEL = paste(first_n_words, first_char_protein, start_end, sep = ", ")  # Combine all components
) %>%
select(-first_n_words, -first_char_protein, -start_end)  # Optionally drop intermediate columns
saveRDS(results_unequal_subset, here::here("results", "eda", "ttest", paste0(BASE_FNAME, "_", "results_nopval.rds")))                          

# Adjust p-values for multiple testing using p.adjust (no package needed)
results_unequal_subset <- results_unequal_subset %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH"))

# Rank peptides by adjusted p-value
results_unequal_subset <- results_unequal_subset %>%
  arrange(p.adjust) %>%
  mutate(rank_p = row_number())

# Rank peptides by absolute t-statistic
results_unequal_subset <- results_unequal_subset %>%
  mutate(abs_t = abs(statistic)) %>%
  arrange(desc(abs_t)) %>%
  mutate(rank_t = row_number())

# Add significance information before selecting top peptides
# Add significance information and process `results_unequal_subset`
results_unequal_subset <- results_unequal_subset %>%
  mutate(significant = ifelse(p.adjust < 0.05, "Yes", "No")) %>%
  mutate(color_group = case_when(
    significant == "Yes" & estimate1 > estimate2 ~ "control",
    significant == "Yes" & estimate1 < estimate2 ~ "covid+",
    TRUE ~ "not significant"
  )) %>%
  mutate(
    PEP_FULL_LABEL = stringr::str_wrap(PEP_FULL_LABEL, width = 30)
  )

# Select top peptides by p-value
top_peptides_p_subset <- results_unequal_subset %>%
  filter(!is.na(p.adjust) & significant == "Yes") %>%  # Exclude NAs and filter for significant peptides
  arrange(p.adjust) %>%
  slice(1:top_n) %>%
  mutate(
    PEP_FULL_LABEL = stringr::str_wrap(PEP_FULL_LABEL, width = 30),
    PEP_FULL_LABEL = factor(PEP_FULL_LABEL, levels = unique(PEP_FULL_LABEL[order(p.adjust)]))
  )

# Define the organism color palette
color_palette <- viridis(n_distinct(top_peptides_p_subset$PEP_FULL_LABEL))

# Map organism levels to colors
organism_colors <- setNames(color_palette, levels(top_peptides_p_subset$PEP_FULL_LABEL))

# Print organism levels and color mapping
print(levels(top_peptides_p_subset$PEP_FULL_LABEL))
print(organism_colors)


# Define top peptides for labeling (adjust number as needed)
num_labels <- min(3, nrow(top_peptides_p_subset))  # Adjust as needed
top_peptides_labels_subset <- top_peptides_p_subset %>%
  arrange(p.adjust) %>%
  slice(1:num_labels)

# Prepare to generate and save plots into a single PDF
pdf(file = pdf_file, width = 11, height = 8.5)  # Landscape

# Print messages for debugging before plot
# Check top peptides data
print(top_peptides_p_subset)
print(unique(top_peptides_p_subset$PEP_FULL_LABEL))
# check all peptides data
print(results_unequal_subset)
print(unique(results_unequal_subset$color_group))
# check color mapping
print(organism_colors)


# 1. Volcano Plot

# Update volcano plot
volcano_plot_subset <- ggplot(results_unequal_subset, 
                              aes(x = statistic, 
                                  y = -log10(p.adjust))) +
  # Default points for all peptides (significance categories)
  geom_point(aes(color = color_group), alpha = 0.8) +
  
  # Points for top 20 peptides with organism-specific colors
  geom_point(data = top_peptides_p_subset, 
             aes(color = PEP_FULL_LABEL), 
             size = 3) +
  
  # Text labels for top 20 peptides with organism-specific colors
  geom_text_repel(data = top_peptides_p_subset, 
                  aes(label = peptide, color = PEP_FULL_LABEL), 
                  size = 3, 
                  max.overlaps = Inf, 
                  force = 2, 
                  nudge_y = 0.5, 
                  segment.color = "grey50") +
  
  # Update color scales for both significance groups and organisms
  scale_color_manual(
    values = c("control" = "blue", "covid+" = "orange", "not significant" = "grey", organism_colors),
    name = "Significance/Organism",
    breaks = c("control", "covid+", "not significant", levels(top_peptides_p_subset$PEP_FULL_LABEL))
  ) +
  
  # Labels and theme adjustments
  labs(
    title = "Volcano Plot (Unequal Variances) - All Peptides",
    x = "t-statistic",
    y = "-log10(Adjusted p-value)",
    color = "Significance/Organism"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",          # Position the legend on the right
    legend.text = element_text(size = 8),  # Adjust legend text size
    plot.margin = unit(c(1, 1, 1, 1), "cm"),  # Adjust plot margins
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)  # Center and bold the title
  )

# Save the volcano plot to the same PDF
print(volcano_plot_subset)




# 2. Ranking Plot
# Generate a larger palette
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(n_distinct(top_peptides_p_subset$PEP_FULL_LABEL))

# Wrap long legend text and order the legend
top_peptides_p_subset <- top_peptides_p_subset %>%
  mutate(
    PEP_FULL_LABEL = stringr::str_wrap(PEP_FULL_LABEL, width = 30),
    PEP_FULL_LABEL = factor(PEP_FULL_LABEL, levels = unique(PEP_FULL_LABEL[order(p.adjust)]))
  )

# Create the ranking plot
ranking_plot_subset <- ggplot(top_peptides_p_subset, 
                         aes(x = reorder(peptide, -p.adjust), 
                           y = -log10(p.adjust), 
                           fill = PEP_FULL_LABEL)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = color_palette, name = "Organism") +
  labs(title = "Top Peptides by Adjusted p-value - All Peptides",
       x = NULL, 
       y = "-log10(Adjusted p-value)") +
  theme_minimal() +
  theme(
    legend.position = "right",  # Legend on the right
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Ensure centered and visible title
    plot.margin = unit(c(1.5, 1, 1, 1), "cm")  # Adjust margins
  )

# Print the plot directly
print(ranking_plot_subset)

# # 3. Heatmap
# # Extract data for the top peptides
# heatmap_data_subset <- peptide_data_subset[, top_peptides_p_subset$peptide, drop = FALSE]

# # Check if all selected columns are numeric
# if(!all(sapply(heatmap_data_subset, is.numeric))) {
#   stop("Not all selected peptide columns are numeric. Please convert non-numeric columns to numeric.")
# }

# # Scale the data
# heatmap_scaled_subset <- t(scale(t(as.matrix(as.numeric(heatmap_data_subset))), center = TRUE, scale = TRUE))

# # Create annotation for disease status
# annotation_col_subset <- data.frame(DiseaseStatus = disease_status_subset)
# rownames(annotation_col_subset) <- rownames(X_y_vir_z)  # Ensure rownames match

# # Plot and print the heatmap to PDF
# pheatmap(
#   heatmap_scaled_subset,
#   annotation_col = annotation_col_subset,
#   show_rownames = FALSE,
#   show_colnames = FALSE,
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   main = paste("Heatmap of Top", top_n, "Peptides by p-value - All Peptides")
# )

# 4. Table of Top Peptides
top_peptides_table_subset <- top_peptides_p_subset %>%
  select(peptide, PEP_FULL_LABEL, p.value, p.adjust, statistic, estimate1, estimate2) %>%
  rename(
    id = peptide, 
    t.statistic = statistic, 
    mean.z.control = estimate1, 
    mean.z.covidpos = estimate2 
  ) %>% # For better human readability  
  arrange(p.adjust) %>%
  mutate(PEP_FULL_LABEL = str_wrap(PEP_FULL_LABEL, width = 30)) %>%
  # Add row numbers for better readability
  mutate(row.num = row_number()) %>%
  select(row.num, everything())
  

# Define number of rows per page
rows_per_page <- 15  # Adjust based on your needs

# Split the table into chunks
table_chunks <- split(
  top_peptides_table_subset, 
  ceiling(seq_len(nrow(top_peptides_table_subset)) / rows_per_page)
)

# Create list of grobs with titles and proper spacing
grob_list <- lapply(table_chunks, function(chunk) {
  
  # Create tableGrob with padding
  table_grob <- tableGrob(chunk, rows = NULL, theme = ttheme_minimal(base_size = 8)) %>%
    gtable::gtable_add_padding(padding = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  # Create title grob
  title <- textGrob(
    "Top Peptides by Adjusted p-value - All Peptides", 
    gp = gpar(fontsize = 16, fontface = "bold")
  )
  
  # Arrange title and table vertically with relative heights
  arrangeGrob(
    grobs = list(title, table_grob),
    ncol = 1,
    heights = unit(c(1, 10), "null")  # Relative allocation to prevent overlap
  )
})

# Arrange grobs across multiple pages
pages <- marrangeGrob(
  grob_list, 
  nrow = 1, 
  ncol = 1
)

# Print the tables to PDF
print(pages)

# Close the PDF device to save the file
dev.off()

# 5. Save Intermediate Results into an Excel File for Test
wb <- createWorkbook()

# Add worksheets for each table
addWorksheet(wb, "Top_Peptides_by_Pval")
addWorksheet(wb, "T_Test_Unequal_Var")
addWorksheet(wb, "Ranked_by_AdjPval")
addWorksheet(wb, "Ranked_by_Abs_t")

# Write data to each sheet
writeData(wb, sheet = "Top_Peptides_by_Pval", top_peptides_p_subset)
writeData(wb, sheet = "T_Test_Unequal_Var", results_unequal_subset)
writeData(wb, sheet = "Ranked_by_AdjPval", results_unequal_subset %>% arrange(p.adjust))
writeData(wb, sheet = "Ranked_by_Abs_t", results_unequal_subset %>% arrange(desc(abs_t)))

# Optionally, add formatting (e.g., bold headers)
addStyle(wb, sheet = "Top_Peptides_by_Pval", 
         style = createStyle(textDecoration = "bold"), 
         rows = 1, cols = 1:ncol(top_peptides_p_subset), gridExpand = TRUE)
addStyle(wb, sheet = "T_Test_Unequal_Var", 
         style = createStyle(textDecoration = "bold"), 
         rows = 1, cols = 1:ncol(results_unequal_subset), gridExpand = TRUE)
addStyle(wb, sheet = "Ranked_by_AdjPval", 
         style = createStyle(textDecoration = "bold"), 
         rows = 1, cols = 1:ncol(results_unequal_subset), gridExpand = TRUE)
addStyle(wb, sheet = "Ranked_by_Abs_t", 
         style = createStyle(textDecoration = "bold"), 
         rows = 1, cols = 1:ncol(results_unequal_subset), gridExpand = TRUE)

# Save the workbook
saveWorkbook(wb, file = excel_file, overwrite = TRUE)