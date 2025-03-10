if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
glmnet,
caret,
dplyr,
ggplot2,
reshape2,
here,
gt,
stringr,
data.table,
purrr
)

# anchor project root directory 
here::i_am("scripts/all_peptides_lasso_scratch.R")

# helper function -- works on any of covscan, virscan, or combination of both so long as
# peptide id is in "c_" or "v_" format
# see formatted data in /n/home01/egraff/hdsi-vector-ai/data/processed/Shrock_vir_cov_combined
get_organism_by_id <- function(id_chr, cov_lib, vir_lib) {
  # Ensure the id columns are character for comparison if they are factors
  if (is.factor(cov_lib$id)) cov_lib$id <- as.character(cov_lib$id)
  if (is.factor(vir_lib$id)) vir_lib$id <- as.character(vir_lib$id)
  
  if (startsWith(id_chr, "c_")) {
    # If cov_lib$id is numeric, extract a numeric value; otherwise, use the full id string.
    if (is.numeric(cov_lib$id)) {
      id_val <- as.numeric(sub("^c_", "", id_chr))
    } else {
      id_val <- id_chr
    }
    row_data <- cov_lib %>% filter(id == id_val)
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
    if (is.numeric(vir_lib$id)) {
      id_val <- as.numeric(sub("^v_", "", id_chr))
    } else {
      id_val <- id_chr
    }
    row_data <- vir_lib %>% filter(id == id_val)
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

# Writing to file ===========================================================
#TEST ON A SUBSET
# Define number of peptides for testing - default empty string
NUM_TEST_PEPTIDES <- ""

# CHECK BEFORE RUNNING - CHANGE TO APPROPRIATE NAMES ####################################################################
LASSO_RESULTS_FNAME <- "lasso-result"
DATASET_NAME <- paste0("combined", NUM_TEST_PEPTIDES)
BASE_FNAME <- paste0(LASSO_RESULTS_FNAME, "_", DATASET_NAME)
EXCEL_FNAME <- paste0(BASE_FNAME, "_", DATASET_NAME, ".xlsx")

# Create directory for this run
RUN_DIR <- here::here("results", "eda", "reglasso", paste0(BASE_FNAME, "_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))
if (!dir.exists(RUN_DIR)) {
    dir.create(RUN_DIR, recursive = TRUE)
}
# CHECK BEFORE RUNNING - CHANGE TO APPROPRIATE NAMES ####################################################################


# load data 
# pull formatted dfs from /n/home01/egraff/hdsi-vector-ai/data/processed/Shrock_vir_cov_combined
DATA_DIR <- here::here("data", "processed", "Shrock_vir_cov_combined")
cat("Step: Loading and preparing data...\n")
com_z_pat <- readRDS(here::here(DATA_DIR, "combined_z_pat.rds")) # combined
vir_z_pat <- readRDS(here::here(DATA_DIR, "vir_z_pat_renamed_ids.rds")) # virscan
cov_z_pat <- readRDS(here::here(DATA_DIR, "cov_z_pat_renamed_ids.rds")) # covscan

cov_lib   <- readRDS(here::here(DATA_DIR, "cov_lib_renamed_ids.rds"))  # if needed
vir_lib   <- readRDS(here::here(DATA_DIR, "vir_lib_renamed_ids.rds"))  # if needed

META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "Hospitalized")

# CHECK BEFORE RUNNING - CHANGE TO APPROPRIATE NAMES ####################################################################
# TODO: separate out data processing steps (pertaining to specific datasets) into a separate script
# ideally, this script should not refer to dataset-specific variables
X_y_com_z <- com_z_pat %>%
  select(-all_of(META_COL_NO_STATUS)) %>%
  mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
  select(COVID19_status, everything())
  
# TEST ON A SUBSET ==========================================================
if (NUM_TEST_PEPTIDES != "") {
  set.seed(12345)  
  TEST_COL_FNAME <- paste0(BASE_FNAME, "_", "test_columns.rds")
  test_peptide_columns <- sample(colnames(X_y_com_z[,-1]), as.numeric(NUM_TEST_PEPTIDES))
  # save in RUN_DIR
  saveRDS(test_peptide_columns, file.path(RUN_DIR, TEST_COL_FNAME))
  X_y_com_z <- X_y_com_z %>% select(all_of(c("COVID19_status", test_peptide_columns)))
}

X_com_z <- X_y_com_z %>% select(-COVID19_status)
y_com_z <- X_y_com_z$COVID19_status

cat("Peptide Data Dimensions:", dim(X_com_z), "\n")
cat("Response Vector Length:", length(y_com_z), "\n")

X <- X_com_z %>%
  as.matrix()

y <- y_com_z %>%
  as.numeric()
# CHECK BEFORE RUNNING - CHANGE TO APPROPRIATE NAMES ####################################################################


# fit model ###############################################################################################################
lasso <- cv.glmnet(x=X, y=y, alpha=1,family="binomial",type.measure = "mse")
LASSO_FNAME <- paste0(BASE_FNAME, "_", "model_fit.rds")
saveRDS(lasso, here::here(RUN_DIR, LASSO_FNAME))

# TEMP - DELETE AFTER TEMP TASK
lasso <- readRDS("/n/home01/egraff/hdsi-vector-ai/results/eda/reglasso/lasso-result_combined_2025-03-10_21-32-43/lasso-result_combined_model_fit.rds")

plot(lasso)
title(main="LASSO covid status vs. all z-scores")

l_min <- lasso$lambda.min
coef <- coef(lasso, s="lambda.min")
nonzero_coef <- coef[coef[,1]!=0,][-1] # also exclude intercept
print(nonzero_coef)


## Relating LASSO coefficients to peptides

nonzero_coef_df <- data.frame( # extract peptide ids
  coef = nonzero_coef, # extract nonzero coefficients
  stringsAsFactors = FALSE
)
COEF_FNAME_RDS <- paste0(BASE_FNAME, "_", "coefs.rds")
saveRDS(nonzero_coef_df, here::here(RUN_DIR, COEF_FNAME_RDS))

TOP_NUM = 50 # some number of top coefficients we're interested in

# Function to get the first N words of a string
get_first_n_words <- function(text, n = 10) {
  words <- unlist(strsplit(text, "\\s+"))  # Split the text into words
  paste(head(words, n), collapse = " ")  # Combine the first N words
}

nonzero_coef_df <- nonzero_coef_df %>% 
  arrange(desc(coef)) %>%
  slice_head(n = TOP_NUM) %>%
  mutate(id = rownames(.)) %>%
  mutate(
    attributes = map(id, get_organism_by_id, cov_lib, vir_lib),
    organism = map_chr(attributes, "Organism"),
    Protein = map_chr(attributes, "Protein"),
    Sequence = map_chr(attributes, "Sequence"),
    Start = map_dbl(attributes, "Start"),
    End = map_dbl(attributes, "End")
  ) %>%
  select(-attributes) %>% # Remove the intermediate list column
  filter(!is.na(organism)) %>%
  distinct(organism, .keep_all = TRUE) %>% 
  mutate(
    first_n_words = sapply(organism, get_first_n_words, n = 10),  # Get the first N words of `organism`
    first_char_protein = substr(Protein, 1, 1),                  # Get the first character of `Protein`
    start_end = paste0(Start, "-", End),                        # Concatenate Start and End with "-"
    PEP_FULL_LABEL = paste(first_n_words, first_char_protein, start_end, sep = ", ")  # Combine all components
  ) %>%
  select(-first_n_words, -first_char_protein, -start_end)

COEF_FNAME_CSV <- paste0(BASE_FNAME, "_", "coefs.csv")
write.csv(nonzero_coef_df, COEF_FNAME_CSV)
# %>%

# format into nice looking table 
nonzero_coef_tb <- nonzero_coef_df %>%
  select(id, PEP_FULL_LABEL, coef) %>% 
  gt() %>%
  tab_header(
    title = paste0("Top ", TOP_NUM, " coefficients from LASSO on ", DATASET_NAME, " dataset(s)")
  )
print(nonzero_coef_tb)
# webshot::install_phantomjs() # for saving to .png but doesn't work -- could be a FAS server thing
TB_FNAME <- paste0(BASE_FNAME,"_table",".html")
gtsave(data = nonzero_coef_tb, filename=TB_FNAME, path = RUN_DIR)

