
# Set-Up & Import #####
library(pacman)
p_load(here
       , readr , readxl
       , dplyr , tidyr
       , data.table # transpose is useful here 
       , tidymodels, parsnip, broom 
)

# anchor project root directory 
here::i_am("scripts/00_1_shrock_data_combining.R")


#==============================================================================
# 25-02-11: combine covidscan and virscan (1) z-scores and (2) libraries
# convention: we want to follow rows - peptides, columns - samples/replicates/patients 
#==============================================================================
WRITE_DIR <- here("data", "processed", "Shrock_vir_cov_combined")
# if directory doesn't exist, create it
if (!dir.exists(WRITE_DIR)) {
  dir.create(WRITE_DIR)
}

# virscan data============================================================================
# load data
vir_z_pat <- readRDS(here::here("data", "processed", "vir_z_pat.rds"))
vir_lib <- readRDS(here::here("data", "rds", "vir_lib.rds"))
vir_z <- readRDS(here::here("data", "rds", "vir_z.rds"))

### data prep
VIR_META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "IP_reagent", 
                   "Hospitalized", "Pull_down_antibody", 
                   "patient", "original_id_column") # excludes COVID19_status which will be a part of the final dataframe
VIR_META_COL <- c(VIR_META_COL_NO_STATUS, "COVID19_status") # includes COVID19_status

# # X_y_vir_z includes z-scores and COVID19_status vs. X_vir_z includes z-scores only
# X_y_vir_z <- vir_z_pat %>%
#     select(-all_of(META_COL_NO_STATUS)) %>%
#     mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
#     # put COVID19_status as the first column
#     select(COVID19_status, everything())


# covscan data============================================================================
# load data
# covscan data
cov_z_pat <- readRDS(here::here("data", "processed", "cov_z_pat.rds"))
cov_lib <- readRDS(here::here("data", "processed", "cov_lib.rds"))
cov_z <- readRDS(here::here("data", "processed", "cov_IgG.rds"))

# Prepare data
COV_META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "Hospitalized")
COV_META_COL <- c(COV_META_COL_NO_STATUS, "COVID19_status")

# X_y_cov_z <- cov_z_pat %>%
#   select(-all_of(META_COL_NO_STATUS)) %>%
#   mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
#   select(COVID19_status, everything()) %>% 
#   mutate(COVID19_status = as.factor(COVID19_status))

# find common rows in "Sequence" column of cov_lib and vir_lib
common_peptides <- intersect(cov_lib$Sequence, vir_lib$Sequence)
common_lib_cols <- intersect(colnames(cov_lib), colnames(vir_lib))

# concatenate string "c" with all id column values in cov_lib dplyr chain
cov_lib <- cov_lib %>% 
  mutate(id = paste("c", id, sep = "_"))
saveRDS(cov_lib, file = here(WRITE_DIR, "cov_lib_renamed_ids.rds"))
# concatenate string "v" with all id column values in vir_lib dplyr chain
vir_lib <- vir_lib %>%
  mutate(id = paste("v", id, sep = "_"))
saveRDS(vir_lib, file = here(WRITE_DIR, "vir_lib_renamed_ids.rds"))

# Rename all column names except COV_META_COL_NO_STATUS in cov_z_pat by concatenating string "c" with the column name
# for example, column name "1" will become "c_1"
cov_z_pat <- cov_z_pat %>%
  rename_with(~paste("c", .x, sep = "_"), -all_of(c(COV_META_COL_NO_STATUS, "COVID19_status"))) 
saveRDS(cov_z_pat, file = here(WRITE_DIR, "cov_z_pat_renamed_ids.rds"))
vir_z_pat <- vir_z_pat %>%
  rename_with(~paste("v", .x, sep = "_"), -all_of(c(VIR_META_COL_NO_STATUS, "COVID19_status")))
saveRDS(vir_z_pat, file = here(WRITE_DIR, "vir_z_pat_renamed_ids.rds"))

# TODO: Remember to save modified data
# select the first ten rows of the dataframes
vir_z_pat <- vir_z_pat %>% head(10)
cov_z_pat <- cov_z_pat %>% head(10)
vir_lib <- vir_lib %>% head(10)
cov_lib <- cov_lib %>% head(10)
dim(vir_lib)
dim(cov_lib)
dim(vir_z_pat)
dim(cov_z_pat)

# full join the two based on "id column"
combined_lib <- full_join(vir_lib, cov_lib, by = "id")
saveRDS(combined_lib, file = here(WRITE_DIR, "combined_lib.rds"))
print("saved combined_lib")

# full join the two based on "rep_id" column -- the joined dataframe will have more columns
shared_meta_cols <- intersect(COV_META_COL, VIR_META_COL)
non_shared <- setdiff(VIR_META_COL, shared_meta_cols)

combined_z_pat <- cov_z_pat %>% 
  left_join(vir_z_pat %>% select(-all_of(non_shared)), by = shared_meta_cols) 
print(head(combined_z_pat))
saveRDS(combined_z_pat, file = here(WRITE_DIR, "combined_z_pat.rds"))

print("saved combined_z_pat")

# num_vir_pep <- ncol(vir_z_pat) - length(VIR_META_COL)
# num_cov_pep <- ncol(cov_z_pat) - length(COV_META_COL)

# dim(vir_z_pat)
# dim(cov_z_pat)
# dim(combined_z_pat)

# load combined datasets
com_lib <- readRDS(here(WRITE_DIR, "combined_lib.rds"))
com_z_pat <- readRDS(here(WRITE_DIR, "combined_z_pat.rds"))

library(dplyr)

get_peptide_attributes_by_id <- function(id_chr, cov_lib, vir_lib) {
  # If the id comes from covidscan data (prefix "c_")
  if (startsWith(id_chr, "c_")) {
    # Remove the "c_" prefix and convert to numeric (assuming cov_lib$id is numeric)
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
      Start    = row_data %>% pull(start) %>% first(),
      End      = row_data %>% pull(end) %>% first()
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
      # For virscan, we pull from "Protein names" but return it as "Protein"
      Protein  = row_data %>% pull(`Protein names`) %>% first(),
      Sequence = row_data %>% pull(Sequence) %>% first(),
      Start    = row_data %>% pull(start) %>% first(),
      End      = row_data %>% pull(end) %>% first()
    ))
    
  } else {
    stop("Unrecognized peptide id prefix. It must start with either 'c_' for covidscan or 'v_' for virscan.")
  }
}

# test works
# att_test <- c("c_1", "c_2", "c_3", "v_1", "v_2", "v_3") %>% map(get_peptide_attributes_by_id, cov_lib, vir_lib)



# mutate(
#     attributes = map(peptide, get_cov_organism_by_id, cov_lib),
#     organism = map_chr(attributes, "Organism"),
#     Protein = map_chr(attributes, "Protein"),
#     Sequence = map_chr(attributes, "Sequence"),
#     Start = map_dbl(attributes, "Start"),
#     End = map_dbl(attributes, "End")
#   ) %>%