
# Set-Up & Import #####
library(pacman)
p_load(here
       , readr , readxl
       , dplyr , tidyr
       , data.table # transpose is useful here 
       , tidymodels, parsnip, broom 
)

# anchor project root directory 
here::i_am("scripts/00_shrock_data_processing.R")


SHROCK_DPATH <- here("data", "Shrock_2020")
PROCESSED_DPATH <- here("data", "processed")
SHROCK_EXCEL_FNAME <- "Shrock_2020.xlsx"
DEMO_DF_FNAME <- "demo_df.rds"
COVID_LIB_FNAME <- "cov_lib.rds"
COVID_IgG_FNAME <- "cov_IgG.rds"
COVID_WIDE_FNAME <- "cov_wide.rds"

# 25-01-27-TODO do saveRDS in a loop 
# demo_df = read_excel(path = here(SHROCK_DPATH, SHROCK_EXCEL_FNAME)
#                      , sheet = "Patient metadata")
# # saveRDS(demo_df, file = here(PROCESSED_DPATH, DEMO_DF_FNAME))

# covid_lib = read_excel(path = here(SHROCK_DPATH, SHROCK_EXCEL_FNAME)
#                        , sheet = "Coronavirus library ") %>% 
#   rename(id = `...1`)
# # saveRDS(covid_lib, file = here(PROCESSED_DPATH, COVID_LIB_FNAME))

# covid_IgG = read_excel(path = here(SHROCK_DPATH, SHROCK_EXCEL_FNAME)
#                        , sheet = "Coronavirus screen - IgG Z") %>% 
#   mutate(id = as.character(id)) 
# # saveRDS(covid_IgG, file = here(PROCESSED_DPATH, COVID_IgG_FNAME))

# covid_wide = covid_IgG %>% select(-c(group, input)) %>% 
#   data.table::transpose(make.names = "id", keep.names="SampleID") %>% 
#   as_tibble() 
# # saveRDS(covid_wide, file = here(PROCESSED_DPATH, COVID_WIDE_FNAME))

# Model Fitting #### 

## Sample Splitting #### 
# Note data includes 
# ncol(covid_IgG %>% select(-c(id, group, input)))  = 1098
# nrow(demo_df) = 422 # (844 samples, two replicates each)
# This matches Shrock paper 

#==============================================================================
# 25-01-27-bookmark
#==============================================================================
# demo_df <- readRDS(here(PROCESSED_DPATH, DEMO_DF_FNAME))
# covid_lib <- readRDS(here(PROCESSED_DPATH, COVID_LIB_FNAME))
# covid_IgG <- readRDS(here(PROCESSED_DPATH, COVID_IgG_FNAME))
# covid_wide <- readRDS(here(PROCESSED_DPATH, COVID_WIDE_FNAME))

# covid_wide <- covid_wide %>%
#   rename(rep_id = SampleID)

# # # just to compare 
# # vir_z <- readRDS(here("data", "rds", "vir_z.rds"))
# # cov_lib <- readRDS(here("data", "rds", "cov_lib.rds"))

# # merge rep_1 and rep_2 in demo_df
# demo_tmp <- demo_df %>%
#   pivot_longer(
#     cols = c(rep_1, rep_2), # Columns to pivot
#     names_to = "rep_type",  # Add a column to identify original columns (optional)
#     values_to = "rep_id"    # New column to store the values
#   ) %>%
#   select(-rep_type) %>% # Drop the 'rep_type' column if not needed 
#   select(rep_id, Sample, COVID19_status, everything())
# saveRDS(demo_tmp, file = here(PROCESSED_DPATH, "cov_pat_meta.rds"))
# # Fix typos and standardize column names to match cov_lib

# cov_lib <- covid_lib %>% #renaming covid_lib to cov_lib for consistency
#   rename(Organism = Organsim,
#          Protein = `Protein name`,
#          Sequence = `Peptide sequence`,
#          start = `Start position`,
#          end = `End position`,
#          oligo = `Nucleotide sequence`) %>%
#   # preserve original column order
#   select(id, `Library type`, Organism, everything())
# saveRDS(cov_lib, file = here(PROCESSED_DPATH, "cov_lib.rds"))


# # intersect_cols <- c("id", "Organism", "Protein", "Sequence", "Species", "start", "end", "peptide", "source")
# # -> indicates fix covid_lib, <- indicates fix cov_lib
# # todo fix organsim to organism in covid_lib, "Protein name" -> "Protein", "Peptide sequence" <- "Sequence", "Start position" -> "start", "End position" -> "end", "Nucleotide sequence" -> "oligo"? 
# # for now, we only care about "id", "Organism," "Protein"

# # # vir_z_pat columns
# # VIR_Z_META_COL <- c("rep_id", "Sample", "Library", "IP_reagent", 
# #                    "Hospitalized", "Pull_down_antibody", 
# #                    "patient", "COVID19_status")

# # TODO modify helper funcs to account for "standardized, preprocessed" datasets


# labels = rbind(demo_df %>% select(SampleID=rep_1, Sample, COVID19_status) 
#              , demo_df %>% select(SampleID=rep_2, Sample, COVID19_status) 
# ) %>% arrange(Sample)

# # Trying tidymodels split objects 
# # 25-01-27-TODO: make naming schemes consistent with virscan data
# cov_z_pat <- covid_wide %>% 
#   inner_join(demo_tmp %>% select(rep_id, Sample, Library, Hospitalized, COVID19_status), "rep_id") %>% 
#   select(rep_id, Sample, Library, COVID19_status, Hospitalized, everything()) %>%
#   mutate(COVID19_status = as.factor(COVID19_status)) 
# saveRDS(cov_z_pat, file = here(PROCESSED_DPATH, "cov_z_pat.rds"))

# load saved data
cov_z_pat <- readRDS(here(PROCESSED_DPATH, "cov_z_pat.rds"))
cov_lib <- readRDS(here(PROCESSED_DPATH, "cov_lib.rds"))


COV_META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "Hospitalized") # no COVID19_status
COV_META_COL <- c(COV_META_COL_NO_STATUS, "COVID19_status") # includes COVID19_status

X_y_cov_z <- cov_z_pat %>%
  select(-all_of(META_COL_NO_STATUS)) %>%
  mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
  # put COVID19_status as the first column
  select(COVID19_status, everything()) %>% 
  mutate(COVID19_status = as.factor(COVID19_status)) # convert to factor

cat("Step: Data Preparation\n")
cat("Dimensions of X_y_cov_z:", dim(X_y_cov_z), "\n")
print(head(X_y_cov_z))  # Inspect a small portion
X_cov_z <- X_y_cov_z %>%
    select(-COVID19_status)

y_cov_z <- X_y_cov_z$COVID19_status

# Confirm dimensions
cat("Peptide Data Dimensions:", dim(X_cov_z), "\n")
cat("COVID19_status Dimensions:", length(y_cov_z), "\n")

# check if covid_lib has any NAs -- none found
# cat("Number of NAs in covid_lib:", sum(is.na(cov_lib)), "\n")

# modify get_organism_from_id to account for covid_lib
# TODO change get_organism_from_id to get_vir_organism_by_id
# Need to extract: Protein, start, end, (peptide) sequence
# also cannot query on the headnode for rhinovirus A, Human Herpesvirus4 and HIV-1
get_cov_organism_by_id <- function(id_chr, cov_lib) {
  id_numeric <- as.numeric(id_chr)
  
  # Filter the row corresponding to the given id
  row_data <- cov_lib %>% 
    filter(id == id_numeric)
  
  # Check if row_data is empty
  if (nrow(row_data) == 0) {
    return(NA)  # Return NA if no matching id is found
  }
  # Try to get the Organism, Protein, Sequence, start, and end
  label <- row_data %>%
    pull(Organism) %>%
    first()
  
  protein <- row_data %>%
    pull(Protein) %>%
    first()
  
  sequence <- row_data %>%
    pull(Sequence) %>%
    first()
  
  start_pos <- row_data %>%
    pull(start) %>%
    first()
  
  end_pos <- row_data %>%
    pull(end) %>%
    first()
  
  # Combine the extracted information into a single string
  label <- paste("Organism:", label, 
                 "Protein:", protein, 
                 "Sequence:", sequence, 
                 "Start:", start_pos, 
                 "End:", end_pos)
  
  return(label)
}

test_cov_ids <- colnames(cov_z_pat)[(7:13)] # first 6 columns are metadata
test_org_names <- sapply(test_cov_ids, get_cov_organism_by_id, cov_lib = cov_lib)
# 25-01-27 TODO figure out meaningful "organism" resolution in covid_library

