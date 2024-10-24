library(Peptides)
library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(data.table)

# set working directory to where script is
here::i_am("scripts/peptides.R")
library(here)

# Helper function to save intermediate results
# dir: specify directory location relative to root e.g. "scripts/processed/"
# TODO: specify file format e.g. RDS, dataframe, csv, JSON, etc. 
# TODO: move processed data that's currently saved under results/processed to data/processed/
save_as_rds <- function(object, name, dir) {
  # Create subdirectory if it doesn't exist
  processed_dir <- paste("./", dir, sep = "")
  # processed_dir <- here::here("./results/processed/")
  if (!dir.exists(processed_dir)) {
    dir.create(processed_dir, recursive = TRUE)
  }
  # Generate the filename and save the object
  filename <- paste0(name, "_", Sys.Date(), ".rds")
  saveRDS(object, file = file.path(processed_dir, filename))
}


# load data
shrock_data_dir <- here("data", "Shrock_2020")
pat_meta <- read_csv(file.path(shrock_data_dir,
                               "Shrock_2020_patient-metadata.csv"))
cov_z <- read_csv(file.path(shrock_data_dir,
                            "Shrock_2020_Coronavirus-screen-IgG-Zscores.csv"))
vir_z <- read_csv(file.path(shrock_data_dir,
                            "Shrock_2020_VirScan-IgG-Zscores.csv"))
cov_lib <- read_csv(file.path(shrock_data_dir,
                              "Shrock_2020_Coronavirus-library.csv"))
vir_lib <- read_csv(file.path(shrock_data_dir,
                              "Shrock_2020_VirScan-library.csv"))
# make column names valid R variables
colnames(pat_meta) <- make.names(gsub(" ", "_", colnames(pat_meta)), unique = TRUE)
colnames(cov_lib) <- make.names(gsub(" ", "_", colnames(cov_lib)), unique = TRUE)
colnames(vir_lib) <- make.names(gsub(" ", "_", colnames(vir_lib)), unique = TRUE)

# note: for some reason covid z score ids start from 0 vs. virscan peptide ids start from 1
colnames(cov_lib)[1] <- "id"
colnames(vir_z)[1] <- "id"

# remove unncessary columns
cov_z <- cov_z[, !colnames(cov_z) %in% c("group", "input")]
vir_z <- vir_z[, !colnames(vir_z) %in% c("group", "input")]

save_as_rds(cov_lib, "cov_lib", "./data/processed")
save_as_rds(vir_lib, "vir_lib", "./data/processed")

# # list of peptides functions that take AA sequences as arguments
# # and whose outputs might be useful for downstream analysis
# # see https://github.com/dosorio/Peptides/tree/master for descriptions
# # exclude functions that require more information than just the peptide sequence
# pep_funcs <- c("aaComp", "aIndex", "boman", "charge", "crucianiProperties",
#                "fasgaiVectors", "hmoment", "hydrophobicity", "instaIndex",
#                "kideraFactors", "lengthpep", "massShift", "membpos",
#                "mswhimScores", "mw", "mz", "pI", "protFP", "stScales",
#                "tScales", "vhseScales", "zScales")

# TODO: save datframes with corrected names



# # sanity check on one peptide sequence
# test_pep <- vir_lib_subset$peptide[1]
# test_results <- sapply(pep_funcs, function(func) {
#     do.call(func, list(test_pep))
# })
# saveRDS(test_results, file = here::here("./results", "test_pep_funcs_results.rds"))
# test_results_tmp <- readRDS(file = here::here("./results", "test_pep_funcs_results.rds"))

# Apply each function to the 'peptide' column of vir_lib_subset and save the outputs into a new dataframe
# pep_funcs_results <- sapply(pep_funcs, function(func) {
#     do.call(func, list(vir_lib$peptide))
# })
# # Save resulting list as a rds file
# saveRDS(pep_funcs_results, file = here::here("./results/", "vir_lib_pep_funcs_results.rds"))
# system("say saved entire virscan library pep func results") # sound alert when done 

# # load the saved rds file
# vir_lib_pep_feats <- readRDS(file = here::here("./results/", "vir_lib_pep_funcs_results.rds"))
# vir_lib_pep_feats <- as.data.frame(vir_lib_pep_feats)
# vir_lib_pep_feats <- cbind(vir_lib, vir_lib_pep_feats)
# # remove duplicate id rows
# vir_lib_pep_feats_u <- vir_lib_pep_feats[unique(vir_lib_pep_feats$id), ]
# 
# vir_z_means <- as.data.frame(rowMeans(vir_z[-1]))
# 
# colnames(vir_z_means) <- "sample_means"
# vir_z_means$id <- vir_z$id
# 
# # remove duplicate ids
# vir_dup_ids <- vir_lib[duplicated(vir_lib$id), ]
# vir_lib_unique_ids <- vir_lib[unique(vir_lib$id), ]
# # Q: what id do I use to merge covid data?
# cov_lib_unique_ids <- cov_lib[unique(cov_lib$id), ]
# 
# # sanity check for duplicate ids 
# dups <- vir_lib_unique_ids[duplicated(vir_lib_unique_ids$id),]
# lib_subset <- vir_lib[c("id", "oligo", "peptide")] %>% distinct()
# vir_lib_pep_feats_sub <- vir_lib_pep_feats[c(
#     "id", "aaComp", "aIndex", "boman", "charge", "crucianiProperties",
#     "fasgaiVectors", "hmoment", "hydrophobicity", "instaIndex",
#     "kideraFactors", "lengthpep", "massShift", "membpos",
#     "mswhimScores", "mw", "mz", "pI", "protFP", "stScales",
#     "tScales", "vhseScales", "zScales"
# )] %>% distinct()
# 
# # make a temporary merged dataframe of mean z scores and peptide library info, peptide features...
# vir_merged <- merge(vir_z_means, lib_subset, by = "id", x.all = TRUE)
# vir_merged <- merge(vir_merged, vir_lib_pep_feats_sub, by = "id", x.all = TRUE)
# saveRDS(vir_merged, file = here::here("./results/", "vir_z_pep_funcs_results.rds"))
# 
# # some covariates are scalars and others are vectors -- how to represent vectors? 
# pep_feats_sc <- c("aIndex", "boman", "charge",  "hmoment",
#                   "hydrophobicity", "instaIndex", "lengthpep", "mw", "mz", "pI")
# pep_feats_vec <- setdiff(pep_funcs, pep_feats_sc)
# vir_merged <- vir_merged %>% mutate_at(all_of(pep_feats_sc), as.numeric)
# 
# lm_str <- paste("sample_means~", paste(pep_feats_sc, collapse="+"), sep = "")
# lm_z_means_feats_sc <- lm(lm_str, data = vir_merged)
# summary(lm_z_means_feats_sc)
# save_as_rds(lm_z_means_feats_sc, "lm_z_means_feats_sc", "./results/")
# 
# ## 2024-10-22-tues 
# ## trying to merge patient metadata with virscan z scores, so we can have patient
# ## sample-level covid and hospitalized status witj their corresponding z-score 
# ## for various peptides
# pat_meta <- pat_meta %>% select(rep_1, rep_2, COVID19_status, Hospitalized)
# 
# pat_meta <- pat_meta %>%
#   # Create the 'sample' column as row number
#   mutate(sample = row_number()) %>%
#   # Reshape the data, turning 'rep_1' and 'rep_2' into a single column 'rep'
#   pivot_longer(cols = c(rep_1, rep_2), 
#                names_to = "rep", 
#                values_to = "rep_value") 
#   # Optionally remove the 'replicate' column if you no longer need it
#   # select(-replicate)
# 
# # View the result
# print(pat_meta)
# select(vir_z, -any_of(2))
# # swap rows and columns of vir_z
# vir_z_long <- vir_z %>%
#   pivot_longer(cols = matches("S\\d+$"),   # Regex to match columns containing 'S' followed by digits and an underscore
#                names_to = "rep_value",     # Store column names (rep_values) in 'rep_value'
#                values_to = "z_score")
# 
# 
# # this line made R run out of memory... need to find more efficient ways to merge
# # trying data.table package which apparently makes working with large 
# # datasets memory-efficient
# # tmp <- pat_meta %>% inner_join(vir_z_long, by="rep_value")
# 
# # convert data.frames to data.tables, merge data.tables, and convert back to data.frames
# pat_meta_dt <- as.data.table(pat_meta)
# vir_z_long_dt <- as.data.table(vir_z_long)
# 
# merged_data_dt <- vir_z_long_dt[pat_meta_dt, on = "rep_value", allow.cartesian = TRUE]
# save_as_rds(merged_data_dt, "pat_meta_vir_z", "./data/")

#TODO rename merged_data_dt to something more descriptive


