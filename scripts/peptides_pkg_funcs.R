
library(Peptides)
library(dplyr)
library(readr)
library(tidyr)
library(tidyverse)
library(tibble)
library(data.table)
library(here)

# configs -- because here package is unable to find project root directory for some reason
# adjust to your directory settings as necessary 
LIZ_DATA_DIR <- "/n/home01/egraff/hdsi-vector-ai/data/rds/"
LIZ_RESULTS_DIR <- "/n/home01/egraff/hdsi-vector-ai/results/"
COV_LIB_FNAME <- "cov_lib.rds"
VIR_LIB_FNAME <- "vir_lib.rds"
VIR_LIB_PFEATS_FNAME <- "vir_lib_pep_funcs_results.rds"

# # load data
cov_lib <- readRDS(paste0(LIZ_DATA_DIR, COV_LIB_FNAME))
vir_lib <- readRDS(paste0(LIZ_DATA_DIR, VIR_LIB_FNAME))

vir_lib_pep_funcs <- readRDS(paste0(LIZ_RESULTS_DIR, VIR_LIB_PFEATS_FNAME))

# list of peptides functions that take AA sequences as arguments
# and whose outputs might be useful for downstream analysis
# see https://github.com/dosorio/Peptides/tree/master for descriptions
# exclude functions that require more information than just the peptide sequence
pep_funcs <- c("aaComp", "aIndex", "boman", "charge", "crucianiProperties",
               "fasgaiVectors", "hmoment", "hydrophobicity", "instaIndex",
               "kideraFactors", "lengthpep", "massShift", "membpos",
               "mswhimScores", "mw", "mz", "pI", "protFP", "stScales",
               "tScales", "vhseScales", "zScales")


# # sanity check on one peptide sequence
# test_pep <- vir_lib_subset$peptide[1]
# test_results <- sapply(pep_funcs, function(func) {
#     do.call(func, list(test_pep))
# })
# saveRDS(test_results, file = here::here("./results", "test_pep_funcs_results.rds"))
# test_results_tmp <- readRDS(file = here::here("./results", "test_pep_funcs_results.rds"))

# Apply each function to the 'peptide' column of vir_lib_subset and save the outputs into a new dataframe
pep_funcs_results <- sapply(pep_funcs, function(func) {
    do.call(func, list(vir_lib$peptide))
})
# Save resulting list as a rds file
saveRDS(pep_funcs_results, file = here::here("./results/", "vir_lib_pep_funcs_results.rds"))
system("say saved entire virscan library pep func results") # sound alert when done

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