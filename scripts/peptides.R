library(Peptides)
library(readxl)
library(tidyverse)
library(here)

# set working directory to where script is
current_file_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_file_dir)

# load data
shrock_data_dir <- "../data/shrock_2020/"
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

# remove unncessary columns 
cov_z <- cov_z[, !colnames(cov_z) %in% c("group", "input")]
vir_z <- vir_z[, !colnames(vir_z) %in% c("group", "input")]

# subset of virscan library whose peptides were detected in samples
vir_lib_subset <- vir_lib[vir_lib$id %in% rownames(vir_z), ]

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

# load the saved rds file
vir_lib_pep_feats <- readRDS(file = here::here("./results/", "vir_lib_pep_funcs_results.rds"))
vir_lib_pep_feats <- as.data.frame(vir_lib_pep_feats)