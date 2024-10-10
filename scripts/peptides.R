library(Peptides)
library(readxl)
library(dplyr)
library(here)

# set working directory to where script is
current_file_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_file_dir)
# set root directory
here::set_here("/Users/ecyoo/github/elizabethyoo/llm-malaria")

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

# note: for some reason covid z score ids start from 0 vs. virscan peptide ids start from 1
colnames(cov_lib)[1] <- "id"
# colnames(vir_z)[1] <- "id"

# remove unncessary columns 
cov_z <- cov_z[, !colnames(cov_z) %in% c("group", "input")]
vir_z <- vir_z[, !colnames(vir_z) %in% c("group", "input")]

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
# pep_funcs_results <- sapply(pep_funcs, function(func) {
#     do.call(func, list(vir_lib$peptide))
# })
# # Save resulting list as a rds file
# saveRDS(pep_funcs_results, file = here::here("./results/", "vir_lib_pep_funcs_results.rds"))
# system("say saved entire virscan library pep func results") # sound alert when done 

# load the saved rds file
vir_lib_pep_feats <- readRDS(file = here::here("./results/", "vir_lib_pep_funcs_results.rds"))
vir_lib_pep_feats <- as.data.frame(vir_lib_pep_feats)
vir_lib_pep_feats <- cbind(vir_lib, vir_lib_pep_feats)
# remove duplicate id rows
vir_lib_pep_feats_u <- vir_lib_pep_feats[unique(vir_lib_pep_feats$id), ]

vir_z_means <- as.data.frame(rowMeans(vir_z[-1]))
colnames(vir_z_means) <- "sample_means"
vir_z_means$id <- vir_z$id

# remove duplicate ids
vir_dup_ids <- vir_lib[duplicated(vir_lib$id), ]
vir_lib_unique_ids <- vir_lib[unique(vir_lib$id), ]
# Q: what id do I use to merge covid data?
cov_lib_unique_ids <- cov_lib[unique(cov_lib$id), ]

# sanity check for duplicate ids 
dups <- vir_lib_unique_ids[duplicated(vir_lib_unique_ids$id),]
lib_subset <- vir_lib[c("id", "oligo", "peptide")] %>% distinct()
vir_lib_pep_feats_sub <- vir_lib_pep_feats[c(
    "id", "aaComp", "aIndex", "boman", "charge", "crucianiProperties",
    "fasgaiVectors", "hmoment", "hydrophobicity", "instaIndex",
    "kideraFactors", "lengthpep", "massShift", "membpos",
    "mswhimScores", "mw", "mz", "pI", "protFP", "stScales",
    "tScales", "vhseScales", "zScales"
)] %>% distinct()


# make a temporary merged dataframe of mean z scores and peptide library info, peptide features...
vir_merged <- merge(vir_z_means, lib_subset, by = "id", x.all = TRUE)
vir_merged_tmp <- merge(vir_merged, vir_lib_pep_feats_sub, by = "id", x.all = TRUE)
sum(duplicated(vir_merged$id))


