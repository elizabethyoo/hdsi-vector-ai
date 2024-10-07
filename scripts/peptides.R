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
cor_z <- read_csv(file.path(shrock_data_dir,
                            "Shrock_2020_Coronavirus-screen-IgG-Zscores.csv"))
vir_z <- read_csv(file.path(shrock_data_dir,
                            "Shrock_2020_VirScan-IgG-Zscores.csv"))
cor_lib <- read_csv(file.path(shrock_data_dir,
                              "Shrock_2020_Coronavirus-library.csv"))
vir_lib <- read_csv(file.path(shrock_data_dir,
                              "Shrock_2020_VirScan-library.csv"))


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
    do.call(func, list(vir_lib_subset$peptide))
})
# Save resulting list as a rds file
saveRDS(pep_funcs_results, file = here::here("./results/", "vir_lib_subset_pep_funcs_results.rds"))
# system("say saved pep func results") # sound alert when done 

# load the saved rds file
pep_funcs_results <- readRDS(file = here::here("./results/", "vir_lib_subset_pep_funcs_results.rds"))
rownames(pep_funcs_results) <- vir_lib_subset$id

# Calculate means of all rows except for values in 'group' and 'id' columns
vir_z$means <- rowMeans(vir_z[, -c(1, 2)], na.rm = TRUE)
vir_z_aIndex <- as.matrix(cbind(c(vir_z_means), c(pep_funcs_results[,"aIndex"])))
ggplot()

# Identify the row indices of the vir_z dataframe
all_indices <- 1:nrow(vir_z)

# Identify the row indices of the 'mean' column (non-NA values)
mean_indices <- which(is.na(vir_z$means))

# Find the missing row index
missing_index <- setdiff(all_indices, mean_indices)

# Print the missing row index
print(paste("The 'mean' column is missing a row at index:", missing_index))