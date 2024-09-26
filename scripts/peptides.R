library(Peptides)
library(readxl)
library(tidyverse)

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
dim(which(is.na(vir_lib)))


# subset of virscan library whose peptides were detected in samples
vir_lib_subset <- vir_lib[vir_lib$id %in% vir_z$id, ]

# list of peptides functions that take AA sequences as arguments
# and whose outputs might be useful for downstream analysis
# see https://github.com/dosorio/Peptides/tree/master for descriptions
pep_funcs <- c("aaComp", "aIndex")
# Apply each function to the 'peptide' column of vir_lib_subset and save the outputs into a new dataframe
results <- sapply(pep_funcs, function(func) {
    do.call(func, list(vir_lib_subset$peptide))
})

# Convert the results to a dataframe
results_df <- as.data.frame(results)

# Set the column names to match the function names
colnames(results_df) <- pep_funcs

# Print the resulting dataframe
print(results_df)

# ... sowhat  