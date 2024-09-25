library(Peptides)
library(readxl)
library(readr)

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
vir_lib_subset <- vir_lib[match(rownames(vir_z), vir_lib$id), ]
test <- vir_lib_subset[1, ]
test_pep <- test$peptide

# test some functions of peptides package on a test peptide sequence 


# Q: the peptides package has a variety of features -- which are the most relevant? 
