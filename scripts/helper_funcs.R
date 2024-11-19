library(Peptides)
library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(tidyverse)
library(tibble)
library(data.table)
library(styler)

# set working directory to where script is
SCRIPT_PATH <- rstudioapi::getSourceEditorContext()$path
setwd(dirname(SCRIPT_PATH))


# function to map peptide id to corresponding organism names
# 
# param library: dataframe that contains both organism name and id e.g. virscan library
# param id_list: list of peptide ids (numeric)
id_to_org <- function(library, id_list) {
  library %>%
    filter(id %in% id_list) %>%  # Filter for the specified IDs
    select(id, Organism) %>%     # Select the ID and organism columns
    deframe()                     # Convert to a named vector (id -> organism)
}

# function to filter a dataframe based on matches to a pattern
# 
# param data: dataframe to filter 
# patterns: list of strings to search for -- the function will look for exact matches, case-insenstivie
#
filter_hits <- function(data, patterns) {
  data %>%
    filter(if_any(everything(), ~ str_detect(., str_c("(?i)\\b(", str_c(patterns, collapse = "|"), ")\\b"))))
}

# randomly sample one peptide in each list of hits 
# data: dataframe you want to sample from
# query: list of strings you want to query e.g. "SARS", "BtCoV"
sample_hits <- function(data, queries) {
  #  filter each query and sample, then combine results
  sampled_rows <- lapply(queries, function(query) {
    data %>%
      filter(if_any(everything(), ~ str_detect(., query))) %>%
      slice_sample(n = 1)
  })
  
  # combine the sampled rows into a single dataframe
  result <- bind_rows(sampled_rows)
  
  return(result)
}


# get_counts: function to subset a z-score dataframe (where rows = peptides, columns = patient replicates) i.e. Virscan or Covscan that matches subset of the virscan or covscan library
# and get counts of z-scores that are greater than or equal to some set threshold. Preserves patient metadata in the returned dataframe 
# PAT_META_COLS <- c("rep_id", "Sample", "Library", "IP_reagent", "COVID19_status", "Hospitalized", "Pull_down_antibody", "patient", "original_id_column")
# 
# param data: a dataframe of z-scores i.e. vir_z or cov_z
# param hits: a library (i.e. Virscan Library, Covscan Library) subdataframe containing information on virus strains of interest
# param col_or_row: specify "col" if you want a count for each fixed column (peptide) or specify "row" if you want a count for each fixed row (patient replicate)
# param threshold: a number threshold for z-scores. Only entries greater than or equal to this threshold will contribute towards the count
get_counts <- function(data, hits, threshold, col_or_row = "row") {
  PAT_META_COLS <- c("rep_id", "Sample", "Library", "IP_reagent", "COVID19_status", 
                     "Hospitalized", "Pull_down_antibody", "patient", "original_id_column")
  
  # Validate col_or_row input
  if (!col_or_row %in% c("col", "row")) {
    stop("Invalid value for 'col_or_row'. Use 'col' for column-based counts or 'row' for row-based counts.")
  }
  
  # Subset data based on hits
  data <- data %>%
    select(all_of(PAT_META_COLS), matches(hits$id))  # Subset to metadata + columns matching hits
  
  if (col_or_row == "col") {
    # Column-based counts
    ct_df <- data %>%
      select(-all_of(PAT_META_COLS)) %>%  # Exclude metadata columns
      summarise(across(everything(), ~ sum(. >= threshold, na.rm = TRUE)))
  } else if (col_or_row == "row") {
    # Row-based counts
    ct_df <- data %>%
      mutate(count = rowSums(across(-all_of(PAT_META_COLS), ~ . >= threshold, na.rm = TRUE))) %>%
      select(all_of(PAT_META_COLS), count, everything())  # Ensure PAT_META_COLS come first
  }
  
  return(ct_df)
}

# Helper function to save intermediate results in one of the three formats: RDS, CSV, or JSON.
# Creates a subdirectory to place the file if it does not exist, and generates a file name based on user input and date
# 
# param object: the r object to be saved e.g., dataframe, list
# param save_path: string specifying absolute path to the directory you wish to save your file e.g. "~/some_dir/another_dir/file_in_this_dir/"
# param format: string specifying format to save the object. Currently supported: RDS, CSV, JSON, Feather (interoperable with Python)
# name: string specifying name of the file.
# 
# 
# TODO: move processed data that's currently saved under results/processed to data/processed/
save_as <- function(object, save_path, name, format) {
    # create save directory if it doesn't exist
    if (!dir.exists(save_path)) {
        dir.create(save_path, recursive = TRUE)
    }
    
    # concatenate object name and file extension to generate file name
    file_extension <- switch(format,
                             "rds" = ".rds",
                             "csv" = ".csv",
                             "json" = ".json",
                             "feather" = ".feather",
                             "parquet" = ".parquet",
                             stop("Unsupported format. Choose 'rds', 'csv', 'json', 'feather', or 'parquet'."))
    # filename <- paste0(name, "_", Sys.Date(), file_extension)
    filename <- paste0(name, file_extension)
    full_path <- file.path(save_path, filename)
    
    # check for necessary packages depending on format
    if (format == "json" && !requireNamespace("jsonlite", quietly = TRUE)) {
        install.packages("jsonlite")
        library(jsonlite)
    }
    if (format %in% c("feather", "parquet") && !requireNamespace("arrow", quietly = TRUE)) {
        install.packages("arrow")
        library(arrow)
    }
    
    # write object based on format
    if (format == "rds") {
        saveRDS(object, file = full_path)
    } else if (format == "csv") {
        write.csv(object, file = full_path, row.names = FALSE)
    } else if (format == "json") {
        write_json(object, path = full_path, pretty = TRUE)
    } else if (format == "feather") {
        write_feather(object, full_path)
    } else if (format == "parquet") {
        write_parquet(object, full_path)
    }
    
    message("Object saved as ", format, " at ", full_path)
}



#TODO: write docstring for the function
#
load_filtered_data <- function(directory, format, keywords) {

# define file extensions
file_extension <- switch(format,
                         "rds" = "\\.rds$",
                         "csv" = "\\.csv$",
                         "json" = "\\.json$",
                         "feather" = "\\.feather$",
                         "parquet" = "\\.parquet$",
                         stop("Unsupported format. Choose 'rds', 'csv', 'json', 'feather', or 'parquet'."))

# all files in the directory with specified extension
all_files <- list.files(directory, pattern = file_extension, full.names = TRUE)

# filter files that contain at least one keyword in their name
selected_files <- all_files[sapply(all_files, function(file) {
    any(sapply(keywords, grepl, file))
})]

# function to load a single file based on format
load_file <- function(file) {
    switch(format,
           "rds" = readRDS(file),
           "csv" = read.csv(file),
           "json" = jsonlite::fromJSON(file),
           "feather" = arrow::read_feather(file),
           "parquet" = arrow::read_parquet(file),
           stop("Unsupported format"))
}

# Load each selected file as a data frame and store in a list
data_list <- lapply(selected_files, load_file)
names(data_list) <- basename(selected_files)  # Name list elements by file names

return(data_list)
}


#
#
# # load data
# shrock_data_dir <- here("data", "Shrock_2020")
# pat_meta <- read_csv(file.path(shrock_data_dir,
#                                "Shrock_2020_patient-metadata.csv"))
# cov_z <- read_csv(file.path(shrock_data_dir,
#                             "Shrock_2020_Coronavirus-screen-IgG-Zscores.csv"))
# vir_z <- read_csv(file.path(shrock_data_dir,
#                             "Shrock_2020_VirScan-IgG-Zscores.csv"))
# cov_lib <- read_csv(file.path(shrock_data_dir,
#                               "Shrock_2020_Coronavirus-library.csv"))
# vir_lib <- read_csv(file.path(shrock_data_dir,
#                               "Shrock_2020_VirScan-library.csv"))
# # make column names valid R variables
# colnames(pat_meta) <- make.names(gsub(" ", "_", colnames(pat_meta)), unique = TRUE)
# colnames(cov_lib) <- make.names(gsub(" ", "_", colnames(cov_lib)), unique = TRUE)
# colnames(vir_lib) <- make.names(gsub(" ", "_", colnames(vir_lib)), unique = TRUE)
# 
# # note: for some reason covid z score ids start from 0 vs. virscan peptide ids start from 1
# colnames(cov_lib)[1] <- "id"
# colnames(vir_z)[1] <- "id"
# 
# # remove unncessary columns
# cov_z <- cov_z[, !colnames(cov_z) %in% c("group", "input")]
# vir_z <- vir_z[, !colnames(vir_z) %in% c("group", "input")]
# 
# save_as_rds(cov_lib, "cov_lib", "./data/processed")
# save_as_rds(vir_lib, "vir_lib", "./data/processed")

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


# TODO
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


