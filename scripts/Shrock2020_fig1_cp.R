################## TODO LOGS ###################################################
# 24-11-18-TODO -- create separate script for various regression e.g. logistic regression, random forest, LASSO/GROUP LASSO, RIDGE 

# 24-11-18 
# TODO NO RANDOM SLICING FIRST THRESHOLD AND TAKE COUNT 
# AND THEN MERGE WITH PAT METADATA
# for the novel coronavirus just randomly select a row from the covlib hits (obviously virscan library doesn't contain novel coronavirus info)

# 24-11-19-TODO: debug get_counts; use get_per_rep_hits for now -- row-based method 


################################################################################



# set working directory to where script is
SCRIPT_PATH <- rstudioapi::getSourceEditorContext()$path
setwd(dirname(SCRIPT_PATH))

# load packages 
library(Peptides)
library(readxl)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tidyverse)
library(tibble)
library(data.table)
library(styler)
library(ggplot2)
source("helper_funcs.R")










DF_NAMES <- c("pat_meta", "cov_z", "vir_z", "cov_lib", "vir_lib")
FORM_PATH <- "../data/processed/Shrock_2020_formatted"
FORMAT <- "rds"
formatted_df_ls <- load_filtered_data(directory = FORM_PATH, format = FORMAT, keywords = DF_NAMES)
cov_lib <- formatted_df_ls$cov_lib.rds
cov_z <- formatted_df_ls$cov_z.rds
pat_meta <- formatted_df_ls$pat_meta.rds
vir_lib <- formatted_df_ls$vir_lib.rds
vir_z <- formatted_df_ls$vir_z.rds



# PRO_DAT_PATH <- dirname(FORM_PATH) # data/processed
# VIZ_PATH <- "../results/eda/viz"
# MERGED_NAME <- "pat_meta_vir_z_2024-10-22.rds"
# # load merged patient metadata--virscan z score dataframe 
# pat_vir_z <- readRDS(file.path(PRO_DAT_PATH, MERGED_NAME)) # is a giant file btw will take a while to load

# lists of strings to query for 
# S_HUM_COR: severe human coronaviruses
# BAT_COR: bat coronaviruses 
# C_HUM_COR: common human coronaviruses
NEW_COR <- c("SARS-Cov-2", "Severe acute respiratory syndrome coronavirus 2")
S_HUM_COR <- c("SARS-CoV", "MERS-CoV", "Middle East respiratory")
BAT_COR <- c("BtCoV-279", "BtCoV-Rp3", "BtCoV-HKU3")
C_HUM_COR <- c("HCoV-HKU1", "HCoV-229E", "HCoV-NL63", "HCoV-OC43")
ALL_KEYWORDS <- c(NEW_COR, S_HUM_COR, BAT_COR, C_HUM_COR)
# ALL_KEYWORDS_SUBSTR <- unlist(str_split(ALL_KEYWORDS,"-"))
# virlib hits
v_cor_hits <- vir_lib %>% 
  filter(if_any(everything(), ~ str_detect(., str_c("(?i)\\b(", str_c(ALL_KEYWORDS, collapse = "|"), ")\\b"))))
# covlib hits 
c_cor_hits <- cov_lib %>% 
  filter(if_any(everything(), ~ str_detect(., str_c("(?i)\\b(", str_c(ALL_KEYWORDS, collapse = "|"), ")\\b")))) # lol organism is misspelled - "Organsim"

# CHECKPOINT
## quick intermission to fix typo in cov lib and overwrite it in the formatted data directory
# cov_lib <- rename(cov_lib, Organism = Organsim)
# COV_LIB_SAVE_PATH <- FORM_PATH
# COV_LIB_NAME <- "cov_lib"
# FORMAT <- "rds"
# save_as(cov_lib, COV_LIB_SAVE_PATH, COV_LIB_NAME, FORMAT)

# # first need to merge vir_z_t and pat metadata before filtering positive and negative -- existing "pat_vir_z" is faulty -- has way too many rows
# # assign patient number to pat_meta (each patient has two replicates (samples) so each row corresponds to a sample not a patient)
pat_meta$patient <- rownames(pat_meta)
# 
# # Create a version of pat_meta with rep_2 and rep_1 stacked into a single column
pat_meta_long <- pat_meta %>%
  pivot_longer(cols = c(rep_1, rep_2), names_to = "original_id_column", values_to = "rep_id") %>%
  select(rep_id, everything())  # remove the helper column

pat_meta_cols <- colnames(pat_meta_long)

# Now perform the left join with tibble_1
vir_z_pat <- vir_z_t %>%
  inner_join(pat_meta_long, by = "rep_id") %>%
  select(!!!pat_meta_cols, everything())

# CHECKPOINT
# VIR_Z_PAT_SAVE_PATH <- "/n/home01/egraff/hdsi-vector-ai/data/processed"
# VIR_Z_NAME <- "vir_z_pat"
# FORMAT <- "rds"
# save_as(vir_z_pat, VIR_Z_SAVE_PATH, VIR_Z_NAME, FORMAT)

# some smaller subsets of the vir_z_pat tibble -- use for testing new features quickly
sm_vir_z_pat <- vir_z_pat %>% select(1:50) 
lg_vir_z_pat <- vir_z_pat %>% select(1:10000)    

# filter z-scores for randomly sampled covid-related peptides of interest from covid+ and pre-covid subgroups 
cov_pat <- vir_z_pat %>% filter(COVID19_status == "positive")
pre_ctr <- vir_z_pat %>% filter(COVID19_status != "positive")

# # extract peptide ids to match covid+/covid- patients with 
# v_shc_sample_ids <- v_shc_sample$id
# c_shc_sample_ids <- c_shc_sample$id
# v_bcor_sample_ids <- v_bcor_sample$id
# c_bcor_sample_ids <- c_bcor_sample$id
# v_chc_sample_ids <- v_chc_sample$id
# c_chc_sample_ids <- c_shc_sample$id
# c_ncor_sample_ids <- c_ncor_sample$id

# preserve pat meta columns when filtering 
PAT_META_COLS <- c("rep_id", "Sample", "Library", "IP_reagent", "COVID19_status", "Hospitalized", "Pull_down_antibody", "patient", "original_id_column")



### Want to transpose vir_z dataframe so rows correspond to patients and columns correspond to peptides 
not_pep_cols <- colnames(vir_z)[!sapply(colnames(vir_z), function(x) grepl("_", x))]
# make a copy of vir_z with only peptide columns i.e. remove group, id columns
vir_z_cp <- as.data.frame(vir_z[, setdiff(colnames(vir_z), not_pep_cols)])
# these column names i.e. rep ids all start with an 'X' -- probably artifact from loading "raw" data -- remove here
formatted_colnames <- unlist(lapply(colnames(vir_z_cp), function(x) substr(x, 2, nchar(x))))


# transpose vir_z dataframe 
vir_z_t <- as_tibble(t(vir_z_cp))

vir_z_t_rownames <- formatted_colnames
vir_z_t_colnames <- c("rep_id", rownames(vir_z_cp))

vir_z_t$rep_id <- formatted_colnames
# colnames(vir_z_t) <- vir_z_t_colnames
vir_z_t <- vir_z_t %>%
  rename_with(
    ~ vir_z_t_colnames[which(vir_z_t_colnames != "rep_id")], 
    .cols = -rep_id
  )

names(vir_z_t) <- make.unique(names(vir_z_t))  # Ensure unique names
names(vir_z_t) <- replace_na(names(vir_z_t), "unnamed_col")  # Replace NA names

vir_z_t <- vir_z_t %>%
  select(rep_id, everything()) %>% # reorder so rep_id comes first
  mutate(rep_id = as.character(rep_id)) # convert to character -- pat_meta has rep_id encoded as character
# CHECKPOINT
# VIR_Z_SAVE_PATH <- "/n/home01/egraff/hdsi-vector-ai/data/processed"
# VIR_Z_NAME <- "vir_z_t"
# FORMAT <- "rds"
# save_as(vir_z_t, VIR_Z_SAVE_PATH, VIR_Z_NAME, FORMAT)

# filter covid-variant-specific peptide information from the Virscan library based on keywords
NEW_COR <- c("SARS-Cov-2", "Severe acute respiratory syndrome coronavirus 2")
S_HUM_COR <- c("SARS-CoV", "MERS-CoV", "Middle East")
BAT_COR <- c("BtCoV-279", "BtCoV-Rp3", "BtCoV-HKU3", "Bat coronavirus 279", "Bat coronavirus Rp3", "Bat coronavirus HKU3")
C_HUM_COR <- c("HCoV-HKU1", "HCoV-229E", "HCoV-NL63", "HCoV-OC43")
ALL_KEYWORDS <- c(NEW_COR, S_HUM_COR, BAT_COR, C_HUM_COR)
# ALL_KEYWORDS_SUBSTR <- unlist(str_split(ALL_KEYWORDS,"-"))
# virlib hits, shc = severe human coronavirus, bcor = bat coronavirus, chc = common human coronavirus, ncor = novel coronavirus
# Apply the function to each pattern list
v_shc_hits <- filter_hits(vir_lib, S_HUM_COR)
v_bcor_hits <- filter_hits(vir_lib, BAT_COR) # no hits
v_chc_hits <- filter_hits(vir_lib, C_HUM_COR)
v_ncor_hits <- filter_hits(vir_lib, NEW_COR) # no hits
# covlib hits
c_shc_hits <- filter_hits(cov_lib, S_HUM_COR)
c_bcor_hits <- filter_hits(cov_lib, BAT_COR) 
c_chc_hits <- filter_hits(cov_lib, C_HUM_COR)
c_ncor_hits <- filter_hits(cov_lib, NEW_COR)

# 24-11-12 Update: Do not need to randomly sample for the purposes of replicating 1D...
# they counted all hits with a z score larger than threshold

Z_THRESH <- 3.5
per_rep_shc_hits <- get_per_rep_hit(vir_z_pat, v_shc_hits, Z_THRESH)
per_rep_bcor_hits <- get_per_rep_hit(vir_z_pat, v_bcor_hits, Z_THRESH)
per_rep_chc_hits <- get_per_rep_hit(vir_z_pat, v_chc_hits, Z_THRESH)
per_rep_ncor_hits <- get_per_rep_hit(vir_z_pat, v_ncor_hits, Z_THRESH)

# per each replicate, need to sum across binary columns to get total number of hits per class of coronavirus peptides n 
per_rep_shc_count <- per_rep_shc_hits %>% mutate(shc_count = rowSums(across(-PAT_META_COLS))) %>% select(PAT_META_COLS, shc_count)
per_rep_bcor_count <- per_rep_bcor_hits %>% mutate(bcor_count = rowSums(across(-PAT_META_COLS))) %>% select(PAT_META_COLS, bcor_count)
per_rep_chc_count <- per_rep_chc_hits %>% mutate(chc_count = rowSums(across(-PAT_META_COLS))) %>% select(PAT_META_COLS, chc_count)
per_rep_ncor_count <- per_rep_ncor_hits %>% mutate(ncor_count = rowSums(across(-PAT_META_COLS))) %>% select(PAT_META_COLS, ncor_count)

count_dfs <- list(per_rep_shc_count, per_rep_bcor_count, per_rep_chc_count, per_rep_ncor_count)
counts_pos <- reduce(count_dfs, ~ full_join(.x, .y, by = PAT_META_COLS))
counts_pos <- arrange(counts_pos, COVID19_status)

ct_cols = c("shc_count", "bcor_count","chc_count", "ncor_count")

#CHECKPOINT 
# PER_REP_SAVE_PATH <- "/n/home01/egraff/hdsi-vector-ai/data/processed/per_rep_hits"
# COUNT_NAMES <- c("per_rep_shc_count", "per_rep_bcor_count", "per_rep_chc_count", "per_rep_ncor_count")
# FORMAT <- "rds"
# lapply(COUNT_NAMES, function(name) {
#   object <- get(name)
#   save_as(object = object, save_path = PER_REP_SAVE_PATH, name = name, format = FORMAT)
# })
# save_as(merged_count_dfs, PER_REP_SAVE_PATH, "merged_counts", "rds")


library(ggplot2)
library(tidyr)

counts_pos <- merged_count_dfs %>% filter(COVID19_status == "positive")
counts_neg <- merged_count_dfs %>% filter(COVID19_status == "negative")

# Reshape the data into long format
heatmap_data_pos <- counts_pos %>%
  select(rep_id, shc_count, bcor_count, chc_count, ncor_count) %>%
  pivot_longer(cols = c(shc_count, bcor_count, chc_count, ncor_count),
               names_to = "Count_Type", 
               values_to = "Count")
# Reshape the data into long format
heatmap_data_neg <- counts_neg %>%
  select(rep_id,shc_count, bcor_count, chc_count, ncor_count) %>%
  pivot_longer(cols = c(shc_count, bcor_count, chc_count, ncor_count),
               names_to = "Count_Type", 
               values_to = "Count")

# Create the heatmaps
ggplot(heatmap_data_pos, aes(x = rep_id , y = Count_Type, fill = Count)) +
  geom_tile(color = "white") +  # Add gridlines
  scale_fill_gradient(low = "white", high = "blue") +  # Customize color gradient
  labs(title = "Covid + patients peptide hits", x = "replicates", y = "virus class", fill = "number of hits") +
  theme_minimal() +
  theme(axis.text.x = element_blank())  # Rotate x-axis labels

ggplot(heatmap_data_neg, aes(x = rep_id , y = Count_Type, fill = Count)) +
  geom_tile(color = "white") +  # Add gridlines
  scale_fill_gradient(low = "white", high = "blue") +  # Customize color gradient
  labs(title = "Pre-Covid controls peptide hits", x = "replicates", y = "virus class", fill = "number of hits") +
  theme_minimal() +
  theme(axis.text.x = element_blank())  # Rotate x-axis labels



### TEST HEATMAP FOR A SMALL NUMBER OF CORONAVIRUSES ##### 

library(tidyr)
library(dplyr)
library(ggplot2)

# Assume your tibble is called cov_pat_v_shc
# Transpose the tibble
org_names <- id_to_org(vir_lib, v_shc_sample$id)
cov_v_shc_z <- cov_pat_v_shc %>% select(c("patient",  all_of(as.character(v_shc_sample$id)))) %>%
  rename_with(~ org_names[.], .cols=all_of(as.character(v_shc_sample$id)))

# Transpose the tibble and keep only the two relevant rows
cov_v_shc_z_transposed <- cov_v_shc_z %>%
  select(patient, `Human SARS coronavirus (SARS-CoV) (Severe acute respiratory syndrome coronavirus)`, 
         `Middle East respiratory syndrome coronavirus`) %>%
  pivot_longer(-patient, names_to = "Virus", values_to = "Value") %>%
  pivot_wider(names_from = patient, values_from = Value)

# Convert to a longer format for plotting
cov_v_shc_z_long <- cov_v_shc_z_transposed %>%
  pivot_longer(cols = -Virus, names_to = "Patient", values_to = "Value") %>%
  mutate(Value = as.numeric(Value))  # Ensure Value is numeric

# Convert Patient to a factor to ensure the correct order on the x-axis
cov_v_shc_z_long <- cov_v_shc_z_long %>%
  mutate(Patient = factor(Patient, levels = unique(cov_v_shc_z$patient)))

# Plot the heatmap
ggplot(cov_v_shc_z_long, aes(x = Patient, y = Virus, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +  # Darker colors for higher values
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Heatmap of Virus Response by Patient",
       x = "Patient Number", y = "Virus")



### DUPLICATE WORKING COPY##### 
library(tidyr)
library(dplyr)
library(ggplot2)

rm(cov_v_shc_z, cov_v_shc_z_transposed, cov_v_shc_z_long) # clear cached variables relevant to this cell

# v_shc_sample_ids <- v_shc_sample$id
# c_shc_sample_ids <- c_shc_sample$id
# v_bcor_sample_ids <- v_bcor_sample$id
# c_bcor_sample_ids <- c_bcor_sample$id
# v_chc_sample_ids <- v_chc_sample$id
# c_chc_sample_ids <- c_shc_sample$id
# c_ncor_sample_ids <- c_ncor_sample$id

v_sample_ids <- c(v_shc_sample$id, v_bcor_sample$id, v_chc_sample$id)
c_sample_ids <- c(c_shc_sample$id, c_bcor_sample$id, c_shc_sample$id, c_ncor_sample$id)


# # Assume your tibble is called cov_pat_v_shc
# # Transpose the tibble
# v_org_names <- id_to_org(vir_lib, v_sample_ids)  # Use the get_organisms function to map IDs to organism names
# c_org_names <- id_to_org(cov_lib, c_sample_ids)
# c_org_names <- c_org_names[1:length(c_org_names)]# contains two entries for "SARS covid 2 so discount the second entry"
# c_sample_ids <- c_sample_ids[c_sample_ids %in% names(c_org_names)] # 56158 does not have corresponding names in c_org_names
# 
# cov_v_z <- vir_z_pat %>% 
#   select(c(PAT_META_COLS, all_of(as.character(v_sample_ids)))) %>%
#   rename_with(~ v_org_names[.], .cols = all_of(as.character(v_sample_ids))) %>%
#   arrange(COVID19_status)
# 
# cov_c_z <- vir_z_pat %>% 
#   select(c(PAT_META_COLS,  all_of(as.character(c_sample_ids)))) %>%
#   rename_with(~ c_org_names[.], .cols = all_of(as.character(c_sample_ids))) %>%
#   arrange(COVID19_status)
# 
# # Step 1: Filter and select the relevant columns
# # Select patient ID, COVID-19 status, and the columns in v_org_names
# cov_v_z_subset <- cov_v_z %>%
#   select(patient, COVID19_status, all_of(v_org_names))
# 
# # Step 2: Transpose the subset and pivot it for plotting, handling duplicates
# # Transpose so each virus becomes a row and each patient becomes a column
# cov_v_z_long <- cov_v_z_subset %>%
#   pivot_longer(cols = -c(patient, COVID19_status), names_to = "Virus", values_to = "Value") %>%
#   pivot_wider(names_from = patient, values_from = Value, values_fn = list(Value = mean))  # Use `mean` to handle duplicates
# 
# # Convert the data to a long format for ggplot, excluding `COVID19_status`
# cov_v_z_long <- cov_v_z_long %>%
#   pivot_longer(cols = -Virus, names_to = "Patient", values_to = "Value") %>%
#   mutate(Value = as.numeric(Value))  # Ensure Value is numeric
# 
# # Step 3: Add COVID-19 status information for plotting
# # Merge to add COVID19_status as a column based on the patient ID
# cov_v_z_long <- cov_v_z_long %>%
#   left_join(cov_v_z_subset %>% select(patient, COVID19_status), by = c("Patient" = "patient"))
# 
# # Step 4: Plot the heatmap
# # Determine positions to add the red separator line
# separator_position <- cov_v_z_long %>%
#   group_by(COVID19_status) %>%
#   summarize(pos = max(as.numeric(factor(Patient)))) %>%
#   filter(COVID19_status == "negative") %>%
#   pull(pos)
# 
# ggplot(cov_v_z_long, aes(x = factor(Patient, levels = unique(Patient)), y = Virus, fill = Value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "blue", high = "red") +  # Darker colors for higher values
#   theme_minimal() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_text(size = 10),
#         axis.title.y = element_blank(),
#         axis.title.x = element_blank()) +
#   labs(title = "Heatmap of Virus Response by Patient") +
#   geom_vline(xintercept = separator_position + 0.5, color = "red", size = 1.5)  # Add red separator line

### logistic regression on just "cov_v" viruses e.g. combination of severe human coronaviruses, bat coronaviruses and common human coronaviruses 
cov_v_z_input <- cov_v_z %>%
  select(-all_of(PAT_META_COLS), COVID19_status) %>%
  mutate(COVID19_status = if_else(COVID19_status == "positive", 1, 0))

cov_v_z_glm <- glm(COVID19_status ~ ., data = cov_v_z_input, family = binomial)
cov_v_z_glm_coef <- summary(model)$coefficients
cov_v_z_glm_coef
