library(dplyr)
library(ggplot2)
library(cowplot)
library(randomForest)
library(reshape2)
library(here)
library(gt)
library(Matrix)
library(stringr)
library(data.table)
library(purrr)
library(broom)
library(parallel) 
library(gridExtra)  # For arranging multiple plots/tables
library(ggrepel)    # For enhanced plot labeling
library(openxlsx)
library(grid)
library(pbapply)    # For parallel processing with progress bar
library(RColorBrewer)
library(patchwork)

# anchor project root directory 
here::i_am("scripts/04_eda_randomforest.R")

#============================================================================

# load data
vir_z_pat <- readRDS(here::here("data", "processed", "vir_z_pat.rds"))
vir_lib <- readRDS(here::here("data", "rds", "vir_lib.rds"))
# vir_z <- readRDS(here::here("data", "rds", "vir_z.rds"))

#============================================================================

### data prep
META_COL_NO_STATUS <- c("rep_id", "Sample", "Library", "IP_reagent", 
                   "Hospitalized", "Pull_down_antibody", 
                   "patient", "original_id_column") # excludes COVID19_status which will be a part of the final dataframe


# X_y_vir_z includes z-scores and COVID19_status vs. X_vir_z includes z-scores only
X_y_vir_z <- vir_z_pat %>%
    select(-all_of(META_COL_NO_STATUS)) %>%
    mutate(COVID19_status = ifelse(is.na(COVID19_status), 0, as.numeric(COVID19_status == "positive"))) %>%
    # put COVID19_status as the first column
    select(COVID19_status, everything()) %>% #factorize covid status
    mutate(COVID19_status = as.factor(COVID19_status))
X_y_vir_z_dt <- as.data.table(X_y_vir_z)
# drop any rows 

X_vir_z <- X_y_vir_z %>%
    select(-COVID19_status)
X_vir_z_dt <- as.data.table(X_vir_z)


# Confirm dimensions
cat("Peptide Data Dimensions:", dim(X_vir_z_dt), "\n")
