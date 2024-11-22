library(Peptides)
library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(data.table)
library(glmnet)
library(gt)

# set working directory to where script is
SCRIPT_PATH <- rstudioapi::getSourceEditorContext()$path
setwd(dirname(SCRIPT_PATH))
source("helper_funcs.R")

# set working directory to where script is
here::i_am("scripts/all_peptides_lasso_scratch.R")
library(here)
source("../scripts/all_peptides_lasso_scratch.R")

# configs -- because here package is unable to find project root directory for some reason
# adjust to your directory settings as necessary 
LIZ_DATA_DIR <- "/n/home01/egraff/hdsi-vector-ai/data/processed/"
PAT_META_FNAME <- "pat_meta_vir_z_2024-10-22.rds"

LIZ_RESULTS_DIR <- "/n/home01/egraff/hdsi-vector-ai/results/"
COV_LIB_FNAME <- "cov_lib_2024-10-24.rds"
VIR_LIB_FNAME <- "vir_lib_2024-10-24.rds"
VIR_LIB_PFEATS_FNAME <- "vir_lib_pep_funcs_results.rds"


pat_meta_vir_z <- readRDS(paste0(LIZ_DATA_DIR,PAT_META_FNAME))
# Convert to a tibble if needed
pat_meta_vir_z <- as_tibble(pat_meta_vir_z)

# Pivot to wide format using pivot_wider
pivot_pat_meta_vir_z <- pat_meta_vir_z %>%
  pivot_wider(names_from = id, 
              values_from = z_score)

pivot_pat_meta_vir_z <- pivot_pat_meta_vir_z %>% mutate(COVID19_status = if_else(COVID19_status == "positive", 1, 0),
                                                        Hospitalized = if_else(Hospitalized == "hospitalized", 1, 0))


# # View the resulting pivoted data
# head(pivot_pat_meta_vir_z)
# 
# ## save pat_meta_vir_z; ready to use as input to model
# save_as_rds(pivot_pat_meta_vir_z, "pat_meta_vir_z_input", "./data/")

# # response:, covariates: all z scores
NUM_PEPTIDES = 115753
covstat_zscores_str <- paste("COVID19_status~", paste(1:NUM_PEPTIDES, collapse="+"), sep = "")
# lm_z_means_feats_sc <- lm(lm_str, data = vir_merged)
# summary(lm_z_means_feats_sc)
# save_as_rds(lm_z_means_feats_sc, "lm_z_means_feats_sc", "./results/")

# lr_covstat_zscores <- glm()
# log_model <- glm()

## Try LR just because, try  LASSO, Elastic Net
## Otherwise, want to try dim reduction first 


set.seed(1031)
## Fit LASSO
## find penalty \lambda by cross-validation by minimizing MSE
z_scores <- as.matrix(pivot_pat_meta_vir_z[, -c(1:5)]) # just the z-score columns
cov_stat <- pivot_pat_meta_vir_z$COVID19_status
lasso_cov_z <- cv.glmnet(x=z_scores, y=cov_stat, alpha=1,family="binomial",type.measure = "mse")
# LASSO_COV_Z_FNAME <- "all_peptides_lasso_fit"
# save_as(object = lasso_cov_z , save_path = LIZ_RESULTS_DIR, name = LASSO_COV_Z_FNAME, format = "rds")

plot(lasso_cov_z)
title(main="LASSO covid status vs. all z-scores")

l_min <- lasso_cov_z$lambda.min
coef <- coef(lasso_cov_z, s="lambda.min")
nonzero_coef <- coef[coef[,1]!=0,][-1] # also exclude intercept
print(nonzero_coef)

## Relating LASSO coefficients to peptides
vir_lib <- readRDS(paste0(LIZ_DATA_DIR, VIR_LIB_FNAME))
cov_lib <- readRDS(paste0(LIZ_DATA_DIR, COV_LIB_FNAME))
nonzero_coef_df <- data.frame(
  id = names(nonzero_coef), # extract peptide ids
  coef = nonzero_coef, # extract nonzero coefficients
  stringsAsFactors = FALSE
)
# make sure id values in nonzero_coef_df and vir_lib match (want both to be numeric)
nonzero_coef_df <- nonzero_coef_df %>% 
  mutate(id = as.numeric(id))

# map back to peptides by descending order of LASSO coefficients
tmp_desc_coef <- vir_lib %>% 
  inner_join(nonzero_coef_df, by = "id") %>%
  select(coef, everything()) %>%
  arrange(desc(coef))


save_as(object = tmp_desc_coef, save_path = LIZ_RESULTS_DIR ,name = "all_virscan_peps_lasso_coeffs", format = "rds")

# print to a tidy looking table 
tmp_desc_coef <- rename(tmp_desc_coef, peptide_id = `...1`)
tmp_desc_coef_tb <- tmp_desc_coef %>% slice(1:20) %>%
  select(coef, id, Organism) %>% 
  gt() %>%
  tab_header(
    title = "Top 20 Coefficients from LASSO on All Peptides"
  )

# webshot::install_phantomjs() # for saving to .png but doesn't work -- could be a FAS server thing
# gtsave(data = tmp_desc_coef_tb, filename = "all_virscan_peps_lasso_coeffs_tb.html", path = LIZ_RESULTS_DIR)

print(tmp_desc_coef)
