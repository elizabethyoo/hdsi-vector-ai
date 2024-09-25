# load packages
library(SingleCellExperiment)
library(Biobase)
library(TOAST)
library(MuSiC)
library(reshape)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(Seurat)
library(readxl)
library(readr)

# NHP and Human data preprocessing code 

# goal: write code to read and preprocess raw count data into appropriate inputs
# for analysis pipielines
# currently testing it on one set of count data 2024_03_11_VirScan_counts.csv 
# and its corresponding sample metadata and peptide metadata. 
# TODO: generalize it to apply to any count data, granted they are all in the same format


setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set working directory to where script is 

read_files_to_dataframes <- function(directory) {
  library(readxl)
  # get list of all .csv and .xlsx files in the directory
  files <- list.files(path = directory, pattern = "\\.(csv|xlsx)$", full.names = TRUE)
  
  # loop through files
  for (file in files) {
    # extract the base file name without the extension
    file_name <- tools::file_path_sans_ext(basename(file))
    
    # make names valid R variable name
    file_name <- make.names(file_name)
    
    # if variable name starts with a number, prepend 'X' to make it valid R variable name
    if (grepl("^[0-9]", file_name)) {
      file_name <- paste0("X", file_name)
    }
    
    # Read the file into a dataframe
    if (grepl("\\.csv$", file)) {
      df <- read.csv(file)
    } else if (grepl("\\.xlsx$", file)) {
      df <- read_excel(file)
    }
    
    # assign variables as global variables
    assign(file_name, df, envir = .GlobalEnv)
  }
}

# load data and metadata
prim_path <- '../data/human-primate_8-25-2024/Primate_Dengue-and-YellowFever' 
read_files_to_dataframes(prim_path)
prim_anno_path <- '../data/human-primate_annotations_8-25-2024'
read_files_to_dataframes(prim_anno_path)

#count data
prim_vir_ct <- X2024_03_11_VirScan_counts
# patient metada 
prim_vir_sample_meta <- Primate.sample.annotations[Primate.sample.annotations$Library == 'VirScan', ]
# virscan 
prim_vir_pep_meta <- VirScan_2_Key

## are we supposed to merge the data? --> don't think so 

##### sanity checks 
# View the first few column names of the count matrix
head(colnames(prim_vir_ct))
# View the first few row names of the sample metadata
head(rownames(prim_vir_sample_meta))
# View the first few row names of the peptide metadata
head(rownames(prim_vir_pep_meta))
# View the first few row names of the count matrix
head(rownames(prim_vir_ct))


# Function to count replicates of PBS; in raw data PBS counts are off by 1 compared
# to other peptides 
adjust_sample_names <- function(sample_names) {
  # Initialize a character vector to store the adjusted names
  adjusted_names <- character(length(sample_names))
  
  # counter for occurrences of string 'PBS'
  pbs_counter <- 0
  
  for (i in seq_along(sample_names)) {
    if (grepl("PBS", sample_names[i])) {
      # Increment the counter each time "PBS" is found
      pbs_counter <- pbs_counter + 1
      # Append the count to the name
      adjusted_names[i] <- paste0("PBS-", pbs_counter)
    } else if (grepl("\\.", sample_names[i])) {
      # If there's a period, replace it with a dash and increment the number
      name_parts <- unlist(strsplit(sample_names[i], "\\."))
      base_name <- name_parts[1]
      replicate_number <- as.numeric(name_parts[2]) + 1
      adjusted_names[i] <- paste0(base_name, "-", replicate_number)
    } else {
      # Append '-1' if there's no period
      adjusted_names[i] <- paste0(sample_names[i], "-1")
    }
  }
  
  return(adjusted_names)
}


# apply adjust_sample_names() to count mx 
adjusted_colnames <- adjust_sample_names(colnames(prim_vir_ct)[-1])
# check that adjusted names look correct
print(data.frame(original = colnames(prim_vir_ct)[-1], adjusted = adjusted_colnames))
# update adjusted column names 
colnames(prim_vir_ct)[-1] <- adjusted_colnames
# set id to rownames to differentiate it from samples 
rownames(prim_vir_ct) <- prim_vir_ct$id
prim_vir_ct$id <- NULL

# Replace any occurrence of "PBS" followed by any characters with just "PBS"
prim_vir_sample_meta$`Sample ID` <- gsub("PBS.*", "PBS", prim_vir_sample_meta$`Sample ID`)
prim_vir_sample_meta$merge_id <- paste(prim_vir_sample_meta$`Sample ID`, prim_vir_sample_meta$`Replicate #`, sep='-')

# Verify if all column names of count matrix now match row names of sample metadata
all(colnames(prim_vir_ct) %in% prim_vir_sample_meta$merge_id)

# rename 'id' column of prim_vir_pep_meta to 'pep_id' to distinguish it from sample id
colnames(prim_vir_pep_meta)[colnames(prim_vir_pep_meta) == 'id'] <- 'pep_id'


#################### end of reading and formatting files into R dataframes #### 


# at this point, the count matrix prim_vir_count and sample metadata prim_vir_sample_meta
# should be mergable based on the concatenation of sample name and replicate number,
# which corresponds to column names of the former and the column 'merge_id' of the latter

# check that the above claim is true  
setdiff(colnames(prim_vir_ct), prim_vir_sample_meta$merge_id)

######## SEURAT stuff ##########################################################
options(ggrepel.max.overlaps = Inf)

prim_vir_s <- CreateSeuratObject(counts = prim_vir_ct, project = 'prim_vir_240311')
prim_vir_s <- NormalizeData(prim_vir_s, normalization.method = "LogNormalize", scale.factor = 10000)
prim_vir_s <- FindVariableFeatures(prim_vir_s, selection.method = "vst", nfeatures = 2000)

# top 10 peptides with the greatest sample-to-sample variation i.e., highly expressed
# in some samples and lowly expressed in some others 
top10 <- head(VariableFeatures(prim_vir_s), 10)
top10_species <- prim_vir_pep_meta[prim_vir_pep_meta$pep_id %in% top10,]$Species
plot1 <- VariableFeaturePlot(prim_vir_s)
plot2 <- LabelPoints(plot = plot1, points = top10, label = top10_species, repel = TRUE, xnudge = 0, ynudge = 0)

saveRDS(prim_vir_s, file = "../results/human-primate/prim_vir_s.rds")

### TODO: add metadata into Seurat object
### TODO: save Seurat object as singleCellExperiment Object
### TODO: reproduce UMAP plots

# # note that you can chain multiple commands together with %>%
# prim_vir_s <- SCTransform(prim_vir_s) %>%
#   RunPCA() %>%
#   FindNeighbors(dims = 1:30) %>%
#   FindClusters() %>%
#   RunUMAP(dims = 1:30)
# DimPlot(prim_vir_s, reduction = "umap")

# naming convention: primate --> pri vs. human --> hum, virscan --> vir vs. 
# vectorscan --> vec , date --> YYMMDD e.g. hum_vec_240311 

#### ##### What melodi did with SEURAT
Plate13_VirScan_Seurat <- CreateSeuratObject(counts = Plate13_counts_Z_threshold_2)

# Extract column (cell) names
column_names <- colnames(Plate13_VirScan_Seurat)


# Determine if each cell is "sample" or "PBS"
#condition_labels <- ifelse(grepl("PBS", column_names), "PBS", "Sample")



# Add condition labels to Seurat object metadata
#Plate13_VirScan_Seurat$condition <- condition_labels

Plate13_VirScan_Seurat <- SCTransform(Plate13_VirScan_Seurat)

# Perform PCA
Plate13_VirScan_Seurat <- RunPCA(Plate13_VirScan_Seurat, features = VariableFeatures(Plate13_VirScan_Seurat), npcs = min(nrow(Plate13_counts_Z_threshold_2), ncol(Plate13_counts_Z_threshold_2)-1))
