# Import the libraries
library(tidyverse)
library(Seurat)
library(velociraptor)
library(SingleCellExperiment)

### Processing Yoshihara data

# Set directory
setwd("../data/Yoshihara data/")

# Import counts H9
counts <- read_tsv("E-MTAB-10581_raw-count.txt")

# Extract genes and eliminate this column
genes <- counts$ID
counts <- counts[, -1]

# Extract col names for 3 (12h_DOX_01), 5(12h_DOX_02) and 6(12h_noDOX_01)
col_names <- c()
for (i in colnames(counts)){
  num <- unlist(str_split(i, "-"))[2]
  if (num == "3"){
    col_names <- c(col_names, i)
  } else if (num == "5"){
    col_names <- c(col_names, i)
  } else if (num == "6"){
    col_names <- c(col_names, i)
  }
  
}

# Create the matrix
Yoshihara <- counts[col_names]
rownames(Yoshihara) <- genes

# Create colData
sample_id <- col_names
cell_type <- c(rep("iBM_01", 6850),
               rep("iBM_02", 6687),
               rep("primed", 4997))
df <- data.frame(sample_id, cell_type)

# Create count matrix general
sce_Y <- SingleCellExperiment(assays = list(counts = Yoshihara), colData= df)
sce_Y <- scuttle::logNormCounts(sce_Y)

# Create Seurat object
Seu_Y <- as.Seurat(sce_Y)
levels(Seu_Y$orig.ident) <- "Yoshihara Data"

# Save seurat as rds
SaveSeuratRds(Seu_Y, file = "Yoshihara_data.rds")



### Processing YANG DATA

# set directory
setwd("../../data/Yang data/")

# import list of cells
lista <- read_tsv("list.txt")

# extract data for dir
dir <- "GSE36552_RAW/"

setwd(dir)
archivos <- list.files()

# extract each dataset in a list df
lista_df <- list()
for (i in archivos){
  
  name <- unlist(strsplit(i, "_"))[1]
  if (name %in% lista$ID){
    df <- read_delim(i, "\t", skip = 1, col_names = c("gene", "read", "len", "coverage", "rpkm"))
    df <- df[, 1:2]
    colnames(df) <- c("genes", name)
    lista_df <- c(lista_df, list(df))
  }
  
}

# Merge df's
df_total <- lista_df[[1]]
for (i in 2:length(lista_df)){
  df_total <- merge(df_total, lista_df[[i]], by = "genes", all = T)
}

# Extract genes names
genes <- df_total$genes
df_total <- df_total[, -1]
rownames(df_total) <- genes

# Replace "NA" by 0
df_total <- mutate_all(df_total, ~replace(., is.na(.), 0))

# Generate ColData
colnames(lista) <- c("sample_id", "cell_type", "sample_num")

# Create SingleCell Experiment object
sce <- SingleCellExperiment(assays = list(counts = df_total), colData= lista)
sce <- scuttle::logNormCounts(sce)

# Create seurat object
Seu_Yang <- as.Seurat(sce)
levels(Seu_Yang$orig.ident) <- "Yang Data"

# Save seurat as rds
SaveSeuratRds(Seu_Yang, file = "../Yand_data.rds")
