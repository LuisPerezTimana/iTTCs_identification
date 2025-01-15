# Import libraries
library(SingleCellExperiment)
library(tidyverse)
library(Seurat)
library(velociraptor)
library(SeuratData)
library(patchwork)

# set seed
set.seed(1)

# set directory
setwd("../Data")


## Import datasets
# Import cluster 39
cluster39 <- readRDS("Cluster39/cluster39.rds")
sce_cluster39 <- as.SingleCellExperiment(cluster39)

# Rename assay
counts <- sce_cluster39@assays@data$counts
sample_id <- rownames(sce_cluster39@colData)
cell_type <- rep("cluster 39", length(sample_id))
coldata <- data.frame(sample_id, cell_type)
sce_cluster39_2 <- SingleCellExperiment(assays = list(counts = counts), colData= coldata)
sce_cluster39_2 <- scuttle::logNormCounts(sce_cluster39_2)
Seu_39 <- as.Seurat(sce_cluster39_2)
levels(Seu_39$orig.ident) <- "cluster39"

# Subset of 50 samples
Seu39 <- subset(x = Seu_39, downsample = 50)


# Import Yoshihara data

Yoshihara <- readRDS("Yoshihara data/Yoshihara_data.rds")

# Subset of 50 samples
Yoshihara_iBM <- subset(Yoshihara, downsample = 25, subset = (cell_type == "iBM_01") | (cell_type == "iBM_02"))
Yoshihara_iBM@meta.data$cell_type = "Yoshihara_iBM"
Yoshihara_primed <- subset(Yoshihara, downsample = 25, subset = (cell_type == "primed"))

# Import Yang data

Yang <- readRDS("Yanga data/Yang_data.rds")

# Rename Morula
Yang@meta.data$cell_type <- ifelse(Yang@meta.data$cell_type == "Morulae", "Morula", Yang@meta.data$cell_type)

table(Yang@meta.data$cell_type)


# Integration dataset
pbmc.combined <- merge(Yang, y = c(Seu_39, Yoshihara_iBM, Yoshihara_primed), add.cell.ids = c("Yang", "cluster39", "Yoshihara_iBM", "Yoshihara_primed"), project = "Integration")

# Tranformar into object assay5
pbmc.combined[["originalexp"]] <- as(object = pbmc.combined[["originalexp"]], Class = "Assay5")

# Normalize data
pbmc.combined <- NormalizeData(pbmc.combined)

# Find Variable Features
pbmc.combined <- FindVariableFeatures(pbmc.combined)

# Scale data
all.genes <- rownames(pbmc.combined)
pbmc.combined <- ScaleData(pbmc.combined, features = all.genes)

# PCA
pbmc.combined <- RunPCA(pbmc.combined, features = VariableFeatures(object = pbmc.combined))

# UMAP
pbmc.combined <- RunUMAP(pbmc.combined, dims = 1:10)


## Create UMAP plot by cell_type

colors <- c("#AA2225", "#F272AE", "#0000AA", "#CCD355", "#2B313F", "#14C789", "#F67452","#2A6698", "#A281B6", "#0FFEF1")

plot <- DimPlot(pbmc.combined, reduction = "umap", group.by = "cell_type", cols = colors)
plot
ggsave("../Results/UMAP_general.pdf", plot, height = 10, width = 15, units = "cm")


# Find clustering
pbmc.combined <- FindNeighbors(pbmc.combined, dims = 1:20, reduction = "pca")
pbmc.combined <- FindClusters(pbmc.combined, resolution = 2, cluster.name = "unintegrated_clusters")


# New UMAP
pbmc.combined <- RunUMAP(pbmc.combined, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# New UMAP plot recluster
plot <- DimPlot(pbmc.combined, reduction = "umap.unintegrated", group.by = c("cell_type", "seurat_clusters"), cols = colors)
ggsave("../Results/UMAP_recluster_celltype.pdf", plot, height = 10, width = 25, units = "cm")


# Find Markers
Idents(pbmc.combined) <- pbmc.combined$cell_type
pbmc.markers <- FindAllMarkers(pbmc.combined, only.pos = T)

# Top 60 genes 
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 60) %>%
  ungroup() -> top60

# Scale data for top60
pbmc.combined <- ScaleData(pbmc.combined, features = unique(top60$gene))

# PCA and UMAP
pbmc.combined <- RunPCA(pbmc.combined, features = VariableFeatures(object = pbmc.combined))
pbmc.combined <- RunUMAP(pbmc.combined, dims = 1:30, reduction = "pca", reduction.name = "umap.top60")

# Create UMAP Top60
plot <- DimPlot(pbmc.combined, reduction = "umap.top60", group.by = c("cell_type", "seurat_clusters"), cols = colors)
plot
ggsave("../Results/UMAP_recluster_celltype_top60.pdf", plot, height = 10, width = 25, units = "cm")


# Perform Integration
pbmc.combined[["originalexp"]] <- split(pbmc.combined[["originalexp"]], f = pbmc.combined$orig.ident)

# Integrate Layers
pbmc.combined <- IntegrateLayers(
  object = pbmc.combined, method = CCAIntegration, assay = "originalexp",
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = F, k.weight = 50
)

# Join Layers
pbmc.combined[["originalexp"]]  <- JoinLayers(pbmc.combined[["originalexp"]] )


# Run UMAP
pbmc.combined <- RunUMAP(pbmc.combined, dims = 1:30, reduction = "integrated.cca")
plot <- DimPlot(pbmc.combined, reduction = "umap", group.by = c("cell_type"), cols = colors)
plot
ggsave("../Results/UMAP_integrate_orig_cell_type.pdf", plot, height = 10, width = 25, units = "cm")

# RUN UMAP split by orig
plot <- DimPlot(pbmc.combined, reduction = "umap", split.by = "orig.ident", cols = colors)
plot
ggsave("Results/UMAP_integrate_split.pdf", plot, height = 10, width = 25, units = "cm")

# Create UMAP and Violin plots by specific genes
genes <- c("DPPA3", "TPRX1", "DUX4", "CCNA1", "KHDC1L", "ALPG", "H3F3A","H3F3B","H3.Y", "MBD3L2", "TRIM43", "TC2N", "SLC34A2", "ALPP", "TRIM64B", "MYCL", "MBD3L2B", "TRIM60", "JPH3", "TMEM254-AS1", "AL662789.1", "PRAMEF12", "SYT2", "LEUTX", "DUXA", "DUXB", "PARP1")
for (i in genes){
  plot <- FeaturePlot(pbmc.combined, i, reduction = "umap")
  ggsave(paste0("../Results/UMAP_",i,".pdf"), plot, height = 10, width = 15, units = "cm")
  plot <- VlnPlot(pbmc.combined, features = i, group.by = "cell_type", pt.size = 0, cols = colors)
  ggsave(paste0("../Results/Vln_",i,".pdf"), plot, height = 10, width = 15, units = "cm")

}

# Create Dot plot by specific genes 
plot <- DotPlot(pbmc.combined, features = genes, group.by = "cell_type", cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()
ggsave(paste0("../Results/DotPlot_genes.pdf"), plot, height = 15, width = 30, units = "cm")
plot


# Top 20 8CLC genes

# Find Markers
Idents(pbmc.combined) <- pbmc.combined$cell_type
pbmc.markers <- FindAllMarkers(pbmc.combined, only.pos = T)
pbmc.markers %>% filter(cluster == "8-cell") %>% slice_head(n = 20)-> top20_8clc


# Create UMAP and Violin plots by top20 
for (i in c(top20_8clc$gene,"LEUTX", "SLC34A2")){
  plot <- FeaturePlot(pbmc.combined, i, reduction = "umap")
  ggsave(paste0("../Results/Top20/UMAP_",i,".pdf"), plot, height = 10, width = 15, units = "cm")
  plot <- VlnPlot(pbmc.combined, features = i, group.by = "cell_type", pt.size = 0, cols = colors)
  ggsave(paste0("../Results/Top20/Vln_",i,".pdf"), plot, height = 10, width = 15, units = "cm")
  
}

# Create Dot plot by top20 genes
plot <- DotPlot(pbmc.combined, 
                features =  c(top20_8clc$gene,"LEUTX", "SLC34A2"), 
                group.by = "cell_type", cols = c("blue", "red"), 
                dot.scale = 8) +
  RotatedAxis()

ggsave(paste0("..Results/Top20/DotPlot_genes.pdf"), plot, height = 12, width = 25, units = "cm")
plot




