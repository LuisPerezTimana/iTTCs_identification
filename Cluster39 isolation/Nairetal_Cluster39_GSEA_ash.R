

#packages 
install.packages('Seurat')
library(Seurat)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DropletUtils")

library(DropletUtils)
library(Matrix)
library(irlba)
library(ggplot2) # Tidyverse is pre-installed, yay!
library(dplyr)
library(scico)
library(ggpointdensity)
library(ggrepel)
library(tidyverse)
library(scales)
library(patchwork)
library(tidyselect)
library(cowplot)
install.packages('devtools')
devtools::install_github('immunogenomics/presto')

library("presto")

BiocManager::install("fgsea")
library(fgsea)

install.packages("msigdbr")
library(msigdbr)
install.packages("SCpubr") #better viz umap plots 
library(SCpubr)
library(viridis)
library(assertthat)
library(ggplotify)

# reading the seruat object by Nair et al depsoited here 
#https://kundajelab.github.io/reprogramming-browser/home.html


Nair_etal_allcluster <- readRDS("/Users/alishariati/Desktop/iTTManuscript/scRNA/seurat.rds")
#update for seurat compatiblity 
Nair_etal_allcluster <- UpdateSeuratObject(object = Nair_etal_allcluster) 
NairetAlUmap <- DimPlot(Nair_etal_allcluster , reduction = "umap") 
NairetAlUmapSCpubr <- SCpubr::do_DimPlot(sample = Nair_etal_allcluster)
Markers <- c("DUXA", "TPRX1", "CCNA1")
NairetAlUmapFeature <- SCpubr::do_FeaturePlot(sample = Nair_etal_allcluster , features = Markers) + 
  scale_colour_gradient(low = "white", high = "magenta", na.value = NA)

#gene enrcihment analysis 
Nair_etal.genes <- wilcoxauc(Nair_etal_allcluster, 'cluster')
# we have all the genes for each cluster
dplyr::count(Nair_etal.genes, group)


#making a gene list from Yan et al and this paper as hallmark of 8C, epiblast and maternal 

#https://www.cell.com/developmental-cell/fulltext/S1534-5807(24)00574-4#mmc3

List1 <- read.table("/Users/alishariati/Desktop/iTTManuscript/tzzhangjuan1-LINE1_hESC-1fcc1c8/RNA-seq analysis/GeneListYan_2023/Maternal_Yan2023.gmt.csv", header = F, stringsAsFactors = F)
List1 <- as.vector(List1$V1)
List2 <- read.table("/Users/alishariati/Desktop/iTTManuscript/tzzhangjuan1-LINE1_hESC-1fcc1c8/RNA-seq analysis/GeneListYan_2023/8cell_genes_Yan.gmt.csv", header = F, stringsAsFactors = F)
List2 <- as.vector(List2$V1)
List3 <- read.table("/Users/alishariati/Desktop/iTTManuscript/tzzhangjuan1-LINE1_hESC-1fcc1c8/RNA-seq analysis/GeneListYan_2023/Epiblast_genes_Yan.gmt.csv", header = F, stringsAsFactors = F)
List3 <- as.vector(List3$V1)
List4 <- read.table("/Users/alishariati/Desktop/iTTManuscript/tzzhangjuan1-LINE1_hESC-1fcc1c8/RNA-seq analysis/GeneListYan_2023/8cell_Morula_genes_Yan.gmt.csv", header = F, stringsAsFactors = F)
List4 <- as.vector(List4$V1)
List5 <- read.table("/Users/alishariati/Desktop/iTTManuscript/tzzhangjuan1-LINE1_hESC-1fcc1c8/RNA-seq analysis/GeneListYan_2023/Primed_hESCs_genes_Yan.gmt.csv", header = F, stringsAsFactors = F)
List5 <- as.vector(List5$V1)
List6 <- read.table("/Users/alishariati/Desktop/iTTManuscript/tzzhangjuan1-LINE1_hESC-1fcc1c8/RNA-seq analysis/GeneListYan_2023/Zygotic_genes_Yan.gmt.csv", header = F, stringsAsFactors = F)
List6 <- as.vector(List6$V1)
## TO CONTINUE the lists (gene sets)....

pathways <- list("Maternal" = List1, "Zygotic"= List6, "8cell" = List2, "8cell_morula"= List4, "epiblast" = List3, "Primed_hESC"=List5)  # Name each gene list, and list them together as pathway
str(head(pathways))

#cluster 39 genes to be used for ranking 

cluster39.markers <- FindMarkers(Nair_etal_allcluster, ident.1 = 39)
head(cluster39.markers, n=10)

#rownames to  a column 

cluster39.markers <- tibble::rownames_to_column(cluster39.markers, "features")

#selcting fold changes and gene name 
cluster39.ranks <-  cluster39.markers %>%  dplyr::select(features, avg_log2FC)

cluster39.ranks <- deframe(cluster39.ranks)

head(cluster39.ranks)


#gene set enrichment 

fgseaRes <- fgsea(pathways, cluster39.ranks, minSize=15, maxSize=3000)
head(fgseaRes) 

#plotting Normalize Enrichment Value 
#Note that zygotic and epiblast genes don't appear in the plot becuase of padj filter (high pvalue)

ggplot(fgseaRes %>% filter(padj < 0.008) %>% head(n= 50), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 0)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

ggplot(fgseaRes, aes(x=pathway, y=padj)) +  geom_boxplot(fill='green')

#heatmap of all markers for all clusters

all.markers <- FindAllMarkers(Nair_etal_allcluster, only.pos = TRUE)
all.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

all.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(Nair_etal_updated, features = top10$gene) + NoLegend()
