# scRNA-seq analysis for identification and characterization of transient totipotent-like cells during somatic cell reprogramming

This project focuses on the preprocessing and integration of single-cell RNA sequencing (scRNA-seq) data. We analyzed data from multiple sources (Table 1) to generate integrated visualizations such as UMAP and Dotplot for specific gene expressions.

## Project Structure

The project consists of two main scripts:
- **preprocessing.R**: Preprocesses and cleans the input data from different sources.
- **Integration.R**: Integrates the preprocessed data, performs dimensional reduction (UMAP), and generates visualizations of specific gene expressions.

## Data

We use the following datasets:
- **Yoshihara Data**
- **Yang Data**
- **Cluster39**

These datasets are loaded, filtered, and preprocessed before integration.

## Requirements

The following R packages are necessary to run this project:
- Seurat
- tidyverse
- SingleCellExperiment
- velociraptor
- SeuratData
- patchwork

## Usage

1. **Preprocessing**:
   - Run `preprocessing.R` to load, filter, and preprocess the scRNA-seq data from Yoshihara, Yang, and Cluster39.
   - The script performs data cleaning, normalization, and transformation.

2. **Integration and Visualization**:
   - Run `Integration.R` to integrate the preprocessed data from the different datasets. The script performs the following tasks:
     - Integrates data into a unified dataset.
     - Reduces the dimensionality of the data using UMAP.
     - Generates UMAP plots to visualize the clustering of cells.
     - Creates a Dotplot to visualize specific gene expressions across clusters.

## Output

The following outputs will be generated:
- **UMAP Plot**: A two-dimensional visualization of the integrated cell clusters.
- **Dotplot**: A visualization showing the expression levels of selected genes across different clusters.
- **Top 20 genes**: UMAP plot and Dotplot for each of the top 20 genes

## Authors and Acknowledgments

This project was conducted by .... . The data utilized in this project originates from the Yoshihara & Kere (2022), Yang (2013), and Cluster39 datasets.

## Reference

Masahito Yoshihara & Juha Kere (2022). Single-cell RNA-sequencing of transient DUX4-expressed human embryonic stem cells. BioStudies, E-MTAB-10581. Retrieved from https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10581

Yan L, Yang M, Guo H, Yang L et al. Single-cell RNA-Seq profiling of human preimplantation embryos and embryonic stem cells. Nat Struct Mol Biol 2013 Sep;20(9):1131-9. PMID: 23934149
