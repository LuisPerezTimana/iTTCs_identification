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

This project was conducted by ESV. Integration analysis was performed by LAPT and guided by ESV. ESV and AS wrote the manuscript. AS supervised the study

Funding Sources: This work was supported by the NIGMS/NIH through a Pathway to Independence Award K99GM126027/ R00GM126027 and Maximizing Investigator Award (R35GM147395), a start-up package from the University of California, Santa Cruz (S.A.S). E.S.V is supported by the National Institutes of Health (NIH) under award number K12GM139185 and the Institute for the Biology of Stem Cells (IBSC) at UC Santa Cruz.. The content is solely the responsibility of the authors and does not necessarily represent the official views of the NIH or the IBSC. Technical support from Benjamin Abrams, UCSC Life Sciences Microscopy Center, RRID: SCR_021135.

## Reference

Masahito Yoshihara & Juha Kere (2022). Single-cell RNA-sequencing of transient DUX4-expressed human embryonic stem cells. BioStudies, E-MTAB-10581. Retrieved from https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10581

Nair, S., Ameen, M., Sundaram, L., Pampari, A., Schreiber, J., Balsubramani, A., Wang, Y. X., Burns, D., Blau, H. M., Karakikes, I., Wang, K. C., & Kundaje, A. (2023). Transcription factor stoichiometry, motif affinity and syntax regulate single-cell chromatin dynamics during fibroblast reprogramming to pluripotency. bioRxiv (Cold Spring Harbor Laboratory). https://doi.org/10.1101/2023.10.04.560808

Yan L, Yang M, Guo H, Yang L et al. Single-cell RNA-Seq profiling of human preimplantation embryos and embryonic stem cells. Nat Struct Mol Biol 2013 Sep;20(9):1131-9. PMID: 23934149

