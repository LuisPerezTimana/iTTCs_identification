# scRNA-seq analysis for identification and characterization of transient totipotent-like cells during somatic cell reprogramming

This project focuses on the preprocessing and integration of single-cell RNA sequencing (scRNA-seq) data. We analyzed data from multiple sources (Table 1) to generate integrated visualizations such as UMAP and Dotplot for specific gene expressions.

## Project Structure

The project consists of the following  main scripts:
- **Nairetal_Cluster39_GSEA_ash.R**: performs dimensional reduction (UMAP) of Nair et al datases, performs sub clustering and isolates the cluster 39 for further integration
- **preprocessing.R**: Preprocesses and cleans the input data from different sources.
- **Integration.R**: Integrates the preprocessed data, performs dimensional reduction (UMAP), and generates visualizations of specific gene expressions.

## Data

| **Sample name**                          | **Description**                                                                 | **Accession number**       | **Citation**                                                                                                                                                                                                                            |
|------------------------------------------------|---------------------------------------------------------------------------------|----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Cluster 39                                    | Totipotent-like cluster obtained from day 14 single cell RNA-seq               | GEO: GSE242424             | Nair et al. (2023). *Transcription factor stoichiometry, motif affinity and syntax regulate single-cell chromatin dynamics during fibroblast reprogramming to pluripotency.* bioRxiv [Preprint]. doi: [10.1101/2023.10.04.560808](https://doi.org/10.1101/2023.10.04.560808). |
| Early preimplantation human development dataset | scRNA-seq of early preimplantation human embryos                                | GEO: GSE36552              | Yan et al. (2013). *Single-cell RNA-Seq profiling of human preimplantation embryos and embryonic stem cells.* Nat Struct Mol Biol. doi: [10.1038/nsmb.2660](https://doi.org/10.1038/nsmb.2660). PMID: 23934149.                                                        |
| induced blastomere-like (iBM) cells [10X_12_possorted_genome_bam.bam] | Single cell RNA seq of iBM at 12h post 15min treatment with doxycycline to induce DUX4 in DUX4-TetOn hESCs | ARRAY_EXPRESS: E-MTAB-10581 | Yoshihara et al. (2022). *Transient DUX4 expression in human embryonic stem cells induces blastomere-like expression program that is marked by SLC34A2.* Stem Cell Reports. doi: [10.1016/j.stemcr.2022.06.002](https://doi.org/10.1016/j.stemcr.2022.06.002). |

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

Nair S, Ameen M, Sundaram L, Pampari A, Schreiber J, Balsubramani A, Wang YX, Burns D, Blau HM, Karakikes I, Wang KC, Kundaje A. Transcription factor stoichiometry, motif affinity and syntax regulate single-cell chromatin dynamics during fibroblast reprogramming to pluripotency. bioRxiv [Preprint]. 2023 Oct 21:2023.10.04.560808. doi: 10.1101/2023.10.04.560808. PMID: 37873116; PMCID: PMC10592962.

Yan L, Yang M, Guo H, Yang L, Wu J, Li R, Liu P, Lian Y, Zheng X, Yan J, Huang J, Li M, Wu X, Wen L, Lao K, Li R, Qiao J, Tang F. Single-cell RNA-Seq profiling of human preimplantation embryos and embryonic stem cells. Nat Struct Mol Biol. 2013 Sep;20(9):1131-9. doi: 10.1038/nsmb.2660. Epub 2013 Aug 11. PMID: 23934149.			

Yoshihara M, Kirjanov I, Nykänen S, Sokka J, Weltner J, Lundin K, Gawriyski L, Jouhilahti EM, Varjosalo M, Tervaniemi MH, Otonkoski T, Trokovic R, Katayama S, Vuoristo S, Kere J. Transient DUX4 expression in human embryonic stem cells induces blastomere-like expression program that is marked by SLC34A2. Stem Cell Reports. 2022 Jul 12;17(7):1743-1756. doi: 10.1016/j.stemcr.2022.06.002. Epub 2022 Jun 30. PMID: 35777358; PMCID: PMC9287684.			

