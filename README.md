# scRNA-seq Analysis using Seurat

## Overview
This repository contains an analysis pipeline for single-cell RNA sequencing (scRNA-seq) data using the Seurat package in R. The workflow includes data loading, quality control, normalization, dimensionality reduction, integration, clustering, and visualization.

---

## Dataset
Two datasets were used in this analysis:

- **DS3**
- **DS1**

Both datasets were processed using the 10X Genomics format.

---

## Workflow

The analysis follows these main steps:

1. **Data Loading**
   - Import count matrices using `Read10X`
   - Create Seurat objects

2. **Data Merging**
   - Combine multiple datasets into a single object

3. **Quality Control**
   - Calculate mitochondrial gene percentage
   - Filter low-quality cells

4. **Normalization**
   - Log-normalization of gene expression data

5. **Feature Selection**
   - Identify highly variable genes

6. **Scaling**
   - Scale the data for downstream analysis

7. **Dimensionality Reduction**
   - Principal Component Analysis (PCA)
   - UMAP visualization

8. **Integration**
   - Batch correction using Harmony

9. **Clustering**
   - Graph-based clustering of cells

10. **Differential Expression**
    - Identification of cluster-specific marker genes

11. **Visualization**
    - Violin plots
    - Scatter plots
    - UMAP plots
    - Feature plots

---

## Key Plots Generated

- QC Violin Plots (`VlnPlot`)
- QC Scatter Plots (`FeatureScatter`)
- PCA Elbow Plot (`ElbowPlot`)
- UMAP Visualization (`DimPlot`)
- Feature Expression Plots (`FeaturePlot`)
- Marker Gene Violin Plots (`VlnPlot`)

---

## Requirements

Install the following R packages before running the code:

```r
install.packages("Seurat")
install.packages("patchwork")
install.packages("dplyr")
install.packages("harmony")
