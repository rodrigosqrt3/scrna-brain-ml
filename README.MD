# Transcriptomic Characterization of CNS Cell Populations in Mus musculus

This repository contains the source code and final monograph developed as a course project for Machine Learning and Statistical Modeling (EST0401), taught by Prof. Dr. Taiane Schaedler Prass at the Graduate Program in Statistics (PPGEst), Federal University of Rio Grande do Sul (UFRGS), in December 2025.

The study analyzes a single-cell RNA-seq (scRNA-seq) dataset containing 3,005 cells from the mouse somatosensory cortex and hippocampus (Zeisel et al., 2015). The main objective of this work was to apply and compare a comprehensive machine learning pipeline—covering dimensionality reduction, unsupervised clustering, differential expression analysis, and supervised classification—for cell-type identification.

> **Note:** The final monograph PDF is written in Portuguese, while the analysis code and repository documentation are structured for broader access.

## Repository Content

* `Villa_2025_Caracterizacao_Transcriptomica.pdf`: Full monograph detailing the mathematical foundations of the methods and biological discussion of the results.
* `analysis.R`: R script containing the complete pipeline, including quality control, normalization, highly variable gene (HVG) selection, dimensionality reduction, clustering algorithms, and classifier training with 5-fold cross-validation.

## Methodology Summary

After selecting the top 2,000 highly variable genes, the analysis compared the performance of several algorithms:

1. Dimensionality Reduction: PCA, t-SNE, UMAP, Kernel PCA (RBF), and Isomap.
2. Clustering: K-means, Hierarchical Clustering (Ward), Gaussian Mixture Models (GMM), Spectral Clustering, and Louvain.
3. Classification: Random Forest, SVM (RBF), XGBoost, Naive Bayes, and KNN.

Detailed performance metrics, confusion matrices, and feature importance analyses (Mean Decrease Gini and Gain) are available in Section 3 of the PDF.

## Requirements

The script was developed using R (version 4.3.0). Main package dependencies include:

* Bioconductor: `scRNAseq`, `SingleCellExperiment`
* Modeling: `xgboost`, `randomForest`, `e1071`, `caret`, `mclust`, `kernlab`, `vegan`, `aricode`, `cluster`
* Visualization: `tidyverse`, `pheatmap`, `Rtsne`, `umap`, `igraph`, `pROC`, `dendextend`, `gridExtra`
