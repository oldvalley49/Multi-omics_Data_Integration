# Benchmarking Algorithms for Multi-omics Single Cell Data Integration Across Unpaired scRNA-seq and scATAC-seq Data.
**Tomoya Furutani**

This repository contains codes and visualizations for benchmarking cross-modality integration across different datasets. Four integration algorithms were evaluated across multiple publicly available datasets. 

## Project Background

Recent advancements in single-cell sequencing technologies have revolutionized biological research; these technologies enable the analysis of cellular heterogeneity and regulatory networks by generating large datasets across various modalities. scRNA-seq provides insights into gene expression by measuring mRNA abundance, while scATAC-seq focuses on open chromatin regions, helping our understanding of regulatory elements. 

Despite the recent advancements in experimental methods that allow for simultaneous measurement of multiple omics modalities at single-cell resolution, most existing datasets only include one modality, presenting challenges for multi-omics integration. A significant issue in the integration of single-modality datasets is that these omics layers have different feature spaces; for scRNA-seq, they are genes while for scATAC-seq, they are peaks corresponding to open-chromatin regions. Therefore they cannot be directly compared. 

However, integration of single-modality data, especially between scRNA-seq and scATAC-seq, allows researchers to identify cell-specific regulatory networks, uncover new cell clusters, and reveal potential biomarkers while taking advantage of the pre-existing single-modality datasets. Thus, growing number of computational tools have been developed to facilitate this process by aligning datasets from two different modalities into a unified latent space. 

In this project, I benchmark four major computational tools developed for scRNA-seq and scATAC-seq integration(Seurat, LIGER, bindSC, and GLUE) across multiple datasets and conditions to evaluate their ability to align the two modalities. Moreover, I will discuss possible ways in which these methods could be improved. 

## Project Walkthrough

### Overview of Benchmarking Strategy

In order to benchmark the performance of these integration algorithms, I used publicly available multimodal datasets that simultaneously profiled gene expression and chromatin accessibility. Then, I treat the datasets for two different modalities as originating from two different single-modality experiments and integrate them using the algorithms. Since these were originally sampled from the same cells, we have the ground truth for cell-to-cell correspondence that allows us to measure the accuracy of the alignment. 

To quantitvely measure the alignment accuracy of the two modalities, I used FOSCTTM(fraction of cells closer than the true match) which was first introduced by [Liu et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8496402/). FOSCTTM measures alignment by first calculating the distance matrix beween the data points originating from scRNA-seq and scATAC-seq in the coembedded space. Then, we calculate the proportion of cells of different modality that are closer to the cell than the ground truth-matched cells and take average across all cells measured. 
### Data Availability
### Pre-processing
### Integration
### Benchmarking
### Conclusion and Future Improvements

