Integrative analysis of RNA-seq and metabolomics data.

This repository contains the code and data required for time series analysis and network analysis of RNA-seq data and metabolomics data. The analysis codes were written in R and contain three parts.

(1) Time_series_analysis.r: perform time series analysis for genes/metabolites using R package maSigpro to identify three clusters.

(2) Pathway_analysis_of_three_clusters.r: perform pathway analysis for three clusters of genes from step #1 and draw dot plots.

(3) Network_analysis.r: combine RNA-seq data and metabolomics data to generate a correlation network, and to calculate a normalized connectivity score to prioritize genes and metabolites.

Dependencies

R packages maSigPro, ComplexHeatmap, gplots, RColorBrewer, pROC, colorspace, ggplot2, stringr are needed for the analyses.

Data

All required data can be found in the 'Data' folder.
