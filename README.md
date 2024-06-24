Integrative analysis of RNA-Seq and metabolic data.

This repository contains the code and data required for time series analysis and network anlaysis of RNA-Seq data and metabolic data. The analyses code was written in R and contain three parts.

(1) Time_series_analysis.r: perform time series analysis for genes/metabolites using R package maSigpro and identify three patterns.

(2) Pathway_analysis_of_three_clusters.r: perform pathway analysis for three clusters of genes from step #1 and draw dot plot.

(3) Network_analysis.r: combine RNA-Seq data and metabolic data to built a correlation network. Calcualte a normalized connectivity score to prioritize genes and metabolites.


Dependencies

R packages maSigPro, ComplexHeatmap, gplots, RColorBrewer, pROC, colorspace, ggplot2, stringr are needed for the analyses.


Data

All required data could be found in the 'Data' folder.
