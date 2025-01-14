# BF591 Final Project RShiny App

This repository contains a Shiny dashboard application developed for the BF591 final project. The app provides a graphical user interface for exploring, visualizing, and analyzing biological data, such as sample metadata, gene count data, differential expression results, and GSEA outcomes.

**Project Introduction**:  
In this project, we integrate and analyze RNA-Seq datasets to understand the gene expression patterns and biological pathways associated with diseases such as Huntington’s Disease. Through interactive dashboards, users can filter genes by variance and non-zero counts, generate heatmaps, perform PCA, visualize DE results with volcano plots, and explore GSEA outcomes. This tool helps researchers quickly gain insights into complex RNA-Seq data and streamline their data exploration process.

**Raw Data Source**:  
*2. Post-mortem Huntington’s Disease prefrontal cortex compared with neurologically healthy controls*  
This dataset profiled gene expression with RNA-Seq in post-mortem human dorsolateral prefrontal cortex from patients who died from Huntington’s Disease and age- and sex-matched neurologically healthy controls.  
*mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington’s Disease and neurologically normal individuals*  
[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810)

**Processed Data**:  
The data provided in this repository has undergone a simple preprocessing pipeline to ensure quality and usability in the Shiny app. Some steps (e.g., filtering, normalization) are demonstrated in `processed.R`.

## Features

- **Samples**:  
  - Load and summarize metadata about samples.  
  - View raw sample tables.  
  - Generate diagnostic and distribution plots.

- **Counts**:  
  - Upload normalized gene count data.  
  - Filter data based on variance and non-zero counts.  
  - Visualize diagnostic plots (median counts vs. variance, median counts vs. zeros).  
  - Generate a hierarchical clustering heatmap.  
  - Perform PCA and visualize principal component scores.

- **Differential Expression (DE)**:  
  - Upload DE results (e.g., from DESeq2).  
  - View DE tables with configurable pagination.  
  - Generate a volcano plot with customizable axes and coloring thresholds.  
  - Display a filtered table of significant genes identified in the volcano plot.

- **Gene Set Enrichment Analysis (GSEA)**:  
  - Upload fgsea result files.  
  - Filter and visualize top GSEA results as bar plots and scatter plots.  
  - Interactively view leadingEdge genes for selected pathways.  
  - Download filtered GSEA results.

## Installation & Usage

1. **Clone this repository:**
   ```bash
   git clone https://github.com/temp0jd/BF591finalproject.git
2. **Install required packages (if not already installed):**
   ```bash
   install.packages(c("shiny", "shinydashboard", "dplyr", "ggplot2", "pheatmap", "plotly", "tidyverse","DT", "colourpicker", "shinyWidgets", "fgsea", "gplots", "RColorBrewer", "rio","ggbeeswarm", "GSEABase", "org.Hs.eg.db", "biomaRt", "stringr"))
3.**Run the Shiny app:**
