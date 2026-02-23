# RNA Profiling pipeline for Saliva and Blood Mixture Analysis

This repository contains the complete bioinformatics and statistical analysis pipeline used for my PhD research. The project focuses on the RNA profiling of pure saliva, pure blood, and mixed physiological samples using Oxford Nanopore sequencing.

## Project Overview
The pipeline processes raw Nanopore sequencing data through basecalling, alignment, transcript quantification, differential expression analysis, and the identification of novel long non-coding RNAs (lncRNAs). It also includes a robust R-based statistical suite for identifying tissue-specific marker genes and deconvoluting mixed samples.

---

## 1. Bioinformatics Pipeline (Bash Scripts)
These scripts are designed to run in a Linux/WSL environment.

* **Basecalling & Quality Control**
  * `run_basecalling.sh`: Automates raw data basecalling using Dorado (v1.0.2).
  * `run_rebasecalling.sh`: Targeted re-basecalling for specific problem samples.
  * `qc_alignment.sh` & `run_analysis.sh`: Performs alignment using Minimap2, sorts/indexes with Samtools, and generates NanoPlot QC reports.
* **Quantification & lncRNA Discovery**
  * `stringtie_analysis.sh`: Assembles and quantifies transcripts using StringTie.
  * `run_lncRNA_prep.sh`: A comprehensive pipeline combining StringTie, FEELnc, and custom filtering to identify and quantify novel non-coding RNAs.
* **Documentation**
  * `how to align sort and index and split bam files`: Text notes and reference commands for manual SAMtools/Minimap2 operations.

---

## 2. Statistical Analysis & Visualization (R Scripts)
These scripts process the TPM (Transcripts Per Million) outputs from StringTie to generate thesis-ready plots and statistical models.

* **Dimensionality Reduction & Clustering**
  * `run_pca_all_samples.R` / `run_pca_mixed_samples.R` / `run_pca_plot.R`: Various implementations of Principal Component Analysis using DESeq2 transformations.
  * `UMAP_analysis.R`: Non-linear dimensionality reduction for sample relationship visualization.
  * `PCA_mixed_pure.R`: Targeted PCA utilizing only top identified marker genes.
* **Differential Expression & Marker Identification**
  * `OneVsRest.R`: A robust Limma-based pipeline for one-vs-rest differential expression to identify strict tissue-specific markers.
  * `wilcox_rank_sum.R`: Non-parametric marker identification comparing pure saliva vs. pure blood.
* **Expression Profiling & Visualization**
  * `run_profiling_heatmap.R`: Generates relative expression heatmaps for top body fluid-specific genes.
  * `run_profiling_analysis.R`: Exploratory analysis and heatmaps based on the top 50 most variable genes.
  * `run_log_tpm_profiles.R` & `visualize_tpm.R`: Bar plot visualizations of highly expressed genes, mapped to human gene symbols.
  * `plot_stringtie_output.R`: Utilities for basic TPM plotting and filtering.
* **Mixture Deconvolution**
  * `mixed_samples_analysis.R`: Calculates normalized signature scores to determine the proportional breakdown of mixed samples based on pure-tissue markers.

---

## System Requirements
* **CLI Tools:** `samtools`, `minimap2`, `stringtie`, `dorado` (v1.0.2), `FEELnc`
* **Python (3.12+):** `NanoPlot==1.44.1`, `kaleido==0.2.1`
* **R (4.0+):** `DESeq2`, `limma`, `pheatmap`, `ggplot2`, `dplyr`, `rtracklayer`, `uwot`, `org.Hs.eg.db`
