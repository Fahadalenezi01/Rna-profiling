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
  * ---

## 3. Exploratory Data Analysis & Visualization (Python Scripts)
These scripts utilize `pandas`, `seaborn`, and `umap-learn` to process the merged count matrices and generate publication-ready visualizations.

* **Dimensionality Reduction**
  * `static_umap_analysis.py`: Generates a static, labeled 2D UMAP scatter plot using seaborn.
  * `umap_analysis.py`: Performs 3D UMAP reduction and generates an interactive HTML scatter plot using Plotly.
* **Expression Heatmaps**
  * `expression_heatmap_known_markers.py`: Filters unwanted genes and plots heatmaps for predefined blood/saliva markers.
  * `expression_heatmap_pure_vs_mix_group.py`: Generates heatmaps utilizing data-driven sample groupings.
  * `final_heatmap_script.py`: Produces dual heatmaps (top 50 pure sample genes vs. mixture sample markers).
* **Sex Determination**
  * `sex_analysis.py`: Analyzes the expression of the XIST lncRNA across all samples to infer biological sex.
  * `sex_analysis_v2.py`: An updated robust version that evaluates XIST alongside a panel of Y-chromosome markers (e.g., RPS4Y1, DDX3Y).

---

## System Requirements
NanoPlot==1.44.1
kaleido==0.2.1
pandas
seaborn
matplotlib
numpy
umap-learn
plotly
