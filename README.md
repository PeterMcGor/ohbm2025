# Brain Segmentation Analysis Tool

[![DOI](https://zenodo.org/badge/DOI/zenodo.14511324.svg)](https://doi.org/10.5281/zenodo.14511324)
[![Shiny](https://img.shields.io/badge/Shiny-shinyapps.io-blue?logo=R)](YOUR-SHINYAPP-URL)

Interactive visualization tool for comparing traditional and deep learning-based brain segmentation methods.

## Overview

This repository contains the code and analysis tools used in our study comparing FreeSurfer and SynthSeg brain segmentation approaches. The analysis focuses on:
- Volumetric comparisons between traditional and DL-based measurements
- Multi-site reliability analysis
- Within-subject reproducibility
- Between-subject discriminability metrics

## Interactive Application

Access the interactive visualization tool at:  https://petermcgor.shinyapps.io/From-Traditional2DL/

The app provides:
- Software comparison visualizations for each brain region
- Coefficient of variation analysis
- Distance metrics visualization
- Discrimination ratios analysis
- Statistical model summaries

## Repository Structure
```
.
├── app.R                  # Shiny application main file
├── data/                  # Data directory
├── R/                     # R functions
│   ├── analysis.R        # Software comparison analysis
│   ├── data_preparation.R # Data preprocessing
│   └── subject_analysis.R # Subject-level analyses
├── figures/              # Generated figures
└── www/                 # Web assets
```

## Installation

To run the analysis locally:

1. Clone this repository
2. Install R dependencies:
```R
install.packages(c("shiny", "tidyverse", "lme4", "emmeans", 
                  "DT", "patchwork", "ggpubr"))
```
3. Run the Shiny app:
```R
shiny::runApp()
```

## Citation

If you use this code or analysis in your research, please cite:
TBD




