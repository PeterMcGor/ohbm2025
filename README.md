# From-Traditional2DL Brain Segmentation analysis

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

@inproceedings{gordaliza2025traditional,
  title     = {From Traditional to Deep Learning: Multi-Site Reliability and Subject Differentiation in Brain Segmentation},
  author    = {Gordaliza, Pedro M. and Alem{\'a}n-G{\'o}mez, Yasser and Banus, Jaume and Bach Cuadra, Meritxell},
  booktitle = {Organization for Human Brain Mapping (OHBM) Annual Meeting},
  address   = {Brisbane, Australia},
  year      = {2025},
  month     = jun,
  doi       = {10.5281/zenodo.15019175},
  url       = {https://doi.org/10.5281/zenodo.15019175}
}




