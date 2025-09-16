# Machine Learning for Crystallography (ml-crystals)

This repository contains the R packages, analysis notebooks, datasets, and manuscript for the `ml-crystals`  project. The primary goal is to analyze crystallographic data custom-built tools and to develop predictive models.

## Repository Navigation

This project is organized into three main directories: `packages`, `notebook`, and `paper`.

### 1. `packages/crystract/`

This directory contains the core software developed for this project: the **`crystract` R package**.

*   **Purpose**: `crystract` is a powerful, self-contained R package for parsing Crystallographic Information Files (`.cif`), performing complex geometric calculations (distances, angles), handling crystallographic symmetry, propagating experimental errors, and preparing data for downstream analysis.
*   **Key Contents**:
    *   `R/`: All the R source code for the package's functions.
    *   `DESCRIPTION`: The package metadata, including dependencies.
    *   `vignettes/`: Detailed documentation and tutorials on how to use `crystract`.
    *   `tests/`: A suite of unit tests to ensure the code is working correctly.
*   **How to Use**: This package is the engine for the analysis performed in the `notebook` directory. For detailed installation and usage instructions, please see the **[crystract README](./packages/crystract/README.md)**.

### 2. `notebook/`

This directory contains the main research and analysis workflow. **If you want to understand or reproduce the scientific findings, this is the best place to start.**

*   **Purpose**: To document the entire research process, from raw data to final results and visualizations.
*   **Key Contents**:
    *   `ml-crystals.Rproj`: The RStudio Project file. Open this to work on the analysis.
    *   `project.Rmd`: The main R Markdown notebook showing preliminary work on the development of the functions.
    *   `data/`: Contains the processed data used in the analysis.
    *   `workflow/`: Contains the primary analysis scripts, raw data, and supplementary materials.
        *   `workflow.Rmd`: A detailed R Markdown document (`.html` and `.pdf` renders are also available) that demonstrates the step-by-step application of the `crystract` package.

### 3. `paper/`

This directory contains the manuscript source files for the scientific publication associated with this work.

*   **Purpose**: To house the files necessary for compiling the academic paper.
*   **Key Contents**:
    *   `paper.md`: The manuscript text written in Markdown.
    *   `paper.bib`: The bibliography and references in BibTeX format.

## Data Availability

Please note that due to licensing restrictions, the Crystallographic Information Files (.cif) obtained from the Inorganic Crystal Structure Database (ICSD) are not included in this repository.

To reproduce the analyses presented in the the package examples, you must have access to the ICSD and download the required files yourself. A complete list of all necessary ICSD entry codes can be found in the file:

-   **`required_ICSD_files.csv`**

This file provides the ICSD codes and, for larger datasets, the subdirectory where each file should be placed to replicate the original analysis structure.
