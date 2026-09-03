

This repository contains the data, code, and reproducible analysis website for the manuscript:

**Declining environmental predictability erodes reproduction in a tropical seabird**



## Overview

This repository documents the full computational workflow used in the study, including:

1. Detection of chlorophyll-a bloom onset.
2. Calculation of phenological mismatch between seabird breeding phenology and marine productivity.
3. Bayesian model fitting for bloom phenology, mismatch, reproductive success, and environmental predictability.
4. Model diagnostics and model comparison.
5. Full model summaries.
6. Sensitivity analyses.
7. Code used to generate manuscript figures.

The analysis is organized as a Quarto reproducibility website. The rendered website can be inspected directly in a web browser, while the source `.qmd` files can be re-rendered locally.

## Repository structure

```text
.
├── book/        # Quarto source files for the reproducible analysis website
├── data/        # Data files used by the analysis workflow
├── docs/        # Rendered website files
├── Plots/       # Output figures
├── icons/       # Website icons
├── tables.xlsx  # Supplementary/output tables
├── LICENSE      # MIT license
└── README.md
```

The main analysis files are:

```text
book/1-bloom-onset.qmd              # Chlorophyll-a bloom onset detection
book/2-phenological-mismatches.qmd  # Phenological mismatch calculation
book/3-model-fitting.qmd            # Bayesian model-fitting workflow
book/4-diagnostics.qmd              # Diagnostics and model selection
book/6-full_summaries.qmd           # Model summaries
book/sensitivity_blooms.qmd         # Bloom-onset sensitivity analyses
book/code_figures.qmd               # Figure code
book/session-info.qmd               # Software/session information
```

## 1. System requirements

### Operating system and tested software environment

The rendered reproducibility website was generated using the following system:

```text
R version: 4.5.2 (2025-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Operating system: Windows 11 x64, build 26200
LAPACK version: 3.12.1
Time zone: America/Edmonton
Locale: English_Guernsey.utf8
Quarto version: 1.8.25
```

The workflow should also run on macOS or Linux if R, Quarto, CmdStan, and the required R packages are installed.

### Main software dependencies

The workflow requires:

- R
- Quarto
- R packages listed below
- CmdStan/CmdStanR only if users want to re-fit the Bayesian models

### Main R package versions

The analysis uses R, Quarto, and the following main R packages.

Packages attached in the bloom-onset workflow include:

```text
DT 0.34.0
zoo 1.8-15
here 1.0.2
patchwork 1.3.2
cowplot 1.2.0
ggspatial 1.1.10
rnaturalearth 1.2.0
sf 1.0-24
lubridate 1.9.4
forcats 1.0.1
stringr 1.6.0
dplyr 1.1.4
purrr 1.2.1
readr 2.1.6
tidyr 1.3.2
tibble 3.3.1
ggplot2 4.0.1
tidyverse 2.0.0
```

Packages attached in the model-fitting workflow include:

```text
DT 0.34.0
loo 2.9.0
brms 2.23.0
Rcpp 1.1.1
here 1.0.2
lubridate 1.9.4
forcats 1.0.1
stringr 1.6.0
dplyr 1.1.4
purrr 1.2.1
readr 2.1.6
tidyr 1.3.2
tibble 3.3.1
ggplot2 4.0.1
tidyverse 2.0.0
```

Additional packages loaded through dependencies include:

```text
bayesplot 1.15.0
posterior 1.6.1
RcppParallel 5.1.11-1
rmarkdown 2.30
knitr 1.51
xfun 0.56
htmlwidgets 1.6.4
yaml 2.3.12
jsonlite 2.0.0
Matrix 1.7-4
nlme 3.1-168
MASS 7.3-65
emmeans 2.0.1
bridgesampling 1.2-1
rstantools 2.6.0
```

The complete package/session information is printed at the end of the rendered analysis chapters.

### Bayesian model-fitting backend

Bayesian models were fitted using `brms` with the `cmdstanr` backend.

The model-fitting workflow sets:

```r
options(mc.cores = parallel::detectCores(), brms.backend = "cmdstanr")
```

A working CmdStan installation and a C++ toolchain are required to re-fit the Bayesian models.

The rendered website records the use of the `cmdstanr` backend, but it does not record the exact CmdStan or `cmdstanr` version used at the time of fitting. Users who re-run the models should record these versions with:

```r
cmdstanr::cmdstan_version()
packageVersion("cmdstanr")
```

To install and check CmdStan from R:

```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

cmdstanr::install_cmdstan()
cmdstanr::check_cmdstan_toolchain()
cmdstanr::cmdstan_version()
```

### Non-standard hardware

No non-standard hardware is required.

However, full Bayesian model re-fitting is computationally intensive. A multi-core desktop or workstation is recommended if users want to re-fit all Bayesian models.

Viewing the rendered website, rendering lightweight chapters, and inspecting the code do not require special hardware.

## 2. Installation guide

### Step 1: Clone or download the repository

Using Git:

```bash
git clone https://github.com/Santiago-0rtega/Collapse_environmental_predictability_erodes_reproductive_success.git
cd Collapse_environmental_predictability_erodes_reproductive_success
```

Alternatively, download the repository as a `.zip` file from GitHub and unzip it.

### Step 2: Install R and Quarto

Install R from:

https://cran.r-project.org/

Install Quarto from:

https://quarto.org/docs/get-started/

This project was rendered with:

```text
Quarto 1.8.25
```

### Step 3: Install R packages

From R, run:

```r
install.packages(c(
  "tidyverse",
  "here",
  "DT",
  "zoo",
  "sf",
  "rnaturalearth",
  "ggspatial",
  "cowplot",
  "patchwork",
  "brms",
  "loo",
  "posterior",
  "bayesplot",
  "sessioninfo"
))
```

Install `cmdstanr` separately if you plan to re-fit Bayesian models:

```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

cmdstanr::install_cmdstan()
cmdstanr::check_cmdstan_toolchain()
```

### Typical installation time

On a normal desktop computer:

```text
R package installation: approximately 10–30 minutes
CmdStan installation: approximately 20–60 minutes
Total installation time: approximately 30–90 minutes
```

Installation time depends on the operating system, internet connection, and whether R packages need to be compiled from source.

## 3. Demo

### Demo option A: view the rendered website

The fastest demo is to view the rendered reproducibility website:

https://santiago-0rtega.github.io/Collapse_environmental_predictability_erodes_reproductive_success/

Expected output:

```text
A rendered website containing the full analysis workflow, including bloom-onset detection, phenological mismatch calculation, model-fitting code, model diagnostics, model summaries, sensitivity analyses, and figure code.
```

Expected run time:

```text
Immediate in a web browser.
```

### Demo option B: render the website locally

From the repository root, run:

```bash
quarto render book
```

Expected output:

```text
A local version of the reproducible analysis website.
```

Expected run time:

```text
A few minutes if using the existing rendered/frozen outputs.
Longer if computational chunks are re-executed.
```

### Demo option C: run a lightweight executable chapter

A lightweight executable demo is the bloom-onset chapter:

```bash
quarto render book/1-bloom-onset.qmd
```

Expected output:

```text
An HTML file documenting chlorophyll-a bloom onset detection, including data loading, seasonal-window definitions, threshold-based onset estimation, summary tables, and bloom-onset visualizations.
```

Expected run time on a normal desktop computer:

```text
Approximately 1–10 minutes.
```

### Demo option D: inspect the model-fitting workflow

The Bayesian model-fitting workflow can be inspected in:

```text
book/3-model-fitting.qmd
```

or through the rendered page:

https://santiago-0rtega.github.io/Collapse_environmental_predictability_erodes_reproductive_success/3-model-fitting.html

Expected output:

```text
A documented model-fitting workflow showing the models used for bloom phenology, phenological mismatch, reproductive success, and environmental predictability.
```

Expected run time:

```text
Immediate when viewing the rendered website.
Potentially hours or longer if users manually re-enable and re-fit the Bayesian models.
```

## 4. Instructions for use

### Running the supplied workflow

The workflow is organized as Quarto chapters in the `book/` folder. The recommended order is:

```text
1. book/1-bloom-onset.qmd
2. book/2-phenological-mismatches.qmd
3. book/3-model-fitting.qmd
4. book/4-diagnostics.qmd
5. book/6-full_summaries.qmd
6. book/sensitivity_blooms.qmd
7. book/code_figures.qmd
```

To render the full website:

```bash
quarto render book
```

To render one chapter at a time:

```bash
quarto render book/1-bloom-onset.qmd
quarto render book/2-phenological-mismatches.qmd
quarto render book/3-model-fitting.qmd
quarto render book/4-diagnostics.qmd
quarto render book/6-full_summaries.qmd
quarto render book/sensitivity_blooms.qmd
quarto render book/code_figures.qmd
```

### Model-fitting note

The model-fitting chapter documents the model-fitting workflow used in the study. It includes Gaussian location-scale models for bloom phenology and hierarchical Bayesian models for phenological mismatch and reproductive performance.

Most model-fitting code chunks are included as a transparent workflow record and are not automatically executed during rendering of the supplementary website. This avoids unintentionally re-fitting computationally intensive Bayesian models when the website is rendered.

Users who want to re-fit the Bayesian models should inspect:

```text
book/3-model-fitting.qmd
```

Before re-fitting models, users should confirm that CmdStan is installed correctly:

```r
cmdstanr::check_cmdstan_toolchain()
cmdstanr::cmdstan_version()
```

Then, users can manually enable or run the relevant model-fitting chunks in `book/3-model-fitting.qmd`.

### Running the workflow on new data

To apply the workflow to new data:

1. Format the new data to match the column names and variable definitions used in the files in `data/`.
2. Place the new data files in the `data/` folder, or modify the file paths in the relevant `.qmd` files.
3. Check that paths are relative to the project root and use `here::here()` where appropriate.
4. Re-run the relevant Quarto chapter.
5. Inspect model diagnostics before interpreting model outputs.

The main files to adapt are:

```text
book/1-bloom-onset.qmd              # Chlorophyll-a bloom detection
book/2-phenological-mismatches.qmd  # Mismatch calculation
book/3-model-fitting.qmd            # Bayesian model fitting
book/4-diagnostics.qmd              # Diagnostics and model checking
book/code_figures.qmd               # Figure generation
```

## Optional: full reproduction instructions

To reproduce the full analysis:

1. Clone or download this repository.

2. Install R, Quarto, CmdStan, and the required R packages.

3. Open the repository root as the working directory.

4. Confirm that the data files are available in `data/`.

5. Render the Quarto website:

```bash
quarto render book
```

6. If full Bayesian model re-fitting is required, inspect `book/3-model-fitting.qmd`, check that `cmdstanr` is correctly installed, and manually enable or run the relevant model-fitting chunks.

7. After model fitting, render the diagnostics and model-summary chapters:

```bash
quarto render book/4-diagnostics.qmd
quarto render book/6-full_summaries.qmd
```

8. Recreate figures using:

```bash
quarto render book/code_figures.qmd
```

9. Record the final software environment:

```r
utils::sessionInfo()
```

or, for a more detailed report:

```r
sessioninfo::session_info()
```

For Bayesian model re-fitting, also record:

```r
cmdstanr::cmdstan_version()
packageVersion("cmdstanr")
```

## Expected outputs from full reproduction

The full workflow produces:

```text
1. Rendered Quarto HTML pages documenting the analysis.
2. Bloom-onset estimates.
3. Phenological mismatch estimates.
4. Bayesian model objects or summaries, depending on which chunks are enabled.
5. Model diagnostics and model-comparison outputs.
6. Full model-summary tables.
7. Sensitivity-analysis outputs.
8. Manuscript figures.
```

## Notes on computational reproducibility

The rendered website is intended to provide a transparent and inspectable record of the full analysis. Some computationally intensive Bayesian model-fitting chunks are not automatically executed during website rendering. This design prevents accidental re-fitting of long-running models while preserving the model code and workflow documentation.

Users who want exact full re-fitting should use the code in `book/3-model-fitting.qmd`, confirm their CmdStan installation, and record their local CmdStan and `cmdstanr` versions.

## License

This repository is released under the MIT License. See `LICENSE` for details.
