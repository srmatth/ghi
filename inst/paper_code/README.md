# Paper Code for GHI Regression Analysis

This directory contains all code used to produce the analyses and results presented in our paper on Golden Health Index (GHI) Regression. We provide two pathways for reproducing our work:

1. **For reviewers and users without NACC data access**: Use our simulated data walkthrough
2. **For users with NACC data access**: Full reproduction using real data

---

## Quick Start: Reproducing Results with Simulated Data

**This is the recommended path for reviewers.**

We provide simulated data that matches the structure and statistical properties of the original NACC dataset. This allows you to verify our code and methodology without requiring data access.

### Instructions

1. Open **[Paper_Reproduction_Simulated.Rmd](Paper_Reproduction_Simulated.Rmd)** in RStudio
2. Scroll down to the `configure-paths` chunk and make sure that the paths point to the simulated dataset (`NACC_mod_dat_simulated.csv`)
3. Click "Knit" to run the entire analysis, or run chunks interactively
4. The notebook walks through every table and figure in the paper

The `Paper_Reproduction_Simulated.Rmd` notebook:
- Uses simulated data generated from our fitted models
- Reproduces Table 1 (descriptive statistics), Table 2 (model coefficients), and all figures
- Includes model diagnostics and sensitivity analyses
- Can be run locally on a standard laptop (no server required)
- Takes approximately 30 minutes to run completely

**Note:** Results will differ numerically from the paper (which uses real NACC data), but the methods and code structure are identical.

---

## Full Reproduction with Real NACC Data

**For researchers who have obtained their own NACC data access.**

See **[Paper_Reproduction_NACC.md](Paper_Reproduction_NACC.md)** for complete, step-by-step instructions including:
- Prerequisites and required R packages
- Data preparation workflow
- Server/HPC requirements for computationally intensive steps
- Execution order and dependencies between scripts

### Summary

1. Obtain NACC data through the official request process at [naccdata.org](https://naccdata.org)
2. Run `03_create_NACC_model_data.R` to prepare the analysis dataset
3. Run `05_fit_NACC_model.R` to fit the two-part model
4. Run analysis scripts (08-14) to generate results
5. Optionally run simulation study (01-02) and bootstrap variance estimation (06-07) on a server

---

## File Descriptions

### Data Simulation

| File | Description | Requires Real Data |
|------|-------------|-------------------|
| `15_simulate_from_fitted_models.R` | Generates simulated data from fitted NACC models. Creates `NACC_mod_dat_simulated.csv` with covariate distributions matching the original data. | Yes (to generate; not to use pre-generated data) |

### Simulation Study (Server)

| File | Description | Requires Real Data |
|------|-------------|-------------------|
| `01_code_for_simulation_on_server.R` | Main simulation study script. Generates simulation study data, fits two-part model, computes bootstrap variance estimates. Designed to run on HPC cluster. | No |
| `02_run_simulation_parallel.sh` | Bash wrapper to run simulation 50 times in parallel on a server. Creates output in `ghi_sim_res/`. | No |

### NACC Data Preparation

| File | Description | Requires Real Data |
|------|-------------|-------------------|
| `03_create_NACC_model_data.R` | Subsets, cleans, and prepares raw NACC database for analysis. Creates `NACC_mod_dat.csv`. | Yes |

### Data Exploration

| File | Description | Requires Real Data |
|------|-------------|-------------------|
| `04_NACC_data_exploration.R` | Exploratory data analysis and descriptive statistics. Generates summary tables and figures. Falls back to simulated data if real data unavailable. | No (uses simulated if available) |

### Model Fitting

| File | Description | Requires Real Data |
|------|-------------|-------------------|
| `05_fit_NACC_model.R` | Fits the two-part model: GLM for dementia probability (Part 1) and LTRC survival model (Part 2). Saves fitted models to `.rds` files. | No (uses simulated if available) |

### Bootstrap Variance Estimation (Server)

| File | Description | Requires Real Data |
|------|-------------|-------------------|
| `06_NACC_bootstrap_on_server.R` | Bootstrap resampling for variance estimation. Computationally intensive; designed for HPC. | Yes |
| `07_run_NACC_bootstrap_parallel.sh` | Bash wrapper to run bootstrap 50 times in parallel on a server. | Yes |

### Model Analysis

| File | Description | Requires Real Data |
|------|-------------|-------------------|
| `08_analyze_NACC_model.R` | Analyzes fitted models. Generates coefficient tables (Table 2), prediction figures, and model summaries. | No (uses simulated if available) |
| `09_additional_NACC_analyses.R` | Supplementary analyses including marital status changes. Some sections require full NACC data. | Partially |
| `10_reviewer_response_code.R` | Additional analyses responding to reviewer feedback. Requires full longitudinal NACC data. | Yes |
| `11_historical_analysis_function.R` | Helper functions for historical/retrospective analysis (analyzing data as if at a past date). | Yes |
| `12_weighted_analysis.R` | Weighted sampling analysis to adjust for APOE genotype oversampling in NACC. | Yes |
| `13_simulation_for_weighted_sampling.R` | Simulation study validating the weighted sampling approach. | No |
| `14_additional_model_checking.R` | Model diagnostics including surrogate residuals and calibration plots. | No (uses simulated if available) |

### Documentation

| File | Description |
|------|-------------|
| `Paper_Reproduction_Simulated.Rmd` | Interactive R Markdown notebook reproducing all paper results with simulated data |
| `Paper_Reproduction_NACC.md` | Detailed workflow documentation for full reproduction with real data |
| `README.md` | This file |

---

## Troubleshooting

### Package Installation

The following R packages are required:

```r
# Core packages
install.packages(c("dplyr", "readr", "ggplot2", "patchwork", "fs", "data.table"))

# Survival analysis
install.packages("survival")

# Additional packages for analyses
install.packages(c("lubridate", "stringr", "xtable", "lme4",
                   "truncreg", "truncnorm", "statmod", "sandwich",
                   "lmtest", "ResourceSelection"))

# ltrc package (for left-truncated right-censored survival models)
devtools::install_github("srmatth/ltrc")
```

### Common Issues

**"NACC data file not found"**
- If you're running with simulated data, this is expected for scripts that require real data
- Scripts 04, 05, 08, and 14 will automatically fall back to simulated data
- Scripts 03, 06, 07, 09, 10, 11, and 12 require real NACC data and will skip or error

**"Object 'get_clean_model' not found"**
- This function is provided by the `ltrc` package
- Ensure `library(ltrc)` is loaded before running model fitting scripts

**Memory errors during bootstrap**
- Bootstrap scripts (01, 06) are designed for servers with 16+ GB RAM
- Reduce `n_bootstrap` parameter for local testing
- Consider running fewer parallel jobs

**Long runtime**
- Scripts 01, 02, 06, 07 are computationally intensive (hours on server)
- For quick validation, skip these scripts and use pre-computed results
- The `Paper_Reproduction_Simulated.Rmd` notebook avoids intensive computations

### Getting Help

If you encounter issues:
1. Check that all required packages are installed
2. Verify file paths match your directory structure
3. For NACC data issues, consult the NACC data dictionary
4. Open an issue in the repository with error messages and session info (`sessionInfo()`)

---

## Citation

If you use this code, please cite our paper:

[Citation information to be added upon publication]
