# Reproducibility Workflow for Paper Results

This document maps the complete workflow for reproducing the paper results, including dependencies, inputs, outputs, and computational requirements.

---

## Quick Start: Choose Your Path

**Are you a reviewer without access to real NACC data?**
→ See **[Paper_Reproduction_Simulated.Rmd](Paper_Reproduction_Simulated.Rmd)** for step-by-step instructions to test our code with simulated data.

**Do you have access to the original NACC database?**
→ Continue reading this document for full reproduction instructions.

---

## Overview

The workflow consists of three main sections:
1. **Simulation Study** (Local or Server) - Generates simulated data and estimates
2. **NACC Real Data Analysis** (Server-intensive) - Fits models and bootstraps variance estimates
3. **Analysis & Visualization** (Local) - Generates figures and tables

---

## SECTION 1: SIMULATION STUDY

### Stage 1a: Simulation with Bootstrap Variance (Server-Required)

#### Script: `01_code_for_simulation_on_server.R`

**Purpose:** Generates simulated data, fits two-part model, and computes bootstrap variance estimates.

**Input:**
- Command-line argument: `FILE_NUM` (integer, e.g., 1, 2, 3, ...)
- No external data files required

**Output:**
- `ghi_sim_res/error_a_{SAMPLE_SIZE}_{FILE_NUM}.rds`
  - Large R object containing:
    - Simulated data
    - Fitted Part 1 model (GLM for disease probability)
    - Fitted Part 2 model (LTRC survival model)
    - Prediction results
    - Bootstrap results (300 iterations per FILE_NUM)

**Parameters (hard-coded in script):**
- `N_SIM` = 5 (simulation iterations per FILE_NUM)
- `SAMPLE_SIZE` = 400
- `n_bootstrap` = 300 (per iteration)
- True parameters: `alpha`, `beta`, `omega`
- Covariate grid for predictions: t, x_1, x_2, l

**Computational Requirements:**
- ⚠️ **VERY HIGH** - Must run on server/HPC
- Bootstrap section is especially computationally intensive
- Estimated time: 1-2 hours per FILE_NUM on modern CPU

**How to Run (Server):**
```bash
# Run a single file
Rscript 01_code_for_simulation_on_server.R 1

# Or use the batch runner script (see 1b below)
```

---

### Stage 1b: Batch Parallel Execution (Server-Required)

#### Script: `02_run_simulation_parallel.sh`

**Purpose:** Wrapper script to run simulation 50 times in parallel for computational efficiency.

**Input:**
- `01_code_for_simulation_on_server.R` (the main simulation script)

**Output:**
- 50 `.rds` files: `ghi_sim_res/error_a_400_1.rds` through `ghi_sim_res/error_a_400_50.rds`
- 50 log files: `out_1.txt` through `out_50.txt` (for error tracking)

**Parameters:**
- Loop: `i in {1..50}`
- Each iteration runs in background (`&`)
- Output logs: `out_{i}.txt`

**Computational Requirements:**
- ⚠️ **VERY HIGH** - Designed for cluster/server with parallel job support
- Uses 50 parallel jobs (adjust loop if resources limited)
- Total runtime: ~1-2 hours wall time (depending on node availability)

**How to Run (Server):**
```bash
cd inst/paper_code
bash 02_run_simulation_parallel.sh
# Monitor: tail -f out_*.txt
# Wait for all background jobs to complete: wait
```

**Output Location:** Creates directory `ghi_sim_res/` if not exists

---

## SECTION 2: NACC REAL DATA ANALYSIS

### Stage 2a: Data Preparation

#### Script: `03_create_NACC_model_data.R`

**Purpose:** Subsets, cleans, and prepares NACC database for analysis.

**Input:**
- External file: `~/Documents/Data/NACC/2025-01-22_UDS_genotype.csv`
  - This file is NOT in the repository (external data source)
  - Must be obtained directly from NACC

**Output:**
- `inst/extdata/NACC_mod_dat.csv`
  - Clean, analysis-ready dataset
  - **Key variables:** id, got_dementia, t (age at death - 65), l (age at entry - 65), s (age at dementia - 65), y (s/t if dementia else 1), mod_y, mod_l, covariates (years_education, is_female, is_married, comorbidity, is_race_black, is_race_other, num_e4)
  - Filters: Age ≥ 65 at baseline, dementia-free at entry, completed follow-up

**Data Filtering Steps:**
1. Subset columns from full NACC database
2. Filter for deaths only
3. Filter for age ≥ 65
4. Filter for dementia-free at baseline (first visit)
5. Drop missing education values
6. Create derived variables (dates, ages, indicators)

**Computational Requirements:**
- ✅ **LOW** - Can run locally
- Estimated time: < 5 minutes

**How to Run (Local):**
```bash
# Update the file path if needed
Rscript 03_create_NACC_model_data.R
```

---

### Stage 2b: Data Exploration (Optional)

#### Script: `04_NACC_data_exploration.R`

**Purpose:** Exploratory data analysis and descriptive statistics.

**Input:**
- `inst/extdata/NACC_mod_dat.csv` (from Stage 2a)

**Output:**
- Figures and tables (graphics device output)
- No data files generated

**Computational Requirements:**
- ✅ **LOW** - Purely exploratory, no modeling

---

### Stage 2c: Model Fitting

#### Script: `05_fit_NACC_model.R`

**Purpose:** Fits the two-part model (GLM + LTRC) on full NACC data.

**Input:**
- `inst/extdata/NACC_mod_dat.csv` (from Stage 2a)

**Output:**
- `inst/extdata/nacc_pt_1_mod.rds`
  - Part 1 model: GLM for dementia probability
  - Formula: `got_dementia ~ t + years_education + is_female + is_married + comorbidity + as.factor(num_e4) + is_race_black + is_race_other + l`

- `inst/extdata/nacc_pt_2_mod.rds`
  - Part 2 model: LTRC survival model
  - Formula: `Surv(mod_y, got_dementia) ~ t + years_education + is_female + is_married + comorbidity + num_e4_1 + num_e4_2 + is_race_black + is_race_other`
  - Fitted on: Subset of never-demented at baseline who developed dementia

**Computational Requirements:**
- 🟡 **MEDIUM** - Model fitting can be slow
- LTRC model fitting is the bottleneck
- Estimated time: 1-5 minutes depending on system

**How to Run (Local or Server):**
```bash
Rscript 05_fit_NACC_model.R
```

---

### Stage 2d: NACC Bootstrap for Variance (Server-Required)

#### Script: `06_NACC_bootstrap_on_server.R`

**Purpose:** Bootstrap resampling to estimate variance of predictions.

**Input:**
- Command-line argument: `FILE_NUM` (integer)
- `inst/extdata/NACC_mod_dat.csv` (from Stage 2a)
- `inst/extdata/nacc_pt_1_mod.rds` (from Stage 2c) - ASSUMED TO EXIST
- Helper function: `get_clean_model()` (must be available in environment)

**Output:**
- `NACC_mod_bs/file_{FILE_NUM}.csv`
  - Bootstrap results with predictions across prediction grid
  - Columns: t, educ, is_female, is_married, has_comorbidity, e4_alleles_1, e4_alleles_2, is_race_black, is_race_other, l, lp_mu, mu, lp_y, expected_y, ghi_estimate
  - 6 bootstrap iterations per FILE_NUM

**Prediction Grid:**
- t: 5-25 (age relative to 65)
- educ: 16
- is_female: 0, 1
- is_married: 0, 1
- has_comorbidity: 0, 1
- e4_alleles_1: 0, 1
- e4_alleles_2: 0, 1
- is_race_black: 0, 1
- is_race_other: 0
- l: 1

**Computational Requirements:**
- ⚠️ **VERY HIGH** - Bootstrap resampling is intensive
- 6 bootstrap iterations per FILE_NUM, each refitting LTRC model
- Must run on server/HPC

**Dependencies:**
- ✅ Must complete Stage 2c first (model files required)
- Custom `get_clean_model()` function must be defined (check for imports/sources)

**How to Run (Server):**
```bash
Rscript 06_NACC_bootstrap_on_server.R 1
```

---

### Stage 2e: Batch Parallel Bootstrap (Server-Required)

#### Script: `07_run_NACC_bootstrap_parallel.sh`

**Purpose:** Wrapper to run bootstrap 50 times in parallel.

**Input:**
- `06_NACC_bootstrap_on_server.R`

**Output:**
- 50 CSV files: `NACC_mod_bs/file_1.csv` through `NACC_mod_bs/file_50.csv`
- 50 log files: `out_1.txt` through `out_50.txt`

**Computational Requirements:**
- ⚠️ **VERY HIGH** - 50 parallel bootstrap jobs
- Total runtime: 1-4 hours wall time

**How to Run (Server):**
```bash
cd inst/paper_code
bash 07_run_NACC_bootstrap_parallel.sh
wait  # Wait for all jobs to complete
```

---

## SECTION 3: ANALYSIS & RESULTS

### Stage 3a: Model Analysis

#### Script: `08_analyze_NACC_model.R`

**Purpose:** Analyzes fitted models and extracts key results.

**Input:**
- `inst/extdata/nacc_pt_1_mod.rds` (from Stage 2c)
- `inst/extdata/nacc_pt_2_mod.rds` (from Stage 2c)

**Output:**
- Summary tables, figures, model diagnostics

**Dependencies:**
- ✅ Requires Stage 2c completion

---

### Stage 3b: Additional NACC Analyses

#### Script: `09_additional_NACC_analyses.R`

**Purpose:** Supplementary analyses and sensitivity checks.

**Input:**
- `inst/extdata/NACC_mod_dat.csv` (from Stage 2a)
- `inst/extdata/nacc_pt_1_mod.rds` (from Stage 2c)
- `inst/extdata/nacc_pt_2_mod.rds` (from Stage 2c)

**Output:**
- Analysis results and figures

---

### Stage 3c: Reviewer Response Code

#### Script: `10_reviewer_response_code.R`

**Purpose:** Additional analyses responding to reviewer feedback.

**Input:**
- Various data files (depends on content)

**Output:**
- Results and figures for reviewer response

---

### Stage 3d: Historical Analysis Function

#### Script: `11_historical_analysis_function.R`

**Purpose:** Helper functions for historical/retrospective analysis.

**Output:**
- Function definitions (used by other scripts)

---

### Stage 3e: Weighted Analysis

#### Script: `12_weighted_analysis.R`

**Purpose:** Weighted sampling analysis.

**Input:**
- Data from earlier stages

**Output:**
- Weighted analysis results

---

### Stage 3f: Weighted Sampling Simulation

#### Script: `13_simulation_for_weighted_sampling.R`

**Purpose:** Simulation for weighted sampling investigation.

**Output:**
- Simulation results

---

### Stage 3g: Model Checking

#### Script: `14_additional_model_checking.R`

**Purpose:** Additional model diagnostics (likely in response to reviewer comments).

**Input:**
- Fitted models from Stage 2c

**Output:**
- Diagnostic plots and statistics

---

## DEPENDENCY GRAPH

```
External Data (NACC database)
        ↓
03_create_NACC_model_data.R
        ↓
NACC_mod_dat.csv
        ├─→ 04_NACC_data_exploration.R
        │
        ├─→ 05_fit_NACC_model.R
        │       ↓
        │   nacc_pt_1_mod.rds
        │   nacc_pt_2_mod.rds
        │       ↓
        │   ├─→ 08_analyze_NACC_model.R
        │   ├─→ 09_additional_NACC_analyses.R
        │   ├─→ 10_reviewer_response_code.R
        │   ├─→ 14_additional_model_checking.R
        │   │
        │   └─→ 06_NACC_bootstrap_on_server.R (× 50 parallel)
        │           ↓
        │       NACC_mod_bs/file_*.csv
        │
        └─→ 12_weighted_analysis.R
        └─→ 13_simulation_for_weighted_sampling.R
        └─→ 11_historical_analysis_function.R

01_code_for_simulation_on_server.R (× 50 parallel via 02_run_simulation_parallel.sh)
        ↓
ghi_sim_res/error_a_*.rds
```

---

## STEP-BY-STEP REPRODUCTION GUIDE

### Prerequisites
- R with packages: dplyr, readr, ltrc, survival, ggplot2, lubridate, stringr, patchwork, fs, data.table
- External NACC database file (must be obtained separately)
- Access to HPC/server for Stages 1b, 2d, 2e

### For Full Reproduction (Including Computationally Intensive Parts)

**On Server:**
1. Update file paths in `02_run_simulation_parallel.sh` and `07_run_NACC_bootstrap_parallel.sh` if needed
2. Create output directories:
   ```bash
   mkdir -p ghi_sim_res NACC_mod_bs
   ```
3. Run simulation batch:
   ```bash
   bash 02_run_simulation_parallel.sh
   wait
   ```
4. Wait for all simulation jobs to complete

5. Run NACC bootstrap batch (after model fitting, see next):
   ```bash
   bash 07_run_NACC_bootstrap_parallel.sh
   wait
   ```

**Locally (or after server jobs complete):**
1. Run data preparation:
   ```bash
   Rscript 03_create_NACC_model_data.R
   ```
2. Optional: Exploratory analysis:
   ```bash
   Rscript 04_NACC_data_exploration.R
   ```
3. Fit models:
   ```bash
   Rscript 05_fit_NACC_model.R
   ```
4. Run analyses (in any order):
   ```bash
   Rscript 08_analyze_NACC_model.R
   Rscript 09_additional_NACC_analyses.R
   Rscript 10_reviewer_response_code.R
   Rscript 12_weighted_analysis.R
   Rscript 13_simulation_for_weighted_sampling.R
   Rscript 14_additional_model_checking.R
   ```

### For Quick Validation (Skip Expensive Computations)

1. Run data preparation: `03_create_NACC_model_data.R`
2. Run model fitting: `05_fit_NACC_model.R`
3. Run quick analyses: `08_analyze_NACC_model.R`, `09_additional_NACC_analyses.R`
4. Skip Stages 1b, 2e (simulations/bootstraps)

---

## OUTPUT LOCATIONS & SIZES

| Stage | Output File(s) | Size (approx.) | Notes |
|-------|----------------|----------------|-------|
| 1b | `ghi_sim_res/error_a_400_*.rds` × 50 | 50-100 MB each | Large R objects, keep on server or compress |
| 2a | `inst/extdata/NACC_mod_dat.csv` | 1-5 MB | Clean analysis-ready data |
| 2c | `inst/extdata/nacc_pt_1_mod.rds` | < 1 MB | Serialized GLM object |
| 2c | `inst/extdata/nacc_pt_2_mod.rds` | 5-20 MB | Serialized LTRC object |
| 2e | `NACC_mod_bs/file_*.csv` × 50 | 10-50 MB total | Bootstrap predictions |
| 3x | Various figures/tables | Variable | PDF, PNG, CSV depending on script |

---

## NOTES FOR REPRODUCIBILITY

1. **Random Seeds:** If reproducibility of exact values is needed, check that all scripts use `set.seed()`. Currently, many scripts may not have fixed seeds.

2. **External Dependencies:**
   - NACC data must be obtained directly from NACC
   - ltrc package is custom/specific to your lab

3. **Computational Time:**
   - Simulations + Bootstrap: 2-4 hours on HPC
   - Model fitting: 5 minutes local
   - Full analysis scripts: < 30 minutes total

4. **Helper Functions:**
   - `get_clean_model()` is used in multiple scripts but its definition location should be verified
   - Ensure it's either sourced or available in the package

5. **Parallelization:**
   - Current bash scripts use `&` for background jobs with `wait`
   - For larger clusters (SLURM, SGE), modify to use proper job submission

---

## DATASET REFERENCE INFORMATION

We used the data lock as of December 2024.

From `03_create_NACC_model_data.R`:
- **Source:** NACC UDS v3 database (Uniform Data Set)
- **File used:** `2025-01-22_UDS_genotype.csv`
- **Sample criteria:**
  - Age ≥ 65 at baseline
  - Dementia-free at first visit
  - Completed follow-up (NACCDIED == 1)
  - Valid APOE genotype (NACCAPOE != 9)
- **Final sample size:** 5,506

