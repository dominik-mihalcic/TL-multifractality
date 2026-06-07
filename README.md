# Scale invariance and multifractality in Taylor's law

This repository contains R scripts used to generate figures presented in the paper *Scale invariance and multifractality in Taylor's law*. The code implements simulation and statistical analysis of log-normal multiplicative cascades, including moment scaling analysis and estimation of the scaling exponent $b(q)$.

All scripts are independent entry points and share common functions and constants defined in `helper.R`.

---

## Structure

- `helper.R` - shared functions, constants and plotting utilities
- `cascade_plots.R` - simulation and visualization of 1D and 2D cascades
- `b_plots.R` - estimation of scaling exponents b(q)
- `moment_plots.R` - empirical moment scaling analysis
- `TL_data.xlsx` - empirical datasets used in the analysis

---

## Dependency

All scripts depend on:

- `helper.R` (must be in working directory)

R packages:
- tidyverse
- plot3D
- readxl

---

## How to run

Open each script in R (or RStudio) and execute line by line.

Each script is self-contained and can be run independently, provided that `helper.R` is in the working directory.

Recommended order (optional):

1. `moment_plots.R`
2. `cascade_plots.R`
3. `b_plots.R`
