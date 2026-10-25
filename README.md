# Scale invariance and multifractality in Taylor's law

This repository contains the R code and data for the empirical analysis in the paper *Scale invariance and multifractality in Taylor's law* by D. Grahovac, N. N. Leonenko and D. Mihalčić.

The code estimates the scaling exponents $b(q)$ of the empirical moments, fits the log-normal and log-Poisson cascades, and compares them with the estimates using bootstrap intervals and a simultaneous sup-$t$ band. It is used to reproduce Tables 1–3 and Figures 4 and 5 of the paper.

---

## Files

- `TL_analysis.R` – main script which analyses all datasets and produces the tables and figures
- `TL_helpers.R` – functions used by the main script (empirical moments, estimation of $b(q)$, model fitting, bootstrap and the sup-$t$ band)
- `TL_data.xlsx` – the datasets, one sheet per dataset

---

## Data

Each sheet of `TL_data.xlsx` has two columns. The column `groups` identifies the group, and the column `counts` contains the observations in that group.

| Sheet | Dataset | Source |
|---|---|---|
| `Beall` | Colorado potato beetle | Beall (1938), *Biometrika* 30, 422–439 |
| `Kleczkowski` | Tobacco mosaic virus | Kleczkowski (1949), *Ann. Appl. Biol.* 36, 139–152, grouped as in Bliss and Owen (1958), *Biometrika* 45, 37–58 |
| `Milne` | Sheep tick | Milne (1943), *Ann. Appl. Biol.* 30, 240–250 |
| `Trigal and Ruete` | Invasive flagellate | Trigal and Ruete (2016), *J. Ecol.* 104, 947–956 |
| `Cohen I`, `Cohen IV`, `Cohen V`, `Cohen VII` | Hokkaido vole, subsets I, IV, V and VII | Saitoh and Cohen (2024), *Popul. Ecol.* 66, 125–142 |

The construction of the groups for each dataset is described in Section 4(a) of the paper.

---

## How to run

Requires the R packages `tidyverse` and `readxl`. Place the three files in the working directory and run `TL_analysis.R`, which prints Tables 1–3 and draws Figures 4 and 5.
