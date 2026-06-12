# Macroeconomics – Assignment 2, Exercise 1
## Forecasting the Term Structure of Interest Rates

This repository contains a single, self-contained and fully commented R script
that solves **Exercise 1** of Assignment 2 (Nelson–Siegel yield-curve modelling
and forecasting), following the assignment directives question by question.

> Authors: replace the placeholder names at the top of the script with the
> full name and surname of every group member.

### How to run

```r
# from the repository root, in R or RStudio:
source("R/exercise1_term_structure.R")
```

The script downloads the data directly from **FRED**, prints all requested
numbers to the console, and writes every figure to the `output/` folder.
It requires internet access and the packages `quantmod`, `xts`, `zoo`, `vars`
(installed automatically on first run). The Kalman filter is hand-coded, so no
state-space package is needed.

### What the script does (mapping to the assignment)

| Question | Content |
|---|---|
| **Q1** | Download daily Treasury yields (8 maturities, 3–120 months), convert to monthly averages, report sample, plot the curve at three dates, build & comment empirical level/slope/curvature. |
| **Q2** | Build the Nelson–Siegel loading matrix (λ = 0.0609), estimate the three factors month-by-month by OLS, plot them, compare with the empirical measures, report fitting errors by maturity. |
| **Q3** | Fit a VAR to the factors, produce recursive out-of-sample yield forecasts for h = 1, 6, 12, compare with the random-walk benchmark, report MAE and RMSE by maturity and horizon. |
| **Q4** | One-step state-space dynamic Nelson–Siegel estimated by a hand-coded Kalman filter (diagonal Q and H); plot the filtered factors, forecast recursively, compare against VAR and random walk. |
| **Q5** | Written interpretation (in-code) of short- vs long-horizon accuracy, where the gains concentrate, factor smoothness, and why a low-dimensional factor structure helps. |

### Maturities and data

Eight FRED constant-maturity series are used so that the panel is balanced and
includes the 3-, 24- and 120-month maturities needed by the empirical formulas:

`DGS3MO, DGS6MO, DGS1, DGS2, DGS3, DGS5, DGS7, DGS10`
→ maturities (months) `3, 6, 12, 24, 36, 60, 84, 120`.

### Outputs (`output/`)

- `q1_yield_curves.png`, `q1_empirical_factors.png`
- `q2_ns_loadings.png`, `q2_ns_factors.png`, `q2_factor_vs_empirical.png`, `q2_fitting_errors.png`
- `q4_kalman_factors.png`, `q4_kalman_vs_ols.png`, `q4_model_comparison.png`

Every step, output and graph is interpreted directly inside the R script via
inline comments, as required by the assignment.
