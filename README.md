# Kernel-based regression estimation for spatial data

R code for the application developed in my Bachelor's Thesis (UC3M).
The project applies nonparametric kernel-based regression with
bias-corrected variogram estimation and dependence-corrected bandwidth
selection to log annual accumulated precipitation from AEMET
meteorological stations in mainland Spain (2024).

## Data

The data comes from the AEMET OpenData API (2024 daily precipitation).
The `data/` folder contains the preprocessed inputs for the analysis:

- `aemet_2024_annual_precip_points_sf.rds` — annual accumulated
  precipitation per station (690 stations).
- `aemet_2024_missing_by_station_month.rds` — number of missing days
  per station and month, used in the exploratory and sensitivity analyses.

The analysis reads directly from these files, so the download step is
not required to reproduce the results.

## Requirements

- R (version 4.5.2)
- Packages: `npsp`, `sf`, `sp`, `dplyr`, `gstat`, `sm`, `ggplot2`,
  `stringr`, `mapSpain`, `e1071`

Install them with:

```r
install.packages(c("npsp", "sf", "sp", "dplyr", "gstat", "sm",
                   "ggplot2", "stringr", "mapSpain", "e1071"))
```

On Windows, installing `npsp` may require Rtools.

## How to run

Script `00` documents how the data was obtained from the AEMET API and
is not needed to reproduce the analysis (the data is already provided in
`data/`). The analysis starts from script `01`.

- `00_download_AEMET_data.R` — Downloads the 2024 daily precipitation
  data from the AEMET OpenData API (documentary; requires a personal API key).
- `01_exploratory_analysis.R` — Describes the dataset and explores the
  response on the original and log scales (spatial and univariate plots,
  summaries by province, nearest-neighbour distances, and a missing-day
  diagnostic supporting the rescaling assumption).
- `02_parametric_preview_and_diagnostics.R` — Fits an exploratory linear
  trend on the log scale and diagnoses the spatial dependence in the
  residuals (empirical and directional variograms, independence test),
  motivating the dependence-corrected approach.
- `03_nonparametric_uncorrected_estimation.R` — Local linear trend
  estimation with the classical bandwidth selectors (CV, GCV, MCV1, MCV2).
  Builds the mainland Spain masking window and saves the MCV1 pilot
  bandwidth reused by the corrected procedure.
- `04_nonparametric_corrected_estimation.R` — First iteration of the
  dependence-corrected procedure: starting from the MCV1 pilot, estimates
  the trend, computes the bias-corrected residual variogram, fits a
  Shapiro-Botha model, and uses that single dependence structure for the
  four corrected bandwidth selectors (CCV, CGCV, CMCV1, CMCV2).
- `05_convergence_check.R` — Runs one further iteration of the corrected
  procedure from the bandwidths saved in script 04 and checks whether the
  estimates have converged.
- `06_sensitivity_analysis.R` — Sensitivity analysis under two checks: a
  stricter missing-day filter and changes to the binning grid (30×30 and 20×20).


  ## Academic context

This project was developed as part of the Bachelor's Degree in
Mathematics and Computing at Universidad Carlos III de Madrid.

- **Author:** Gonzalo Barroso
- **Thesis title:** Kernel-based regression estimation for spatial data
