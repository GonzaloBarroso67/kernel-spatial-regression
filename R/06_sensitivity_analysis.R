# =========================================================
# Script 06 - Sensitivity analyses
# =========================================================
# AEMET annual precipitation, mainland Spain, 2024.
# Coordinates: UTM 30N (EPSG:25830) in kilometres.
#
# Reusable function that runs the dependence-corrected pipeline up to
# the second-iteration bandwidths (CCV, CGCV, CMCV1, CMCV2), so that the
# same procedure can be re-run under different conditions. Each test
# below changes one thing (the quality threshold, the binning grid) and
# compares the resulting bandwidths with the base case.
#
# NOTE: the threshold test (Section 4) compares both selected bandwidths
# and estimated surfaces (relative difference in precipitation, %).
# The binning grid test (Section 5) compares selected bandwidths only;
# surfaces are plotted for visual inspection but not numerically compared.
#
# Reads the frozen dataset from script 00. Paths are relative to the
# project root.
# =========================================================


# =========================================================
# 1. Load libraries
# =========================================================

library(sf)
library(dplyr)
library(sp)
library(npsp)


# =========================================================
# 2. Reusable corrected-pipeline function
# =========================================================
# Given a data frame (with columns x, y, log_annual_precip), a selection
# binning grid and a fine binning grid for display, this runs the full
# dependence-corrected procedure from the MCV1 pilot up to the four
# corrected bandwidths, reselecting the pilot inside so that the whole
# pipeline is applied to the given data. The lower-bound rule matches
# script 04: CCV uses 1.5 * lag, CMCV2 uses 3 * lag, and CGCV / CMCV1 use
# the npsp default.
#
# It returns the four bandwidths and the four masked fine-grid surfaces.
# The fine grid is passed in from outside so that two runs share exactly
# the same evaluation points, which is required to subtract the surfaces.

corrected_pipeline <- function(dat, nbin, bin_fine, window) {
  
  coords <- as.matrix(dat[, c("x", "y")])
  logz <- dat$log_annual_precip
  
  bin <- npsp::binning(x = coords, y = logz, nbin = nbin)
  lag <- bin$grid$lag
  
  maxlag <- sqrt(
    diff(range(coords[, 1]))^2 + diff(range(coords[, 2]))^2
  )
  
  # Pilot trend (MCV1)
  h_pilot <- npsp::hcv.data(
    bin = bin, objective = "CV", ncv = 2, degree = 1,
    DEalgorithm = FALSE, warn = FALSE
  )
  lp_pilot <- npsp::locpol(x = bin, h = h_pilot$h, degree = 1, hat.bin = TRUE)
  res <- as.vector(residuals(lp_pilot))
  
  # Dependence structure from the pilot residuals
  sv_bin <- npsp::svar.bin(
    x = coords, y = res,
    maxlag = maxlag / 6, nlags = 50, estimator = "classical"
  )
  g <- npsp::h.cv(bin = sv_bin, loss = "MRSE", degree = 1)$h
  sv_corr <- npsp::np.svariso.corr(
    lp = lp_pilot, maxlag = maxlag / 6, nlags = 50, h = g,
    degree = 1, tol = 0.03, max.iter = 10, plot = FALSE, verbose = FALSE
  )
  fit_sb <- npsp::fitsvar.sb.iso(esv = sv_corr, dk = 0, method = "cressie")
  
  # Four dependence-corrected criteria (same bound rule as script 04)
  h_ccv <- npsp::hcv.data(
    bin = bin, objective = "CV", h.lower = 1.5 * lag, ncv = 1, degree = 1,
    cov.dat = fit_sb, DEalgorithm = FALSE, warn = FALSE
  )
  h_cgcv <- npsp::hcv.data(
    bin = bin, objective = "GCV", ncv = 0, degree = 1,
    cov.dat = fit_sb, DEalgorithm = FALSE, warn = FALSE
  )
  h_cmcv1 <- npsp::hcv.data(
    bin = bin, objective = "CV", ncv = 2, degree = 1,
    cov.dat = fit_sb, DEalgorithm = FALSE, warn = FALSE
  )
  h_cmcv2 <- npsp::hcv.data(
    bin = bin, objective = "CV", h.lower = 3 * lag, ncv = 3, degree = 1,
    cov.dat = fit_sb, DEalgorithm = FALSE, warn = FALSE
  )
  
  bandwidths <- data.frame(
    criterion = c("CCV", "CGCV", "CMCV1", "CMCV2"),
    h_x = c(h_ccv$h[1, 1], h_cgcv$h[1, 1], h_cmcv1$h[1, 1], h_cmcv2$h[1, 1]),
    h_y = c(h_ccv$h[2, 2], h_cgcv$h[2, 2], h_cmcv1$h[2, 2], h_cmcv2$h[2, 2])
  )
  
  # Surfaces on the shared fine grid, masked to mainland Spain
  surf <- function(h) {
    lp <- npsp::locpol(x = bin_fine, h = h, degree = 1, hat.bin = FALSE)
    npsp::mask(lp, window = window, set.NA = TRUE)$est
  }
  
  surfaces <- list(
    CCV = surf(h_ccv$h),
    CGCV = surf(h_cgcv$h),
    CMCV1 = surf(h_cmcv1$h),
    CMCV2 = surf(h_cmcv2$h)
  )
  
  list(bandwidths = bandwidths, surfaces = surfaces)
}


# =========================================================
# 3. Read frozen dataset and build shared objects
# =========================================================

aemet_annual_precip_utm <- readRDS(
  "data/aemet_2024_annual_precip_points_sf.rds"
)

model_data <- aemet_annual_precip_utm %>%
  st_drop_geometry() %>%
  filter(
    !is.na(x),
    !is.na(y),
    !is.na(log_annual_precip)
  )

study_window_sp <- readRDS("data/study_window_sp.rds")

# Shared fine grid, built once on the full domain so that both runs are
# evaluated at exactly the same points (required to subtract surfaces).
coords_full <- as.matrix(model_data[, c("x", "y")])
bin_fine_shared <- npsp::binning(
  x = coords_full,
  y = model_data$log_annual_precip,
  nbin = c(90, 90),
  set.NA = FALSE
)


# =========================================================
# 4. Sensitivity to the quality threshold (>= 300 vs >= 346)
# =========================================================
# The stricter threshold allows at most 20 missing days (366 - 346)
# instead of 66. Both cases use the same 25 x 25 selection grid and the
# same fine grid.

model_data_346 <- model_data %>%
  filter(n_valid_prec >= 346)

cat("Stations with >= 300 valid days:", nrow(model_data), "\n")
cat("Stations with >= 346 valid days:", nrow(model_data_346), "\n")

run_base <- corrected_pipeline(model_data, c(25, 25), bin_fine_shared, study_window_sp)
run_346 <- corrected_pipeline(model_data_346, c(25, 25), bin_fine_shared, study_window_sp)


# ---- 4.1 Bandwidth comparison ----

threshold_bw <- data.frame(
  criterion = run_base$bandwidths$criterion,
  h_x_300 = run_base$bandwidths$h_x,
  h_x_346 = run_346$bandwidths$h_x,
  h_y_300 = run_base$bandwidths$h_y,
  h_y_346 = run_346$bandwidths$h_y
) %>%
  mutate(
    rel_change_x = 100 * (h_x_346 - h_x_300) / h_x_300,
    rel_change_y = 100 * (h_y_346 - h_y_300) / h_y_300
  )

threshold_bw


# ---- 4.2 Surface comparison ----
# Relative difference in PRECIPITATION (mm) at each grid point, in per
# cent: 100 * |exp(est_346) - exp(est_300)| / exp(est_300). The 95th
# percentile summarises the typical change across the region; the
# maximum flags localised changes (relevant for CGCV).

surface_diff <- function(est_300, est_346) {
  p300 <- exp(est_300)
  p346 <- exp(est_346)
  rel <- 100 * abs(p346 - p300) / p300
  rel <- rel[is.finite(rel)]
  c(
    median = median(rel),
    q95 = unname(quantile(rel, 0.95)),
    max = max(rel)
  )
}

threshold_surface <- t(sapply(
  c("CCV", "CGCV", "CMCV1", "CMCV2"),
  function(cr) surface_diff(run_base$surfaces[[cr]], run_346$surfaces[[cr]])
))

threshold_surface


# =========================================================
# 5. Sensitivity to the selection grid (25×25 vs 30×30 and 20×20)
# =========================================================
# Only the coarse selection grid changes; the full dataset and the fine
# display grid are the same. Both alternatives (30×30 and 20×20) are
# compared against the base case (25×25) in terms of selected bandwidths
# only; surfaces are plotted for visual inspection but not numerically
# compared.

run_30 <- corrected_pipeline(model_data, c(30, 30), bin_fine_shared, study_window_sp)

grid_bw <- data.frame(
  criterion = run_base$bandwidths$criterion,
  h_x_25 = run_base$bandwidths$h_x,
  h_x_30 = run_30$bandwidths$h_x,
  h_y_25 = run_base$bandwidths$h_y,
  h_y_30 = run_30$bandwidths$h_y
) %>%
  mutate(
    rel_change_x = 100 * (h_x_30 - h_x_25) / h_x_25,
    rel_change_y = 100 * (h_y_30 - h_y_25) / h_y_25
  )

grid_bw



col_map <- npsp::jet.colors(256)
criterios <- c("CCV", "CGCV", "CMCV1", "CMCV2")

# Reconstruct the 8 surfaces (4 criteria x 2 grids) from their bandwidths
make_surface <- function(bw_row) {
  h <- diag(c(bw_row$h_x, bw_row$h_y))
  npsp::mask(
    npsp::locpol(x = bin_fine_shared, h = h, degree = 1, hat.bin = FALSE),
    window = study_window_sp, set.NA = TRUE
  )
}

surf_25 <- lapply(criterios, function(cr)
  make_surface(run_base$bandwidths[run_base$bandwidths$criterion == cr, ]))
surf_30 <- lapply(criterios, function(cr)
  make_surface(run_30$bandwidths[run_30$bandwidths$criterion == cr, ]))
names(surf_25) <- criterios
names(surf_30) <- criterios

# Common color scale for the 8 surfaces
slim_common <- range(
  c(sapply(surf_25, function(s) s$est), sapply(surf_30, function(s) s$est)),
  na.rm = TRUE
)

# Plot: for each criterion, first 25x25, then 30x30
for (cr in criterios) {
  npsp::simage(surf_25[[cr]], slim = slim_common, main = paste(cr, "- 25x25"),
               col = col_map, xlab = "UTM X coordinate (km)",
               ylab = "UTM Y coordinate (km)", asp = 1)
  npsp::simage(surf_30[[cr]], slim = slim_common, main = paste(cr, "- 30x30"),
               col = col_map, xlab = "UTM X coordinate (km)",
               ylab = "UTM Y coordinate (km)", asp = 1)
}



run_20 <- corrected_pipeline(model_data, c(20, 20), bin_fine_shared, study_window_sp)

grid_bw_20 <- data.frame(
  criterion = run_base$bandwidths$criterion,
  h_x_25 = run_base$bandwidths$h_x,
  h_x_20 = run_20$bandwidths$h_x,
  h_y_25 = run_base$bandwidths$h_y,
  h_y_20 = run_20$bandwidths$h_y
)
grid_bw_20

col_map <- npsp::jet.colors(256)
criterios <- c("CCV", "CGCV", "CMCV1", "CMCV2")

make_surface <- function(bw_row) {
  h <- diag(c(bw_row$h_x, bw_row$h_y))
  npsp::mask(
    npsp::locpol(x = bin_fine_shared, h = h, degree = 1, hat.bin = FALSE),
    window = study_window_sp, set.NA = TRUE
  )
}

surf_25 <- lapply(criterios, function(cr)
  make_surface(run_base$bandwidths[run_base$bandwidths$criterion == cr, ]))
surf_20 <- lapply(criterios, function(cr)
  make_surface(run_20$bandwidths[run_20$bandwidths$criterion == cr, ]))
names(surf_25) <- criterios
names(surf_20) <- criterios

# Common color scale for the 8 surfaces
slim_common <- range(
  c(sapply(surf_25, function(s) s$est), sapply(surf_20, function(s) s$est)),
  na.rm = TRUE
)

# Plot: for each criterion, first 25x25, then 20x20
for (cr in criterios) {
  npsp::simage(surf_25[[cr]], slim = slim_common, main = paste(cr, "- 25x25"),
               col = col_map, xlab = "UTM X coordinate (km)",
               ylab = "UTM Y coordinate (km)", asp = 1)
  npsp::simage(surf_20[[cr]], slim = slim_common, main = paste(cr, "- 20x20"),
               col = col_map, xlab = "UTM X coordinate (km)",
               ylab = "UTM Y coordinate (km)", asp = 1)
}
