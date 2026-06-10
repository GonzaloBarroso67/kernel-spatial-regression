# =========================================================
# Script 04 - Dependence-corrected bandwidth selection
# First iteration from the MCV1 pilot bandwidth
# =========================================================
# AEMET annual precipitation, mainland Spain, 2024.
# Coordinates: UTM 30N (EPSG:25830) in kilometres.
#
# Starting from the MCV1 pilot bandwidth (script 03), this script
# estimates an initial trend, computes the bias-corrected residual
# variogram, fits a Shapiro-Botha model, and uses that single dependence
# estimate to compute the four dependence-corrected bandwidth selectors:
# CCV, CGCV, CMCV1 and CMCV2. The four criteria share the same fitted
# covariance structure so that they are directly comparable.
#
# Reads the frozen dataset from script 00 and the pilot bandwidth and
# masking window from script 03. Paths are relative to the project root.
# =========================================================


# =========================================================
# 1. Load libraries
# =========================================================

library(sf)
library(dplyr)
library(sp)
library(npsp)


# =========================================================
# 2. Read frozen objects and rebuild working objects
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

coords_np <- as.matrix(model_data[, c("x", "y")])
log_z_np <- model_data$log_annual_precip

h_mcv1 <- readRDS("data/h_mcv1_pilot.rds")
study_window_sp <- readRDS("data/study_window_sp.rds")

# Binning grid (same 25 x 25 as in script 03) and maximum lag
bin_data <- npsp::binning(
  x = coords_np,
  y = log_z_np,
  nbin = c(25, 25)
)

maxlag_np <- sqrt(
  diff(range(coords_np[, 1]))^2 + diff(range(coords_np[, 2]))^2
)

col_map <- npsp::jet.colors(256)


# =========================================================
# 3. Pilot trend and residuals (MCV1)
# =========================================================

lp_mcv1_iter1 <- npsp::locpol(
  x = bin_data,
  h = h_mcv1$h,
  degree = 1,
  hat.bin = TRUE
)

model_data$res_mcv1_iter1 <- as.vector(residuals(lp_mcv1_iter1))
model_data$mu_hat_mcv1_iter1 <- log_z_np - model_data$res_mcv1_iter1


# =========================================================
# 4. Pilot binned semivariogram and pair-count check
# =========================================================
# The maximum lag is about one sixth of the maximum distance between
# stations. The pair-count check confirms that the least populated lag
# class stays above the recommended minimum of 30 pairs.

sv_bin_mcv1_iter1 <- npsp::svar.bin(
  x = coords_np,
  y = model_data$res_mcv1_iter1,
  maxlag = maxlag_np / 6,
  nlags = 50,
  estimator = "classical"
)

sv_check <- data.frame(
  distance = as.numeric(npsp::coords(sv_bin_mcv1_iter1)),
  npairs = as.numeric(sv_bin_mcv1_iter1$binw)
)

summary(sv_check$npairs)


# =========================================================
# 5. Bandwidth for variogram smoothing
# =========================================================

h_sv_mcv1_iter1 <- npsp::h.cv(
  bin = sv_bin_mcv1_iter1,
  loss = "MRSE",
  degree = 1
)

g_mcv1_iter1 <- h_sv_mcv1_iter1$h
g_mcv1_iter1


# =========================================================
# 6. Uncorrected and bias-corrected nonparametric semivariograms
# =========================================================
# Both estimates are computed here; they are only plotted together,
# after fitting, in the next section.

sv_ll_mcv1_iter1 <- npsp::np.svariso(
  x = coords_np,
  y = model_data$res_mcv1_iter1,
  maxlag = maxlag_np / 6,
  nlags = 50,
  h = g_mcv1_iter1,
  degree = 1,
  plot = FALSE
)

sv_corr_mcv1_iter1 <- npsp::np.svariso.corr(
  lp = lp_mcv1_iter1,
  maxlag = maxlag_np / 6,
  nlags = 50,
  h = g_mcv1_iter1,
  degree = 1,
  tol = 0.03,
  max.iter = 10,
  plot = FALSE,
  verbose = FALSE
)


# =========================================================
# 7. Shapiro-Botha fits and bias-correction figure
# =========================================================
# The corrected and biased estimates are fitted with the Shapiro-Botha
# model and shown together. This is the bias-correction figure in the
# thesis: the solid line is the corrected estimate, the dashed line the
# original biased one.

fit_sb_biased <- npsp::fitsvar.sb.iso(
  esv = sv_ll_mcv1_iter1,
  dk = 0,
  method = "cressie"
)

fit_sb_mcv1_iter1 <- npsp::fitsvar.sb.iso(
  esv = sv_corr_mcv1_iter1,
  dk = 0,
  method = "cressie"
)

ylim_max <- 1.1 * max(
  max(fit_sb_biased$fit$sv, na.rm = TRUE),
  max(fit_sb_mcv1_iter1$fit$sv, na.rm = TRUE)
)

# Corrected Shapiro-Botha (solid)
plot(
  fit_sb_mcv1_iter1,
  ylim = c(0, ylim_max),
  xlab = "distance",
  ylab = "semivariance",
  legend = FALSE,
  lwd = c(0, 2),  
  pch = NA
)

# Biased Shapiro-Botha on top (dashed, triangles)
plot(
  fit_sb_biased,
  add = TRUE,
  legend = FALSE,
  lty = c(2, 2),
  lwd = c(1, 2),
  pch = 2
)

legend(
  "bottomright",
  legend = c("Corrected S-B", "Biased S-B"),
  lty = c(1, 2),
  lwd = c(2, 2),
  bty = "n",
  cex = 0.85
)


# =========================================================
# 8. Dependence-corrected bandwidths
# =========================================================
# All four criteria use the same fitted dependence structure
# (fit_sb_mcv1_iter1), so the differences between them come from the
# criterion, not from a different assumed covariance. The lower bounds
# match those used for the classical criteria in script 03.

h2_ccv_from_mcv1 <- npsp::hcv.data(
  bin = bin_data,
  objective = "CV",
  h.lower = 1.5 * bin_data$grid$lag,
  ncv = 1,
  degree = 1,
  cov.dat = fit_sb_mcv1_iter1,
  DEalgorithm = FALSE,
  warn = TRUE
)

h2_cgcv_from_mcv1 <- npsp::hcv.data(
  bin = bin_data,
  objective = "GCV",
  ncv = 0,
  degree = 1,
  cov.dat = fit_sb_mcv1_iter1,
  DEalgorithm = FALSE,
  warn = TRUE
)

h2_cmcv1_from_mcv1 <- npsp::hcv.data(
  bin = bin_data,
  objective = "CV",
  ncv = 2,
  degree = 1,
  cov.dat = fit_sb_mcv1_iter1,
  DEalgorithm = FALSE,
  warn = TRUE
)

h2_cmcv2_from_mcv1 <- npsp::hcv.data(
  bin = bin_data,
  objective = "CV",
  h.lower = 3 * bin_data$grid$lag,
  ncv = 3,
  degree = 1,
  cov.dat = fit_sb_mcv1_iter1,
  DEalgorithm = FALSE,
  warn = TRUE
)

# =========================================================
# 9. Store and compare selected bandwidths
# =========================================================
# These are the corrected bandwidths (right half of the comparison
# table in the thesis).

bandwidth_results_corrected <- data.frame(
  criterion = c("CCV", "CGCV", "CMCV1", "CMCV2"),
  h_x = c(
    h2_ccv_from_mcv1$h[1, 1],
    h2_cgcv_from_mcv1$h[1, 1],
    h2_cmcv1_from_mcv1$h[1, 1],
    h2_cmcv2_from_mcv1$h[1, 1]
  ),
  h_y = c(
    h2_ccv_from_mcv1$h[2, 2],
    h2_cgcv_from_mcv1$h[2, 2],
    h2_cmcv1_from_mcv1$h[2, 2],
    h2_cmcv2_from_mcv1$h[2, 2]
  ),
  value = c(
    h2_ccv_from_mcv1$value,
    h2_cgcv_from_mcv1$value,
    h2_cmcv1_from_mcv1$value,
    h2_cmcv2_from_mcv1$value
  )
)

bandwidth_results_corrected

# =========================================================
# 10. Corrected trend surfaces on the fine grid
# =========================================================
# These are the four dependence-corrected surfaces shown in the thesis.

bin_data_fine <- npsp::binning(
  x = coords_np,
  y = log_z_np,
  nbin = c(90, 90),
  set.NA = FALSE
)

lp_ccv_fine <- npsp::locpol(x = bin_data_fine, h = h2_ccv_from_mcv1$h, degree = 1, hat.bin = FALSE)
lp_cgcv_fine <- npsp::locpol(x = bin_data_fine, h = h2_cgcv_from_mcv1$h, degree = 1, hat.bin = FALSE)
lp_cmcv1_fine <- npsp::locpol(x = bin_data_fine, h = h2_cmcv1_from_mcv1$h, degree = 1, hat.bin = FALSE)
lp_cmcv2_fine <- npsp::locpol(x = bin_data_fine, h = h2_cmcv2_from_mcv1$h, degree = 1, hat.bin = FALSE)

lp_ccv_fine_masked <- npsp::mask(lp_ccv_fine, window = study_window_sp, set.NA = TRUE)
lp_cgcv_fine_masked <- npsp::mask(lp_cgcv_fine, window = study_window_sp, set.NA = TRUE)
lp_cmcv1_fine_masked <- npsp::mask(lp_cmcv1_fine, window = study_window_sp, set.NA = TRUE)
lp_cmcv2_fine_masked <- npsp::mask(lp_cmcv2_fine, window = study_window_sp, set.NA = TRUE)

# Common colour scale across the four corrected surfaces
slim_common <- range(c(
  lp_ccv_fine_masked$est,
  lp_cgcv_fine_masked$est,
  lp_cmcv1_fine_masked$est,
  lp_cmcv2_fine_masked$est
), na.rm = TRUE)

npsp::simage(lp_ccv_fine_masked,   slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
npsp::simage(lp_cgcv_fine_masked,  slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
npsp::simage(lp_cmcv1_fine_masked, slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
npsp::simage(lp_cmcv2_fine_masked, slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)


# =========================================================
# 11. Save objects needed by the next iteration (script 05)
# =========================================================
# The second-iteration script starts from these corrected bandwidths
# and the fitted dependence structure.

h2_iter2 <- list(
  ccv = h2_ccv_from_mcv1,
  cgcv = h2_cgcv_from_mcv1,
  cmcv1 = h2_cmcv1_from_mcv1,
  cmcv2 = h2_cmcv2_from_mcv1
)

saveRDS(h2_iter2, "data/h2_corrected_iter2.rds")
saveRDS(g_mcv1_iter1, "data/g_mcv1_iter1.rds")
saveRDS(fit_sb_mcv1_iter1, "data/fit_sb_mcv1_iter1.rds")




