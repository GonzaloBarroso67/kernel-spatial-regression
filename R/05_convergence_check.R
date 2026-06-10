# =========================================================
# Script 05 - Convergence check (third iteration)
# =========================================================
# AEMET annual precipitation, mainland Spain, 2024.
# Coordinates: UTM 30N (EPSG:25830) in kilometres.
#
# This script does NOT produce new final results. The final bandwidths
# are the second-iteration ones from script 04. Here one further
# iteration is run for each criterion to check convergence: each
# criterion is re-estimated from its own second-iteration bandwidth,
# its own dependence structure is updated, and the third-iteration
# bandwidth is computed. The comparison of the second- and third-
# iteration bandwidths is the evidence that the procedure has
# (or, for CGCV, has not) stabilised.
#
# Reads the frozen dataset from script 00, the masking window from
# script 03, and the second-iteration bandwidths from script 04. Paths
# are relative to the project root.
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

study_window_sp <- readRDS("data/study_window_sp.rds")

# Second-iteration bandwidths from script 04
h2_iter2 <- readRDS("data/h2_corrected_iter2.rds")
h2_ccv_from_mcv1 <- h2_iter2$ccv
h2_cgcv_from_mcv1 <- h2_iter2$cgcv
h2_cmcv1_from_mcv1 <- h2_iter2$cmcv1
h2_cmcv2_from_mcv1 <- h2_iter2$cmcv2

# Binning grid (same 25 x 25) and maximum lag
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
# 3. Trend estimators from the second-iteration bandwidths
# =========================================================

lp_ccv_iter2 <- npsp::locpol(x = bin_data, h = h2_ccv_from_mcv1$h, degree = 1, hat.bin = TRUE)
lp_cgcv_iter2 <- npsp::locpol(x = bin_data, h = h2_cgcv_from_mcv1$h, degree = 1, hat.bin = TRUE)
lp_cmcv1_iter2 <- npsp::locpol(x = bin_data, h = h2_cmcv1_from_mcv1$h, degree = 1, hat.bin = TRUE)
lp_cmcv2_iter2 <- npsp::locpol(x = bin_data, h = h2_cmcv2_from_mcv1$h, degree = 1, hat.bin = TRUE)

model_data$res_ccv_iter2 <- as.vector(residuals(lp_ccv_iter2))
model_data$res_cgcv_iter2 <- as.vector(residuals(lp_cgcv_iter2))
model_data$res_cmcv1_iter2 <- as.vector(residuals(lp_cmcv1_iter2))
model_data$res_cmcv2_iter2 <- as.vector(residuals(lp_cmcv2_iter2))


# =========================================================
# 4. Updated dependence structure for each criterion
# =========================================================
# For each criterion: pilot binned semivariogram, smoothing bandwidth,
# bias-corrected nonparametric semivariogram, and Shapiro-Botha fit.
# Each criterion uses its own residuals, so each gets its own dependence
# estimate.

fit_dependence <- function(lp_obj, res_vec) {
  sv_bin <- npsp::svar.bin(
    x = coords_np, y = res_vec,
    maxlag = maxlag_np / 6, nlags = 50, estimator = "classical"
  )
  g <- npsp::h.cv(bin = sv_bin, loss = "MRSE", degree = 1)$h
  sv_corr <- npsp::np.svariso.corr(
    lp = lp_obj, maxlag = maxlag_np / 6, nlags = 50, h = g,
    degree = 1, tol = 0.03, max.iter = 10, plot = FALSE, verbose = FALSE
  )
  npsp::fitsvar.sb.iso(esv = sv_corr, dk = 0, method = "cressie")
}

fit_sb_ccv_iter2 <- fit_dependence(lp_ccv_iter2, model_data$res_ccv_iter2)
fit_sb_cgcv_iter2 <- fit_dependence(lp_cgcv_iter2, model_data$res_cgcv_iter2)
fit_sb_cmcv1_iter2 <- fit_dependence(lp_cmcv1_iter2, model_data$res_cmcv1_iter2)
fit_sb_cmcv2_iter2 <- fit_dependence(lp_cmcv2_iter2, model_data$res_cmcv2_iter2)


# =========================================================
# 5. Third-iteration bandwidths
# =========================================================
# The lower bounds match those used in the second iteration (script 04),
# so that the h2-to-h3 comparison reflects convergence and not a change
# in the search constraint.

h3_ccv_from_mcv1 <- npsp::hcv.data(
  bin = bin_data,
  objective = "CV",
  h.lower = 1.5 * bin_data$grid$lag,
  ncv = 1,
  degree = 1,
  cov.dat = fit_sb_ccv_iter2,
  DEalgorithm = FALSE,
  warn = TRUE
)

h3_cgcv_from_mcv1 <- npsp::hcv.data(
  bin = bin_data,
  objective = "GCV",
  ncv = 0,
  degree = 1,
  cov.dat = fit_sb_cgcv_iter2,
  DEalgorithm = FALSE,
  warn = TRUE
)

h3_cmcv1_from_mcv1 <- npsp::hcv.data(
  bin = bin_data,
  objective = "CV",
  ncv = 2,
  degree = 1,
  cov.dat = fit_sb_cmcv1_iter2,
  DEalgorithm = FALSE,
  warn = TRUE
)

h3_cmcv2_from_mcv1 <- npsp::hcv.data(
  bin = bin_data,
  objective = "CV",
  h.lower = 3 * bin_data$grid$lag,
  ncv = 3,
  degree = 1,
  cov.dat = fit_sb_cmcv2_iter2,
  DEalgorithm = FALSE,
  warn = TRUE
)


# =========================================================
# 6. Compare second- and third-iteration bandwidths
# =========================================================
# This is the central output of the script: small relative changes
# indicate that the procedure has converged.

bandwidth_comparison_h2_h3 <- data.frame(
  criterion = c("CCV", "CGCV", "CMCV1", "CMCV2"),
  h2_x = c(
    h2_ccv_from_mcv1$h[1, 1], h2_cgcv_from_mcv1$h[1, 1],
    h2_cmcv1_from_mcv1$h[1, 1], h2_cmcv2_from_mcv1$h[1, 1]
  ),
  h2_y = c(
    h2_ccv_from_mcv1$h[2, 2], h2_cgcv_from_mcv1$h[2, 2],
    h2_cmcv1_from_mcv1$h[2, 2], h2_cmcv2_from_mcv1$h[2, 2]
  ),
  h3_x = c(
    h3_ccv_from_mcv1$h[1, 1], h3_cgcv_from_mcv1$h[1, 1],
    h3_cmcv1_from_mcv1$h[1, 1], h3_cmcv2_from_mcv1$h[1, 1]
  ),
  h3_y = c(
    h3_ccv_from_mcv1$h[2, 2], h3_cgcv_from_mcv1$h[2, 2],
    h3_cmcv1_from_mcv1$h[2, 2], h3_cmcv2_from_mcv1$h[2, 2]
  )
) %>%
  mutate(
    relative_change_x = 100 * (h3_x - h2_x) / h2_x,
    relative_change_y = 100 * (h3_y - h2_y) / h2_y
  )

bandwidth_comparison_h2_h3


# =========================================================
# 7. Third-iteration trend surfaces on the fine grid
# =========================================================
# Diagnostic plots: kept so that anyone running the code can check that
# the third-iteration surfaces are visually indistinguishable from the
# second-iteration ones.

bin_data_fine <- npsp::binning(
  x = coords_np,
  y = log_z_np,
  nbin = c(90, 90),
  set.NA = FALSE
)

lp_ccv_fine_iter3 <- npsp::locpol(x = bin_data_fine, h = h3_ccv_from_mcv1$h, degree = 1, hat.bin = FALSE)
lp_cgcv_fine_iter3 <- npsp::locpol(x = bin_data_fine, h = h3_cgcv_from_mcv1$h, degree = 1, hat.bin = FALSE)
lp_cmcv1_fine_iter3 <- npsp::locpol(x = bin_data_fine, h = h3_cmcv1_from_mcv1$h, degree = 1, hat.bin = FALSE)
lp_cmcv2_fine_iter3 <- npsp::locpol(x = bin_data_fine, h = h3_cmcv2_from_mcv1$h, degree = 1, hat.bin = FALSE)

lp_ccv_fine_masked_iter3 <- npsp::mask(lp_ccv_fine_iter3, window = study_window_sp, set.NA = TRUE)
lp_cgcv_fine_masked_iter3 <- npsp::mask(lp_cgcv_fine_iter3, window = study_window_sp, set.NA = TRUE)
lp_cmcv1_fine_masked_iter3 <- npsp::mask(lp_cmcv1_fine_iter3, window = study_window_sp, set.NA = TRUE)
lp_cmcv2_fine_masked_iter3 <- npsp::mask(lp_cmcv2_fine_iter3, window = study_window_sp, set.NA = TRUE)

slim_common <- range(c(
  lp_ccv_fine_masked_iter3$est,
  lp_cgcv_fine_masked_iter3$est,
  lp_cmcv1_fine_masked_iter3$est,
  lp_cmcv2_fine_masked_iter3$est
), na.rm = TRUE)

npsp::simage(lp_ccv_fine_masked_iter3,   slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
npsp::simage(lp_cgcv_fine_masked_iter3,  slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
npsp::simage(lp_cmcv1_fine_masked_iter3, slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
npsp::simage(lp_cmcv2_fine_masked_iter3, slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)

