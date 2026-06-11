# =========================================================
# Script 05 - Convergence check (iterations 2 to 5)
# =========================================================
# AEMET annual precipitation, mainland Spain, 2024.
# Coordinates: UTM 30N (EPSG:25830) in kilometres.
#
# This script does NOT produce new final results. The final bandwidths
# are the second-iteration ones from script 04. Here the procedure is
# run for several further iterations for each criterion to check
# convergence: each criterion is re-estimated from its own current
# bandwidth, its own dependence structure is updated, and the next
# bandwidth is computed. Comparing the bandwidths across iterations is
# the evidence that the procedure has (or, for CGCV, has not) stabilised.
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
# 3. Dependence structure for a given trend estimate
# =========================================================
# Pilot binned semivariogram, smoothing bandwidth, bias-corrected
# nonparametric semivariogram, and Shapiro-Botha fit. Each criterion
# uses its own residuals, so each gets its own dependence estimate.

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


# =========================================================
# 4. Criterion settings and one-iteration step
# =========================================================
# Settings per criterion. h.lower_mult is the multiplier of
# bin_data$grid$lag used as lower bound, or NULL to use the npsp default
# (CGCV and CMCV1). These match the second-iteration settings (script
# 04), so the comparison across iterations reflects convergence and not
# a change in the search constraint.

crit_settings <- list(
  CCV   = list(objective = "CV",  ncv = 1, h.lower_mult = 1.5),
  CGCV  = list(objective = "GCV", ncv = 0, h.lower_mult = NULL),
  CMCV1 = list(objective = "CV",  ncv = 2, h.lower_mult = NULL),
  CMCV2 = list(objective = "CV",  ncv = 3, h.lower_mult = 3)
)

# From a current bandwidth matrix, build the trend, take residuals,
# update the dependence structure, and select the next bandwidth.
next_h <- function(h_mat, settings) {
  lp  <- npsp::locpol(x = bin_data, h = h_mat, degree = 1, hat.bin = TRUE)
  res <- as.vector(residuals(lp))
  fit <- fit_dependence(lp, res)
  
  args <- list(
    bin         = bin_data,
    objective   = settings$objective,
    ncv         = settings$ncv,
    degree      = 1,
    cov.dat     = fit,
    DEalgorithm = FALSE,
    warn        = TRUE
  )
  if (!is.null(settings$h.lower_mult)) {
    args$h.lower <- settings$h.lower_mult * bin_data$grid$lag
  }
  do.call(npsp::hcv.data, args)
}

# Summary table for one iteration: h_x, h_y and the objective value.
safe_value <- function(obj) {
  v <- obj$value
  if (is.null(v)) NA_real_ else v
}

summarise_bandwidths <- function(h_list) {
  crit <- names(h_list)
  data.frame(
    criterion = crit,
    h_x   = vapply(crit, function(k) h_list[[k]]$h[1, 1], numeric(1)),
    h_y   = vapply(crit, function(k) h_list[[k]]$h[2, 2], numeric(1)),
    value = vapply(crit, function(k) safe_value(h_list[[k]]),  numeric(1)),
    row.names = NULL
  )
}


# =========================================================
# 5. Run iterations 2 to 5 and print the summary at each step
# =========================================================

h_by_iter <- list()
h_by_iter[["h2"]] <- list(
  CCV   = h2_ccv_from_mcv1,
  CGCV  = h2_cgcv_from_mcv1,
  CMCV1 = h2_cmcv1_from_mcv1,
  CMCV2 = h2_cmcv2_from_mcv1
)

cat("\n--- Iteration 2 (starting point, from script 04) ---\n")
print(summarise_bandwidths(h_by_iter[["h2"]]))

for (step in 2:5) {
  prev    <- h_by_iter[[paste0("h", step - 1)]]
  current <- list()
  for (k in names(crit_settings)) {
    current[[k]] <- next_h(prev[[k]]$h, crit_settings[[k]])
  }
  h_by_iter[[paste0("h", step)]] <- current
  
  cat(sprintf("\n--- Iteration %d ---\n", step))
  print(summarise_bandwidths(current))
}


# =========================================================
# 6. (Optional) Iteration-5 trend surfaces on the fine grid
# =========================================================
# Diagnostic plots of the iteration-5 surfaces. For the stable criteria
# these match the earlier iterations; for CGCV they do not, which is the
# visual counterpart of its non-convergence. Delete this block if only
# the bandwidth tables are needed.

h5 <- h_by_iter[["h16"]]

bin_data_fine <- npsp::binning(
  x = coords_np,
  y = log_z_np,
  nbin = c(90, 90),
  set.NA = FALSE
)

lp_ccv_fine_iter5   <- npsp::locpol(x = bin_data_fine, h = h5$CCV$h,   degree = 1, hat.bin = FALSE)
lp_cgcv_fine_iter5  <- npsp::locpol(x = bin_data_fine, h = h5$CGCV$h,  degree = 1, hat.bin = FALSE)
lp_cmcv1_fine_iter5 <- npsp::locpol(x = bin_data_fine, h = h5$CMCV1$h, degree = 1, hat.bin = FALSE)
lp_cmcv2_fine_iter5 <- npsp::locpol(x = bin_data_fine, h = h5$CMCV2$h, degree = 1, hat.bin = FALSE)

lp_ccv_fine_masked_iter5   <- npsp::mask(lp_ccv_fine_iter5,   window = study_window_sp, set.NA = TRUE)
lp_cgcv_fine_masked_iter5  <- npsp::mask(lp_cgcv_fine_iter5,  window = study_window_sp, set.NA = TRUE)
lp_cmcv1_fine_masked_iter5 <- npsp::mask(lp_cmcv1_fine_iter5, window = study_window_sp, set.NA = TRUE)
lp_cmcv2_fine_masked_iter5 <- npsp::mask(lp_cmcv2_fine_iter5, window = study_window_sp, set.NA = TRUE)

slim_common <- range(c(
  lp_ccv_fine_masked_iter5$est,
  lp_cgcv_fine_masked_iter5$est,
  lp_cmcv1_fine_masked_iter5$est,
  lp_cmcv2_fine_masked_iter5$est
), na.rm = TRUE)

npsp::simage(lp_ccv_fine_masked_iter5,   slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
npsp::simage(lp_cgcv_fine_masked_iter5,  slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
npsp::simage(lp_cmcv1_fine_masked_iter5, slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
npsp::simage(lp_cmcv2_fine_masked_iter5, slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
