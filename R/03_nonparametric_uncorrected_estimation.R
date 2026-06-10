# =========================================================
# Script 03 - Classical (uncorrected) bandwidth selection
# =========================================================
# AEMET annual precipitation, mainland Spain, 2024.
# Coordinates: UTM 30N (EPSG:25830) in kilometres.
#
# This script performs the local linear trend estimation with the four
# CLASSICAL bandwidth selectors (CV, GCV, MCV1, MCV2), which ignore
# spatial dependence. These provide the uncorrected half of the
# bandwidth comparison table, and MCV1 is used as the pilot bandwidth
# for the dependence-corrected procedure in the later scripts. It also
# builds the mainland Spain masking window used throughout the analysis.
#
# Reads the frozen dataset from script 00. Paths are relative to the
# project root.
# =========================================================


# =========================================================
# 1. Load libraries
# =========================================================

library(mapSpain)
library(sf)
library(dplyr)
library(stringr)
library(sp)
library(npsp)


# =========================================================
# 2. Read frozen dataset and prepare data
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

# Spatial coordinates (km) and response on the log scale
coords_np <- as.matrix(model_data[, c("x", "y")])
log_z_np <- model_data$log_annual_precip


# =========================================================
# 3. Mainland Spain window for plotting/masking
# =========================================================
# The province polygons are restricted to mainland Spain, transformed to
# UTM 30N and divided by 1000, so that the window coordinates match the
# x, y variables (in km) used in the model. This window is reused by the
# later scripts to mask the estimated surfaces.

spain_provinces <- mapSpain::esp_get_prov(year = "2024")

excluded_codes <- c("07", "35", "38", "51", "52")

spain_mainland_sf <- spain_provinces %>%
  mutate(cpro = str_pad(as.character(cpro), width = 2, pad = "0")) %>%
  filter(!cpro %in% excluded_codes) %>%
  st_transform(25830) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

spain_mainland_km <- spain_mainland_sf
st_geometry(spain_mainland_km) <- st_geometry(spain_mainland_km) / 1000
st_crs(spain_mainland_km) <- NA

study_window_sp <- as(spain_mainland_km, "Spatial")


# =========================================================
# 4. Binning for local linear estimation
# =========================================================
# A 25 x 25 grid gives approximately one observation per cell for the
# 690 stations. A finer grid is used later only for display.

nbin_main <- c(25, 25)

bin_data <- npsp::binning(
  x = coords_np,
  y = log_z_np,
  nbin = nbin_main
)


# =========================================================
# 5. Classical bandwidth selection
# =========================================================
# CV:   objective = "CV",  ncv = 1
# GCV:  objective = "GCV", ncv = 0
# MCV1: objective = "CV",  ncv = 2
# MCV2: objective = "CV",  ncv = 3
# The lower bound depends on the binning cell size (bin_data$grid$lag).

bin_data$grid$lag

h_cv <- npsp::hcv.data(
  bin = bin_data,
  objective = "CV",
  h.lower = 1.5 * bin_data$grid$lag,
  ncv = 1,
  degree = 1,
  warn = TRUE,
  DEalgorithm = FALSE
)

h_gcv <- npsp::hcv.data(
  bin = bin_data,
  objective = "GCV",
  ncv = 0,
  degree = 1,
  warn = TRUE,
  DEalgorithm = FALSE
)

h_mcv1 <- npsp::hcv.data(
  bin = bin_data,
  objective = "CV",
  ncv = 2,
  degree = 1,
  warn = TRUE,
  DEalgorithm = FALSE
)

h_mcv2 <- npsp::hcv.data(
  bin = bin_data,
  objective = "CV",
  h.lower = 3 * bin_data$grid$lag,
  ncv = 3,
  degree = 1,
  warn = TRUE,
  DEalgorithm = FALSE
)


# =========================================================
# 6. Store and compare selected bandwidths
# =========================================================
# These are the uncorrected bandwidths (left half of the comparison
# table in the thesis).

bandwidth_results <- data.frame(
  criterion = c("CV", "GCV", "MCV1", "MCV2"),
  h_x = c(
    h_cv$h[1, 1],
    h_gcv$h[1, 1],
    h_mcv1$h[1, 1],
    h_mcv2$h[1, 1]
  ),
  h_y = c(
    h_cv$h[2, 2],
    h_gcv$h[2, 2],
    h_mcv1$h[2, 2],
    h_mcv2$h[2, 2]
  ),
  value = c(
    h_cv$value,
    h_gcv$value,
    h_mcv1$value,
    h_mcv2$value
  )
)

bandwidth_results


# =========================================================
# 7. Local linear trend estimation on the main grid
# =========================================================

lp_cv <- npsp::locpol(x = bin_data, h = h_cv$h, degree = 1, hat.bin = TRUE)
lp_gcv <- npsp::locpol(x = bin_data, h = h_gcv$h, degree = 1, hat.bin = TRUE)
lp_mcv1 <- npsp::locpol(x = bin_data, h = h_mcv1$h, degree = 1, hat.bin = TRUE)
lp_mcv2 <- npsp::locpol(x = bin_data, h = h_mcv2$h, degree = 1, hat.bin = TRUE)


# =========================================================
# 8. Mask estimated surfaces
# =========================================================

lp_cv_masked <- npsp::mask(lp_cv, window = study_window_sp, set.NA = TRUE)
lp_gcv_masked <- npsp::mask(lp_gcv, window = study_window_sp, set.NA = TRUE)
lp_mcv1_masked <- npsp::mask(lp_mcv1, window = study_window_sp, set.NA = TRUE)
lp_mcv2_masked <- npsp::mask(lp_mcv2, window = study_window_sp, set.NA = TRUE)


# =========================================================
# 9. Plot main-grid surfaces
# =========================================================
# Diagnostic plots (not included in the thesis): the thesis shows the
# pilot MCV1 surface and the dependence-corrected surfaces, not the four
# classical ones. These are kept for inspection.

col_map <- npsp::jet.colors(256)

slim_common <- range(log_z_np, na.rm = TRUE)

npsp::simage(lp_cv_masked,   slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
npsp::simage(lp_gcv_masked,  slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
npsp::simage(lp_mcv1_masked, slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
npsp::simage(lp_mcv2_masked, slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)


# =========================================================
# 10. Fine grid for display
# =========================================================
# The bandwidths are kept fixed; only the display grid is refined.

nbin_fine <- c(90, 90)

bin_data_fine <- npsp::binning(
  x = coords_np,
  y = log_z_np,
  nbin = nbin_fine,
  set.NA = FALSE
)


# =========================================================
# 11. Re-estimate surfaces on the fine grid
# =========================================================

lp_cv_fine <- npsp::locpol(x = bin_data_fine, h = h_cv$h, degree = 1, hat.bin = FALSE)
lp_gcv_fine <- npsp::locpol(x = bin_data_fine, h = h_gcv$h, degree = 1, hat.bin = FALSE)
lp_mcv1_fine <- npsp::locpol(x = bin_data_fine, h = h_mcv1$h, degree = 1, hat.bin = FALSE)
lp_mcv2_fine <- npsp::locpol(x = bin_data_fine, h = h_mcv2$h, degree = 1, hat.bin = FALSE)


# =========================================================
# 12. Mask fine-grid surfaces
# =========================================================

lp_cv_fine_masked <- npsp::mask(lp_cv_fine, window = study_window_sp, set.NA = TRUE)
lp_gcv_fine_masked <- npsp::mask(lp_gcv_fine, window = study_window_sp, set.NA = TRUE)
lp_mcv1_fine_masked <- npsp::mask(lp_mcv1_fine, window = study_window_sp, set.NA = TRUE)
lp_mcv2_fine_masked <- npsp::mask(lp_mcv2_fine, window = study_window_sp, set.NA = TRUE)


# =========================================================
# 13. Plot fine-grid surfaces
# =========================================================
# Diagnostic plots (not included in the thesis), kept for inspection.

slim_common <- range(
  lp_cv_fine_masked$est,
  lp_gcv_fine_masked$est,
  lp_mcv1_fine_masked$est,
  lp_mcv2_fine_masked$est,
  na.rm = TRUE
)

npsp::simage(lp_cv_fine_masked,   slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
npsp::simage(lp_gcv_fine_masked,  slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
npsp::simage(lp_mcv1_fine_masked, slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)
npsp::simage(lp_mcv2_fine_masked, slim = slim_common, main = "", col = col_map,
             xlab = "UTM X coordinate (km)", ylab = "UTM Y coordinate (km)", asp = 1)


# =========================================================
# 14. Pilot fit (MCV1) residual diagnostics
# =========================================================
# Diagnostic plots (not included in the thesis): histogram and density
# of the MCV1 residuals, and an observed-versus-fitted plot, to check
# the pilot bandwidth fit before it is used in the corrected procedure.


res_mcv1 <- residuals(lp_mcv1)

dens <- density(res_mcv1)

hist(
  res_mcv1,
  freq = FALSE,
  main = "",
  xlab = "Residuals",
  col = "grey80",
  border = "white",
  ylim = c(0, max(dens$y))
)

lines(dens, lwd = 2)



trend_obs <- lp_mcv1$data$y - res_mcv1

cor(trend_obs, lp_mcv1$data$y, use = "complete.obs")
obs_mcv1 <- lp_mcv1$data$y
fit_mcv1 <- obs_mcv1 - residuals(lp_mcv1)

ok <- complete.cases(fit_mcv1, obs_mcv1)

plot(fit_mcv1[ok], obs_mcv1[ok],
     xlab = "Estimates",
     ylab = "Observations",
     pch = 1,
     cex = 0.7)

abline(0, 1, lty = 2, col = "grey70")

# =========================================================
# 15. Save objects needed by the dependence-corrected scripts
# =========================================================
# MCV1 is the pilot bandwidth for the corrected procedure, and the
# masking window is reused there, so both are frozen for reuse.

saveRDS(h_mcv1, "data/h_mcv1_pilot.rds")
saveRDS(study_window_sp, "data/study_window_sp.rds")
saveRDS(bandwidth_results, "data/bandwidth_results_classical.rds")

