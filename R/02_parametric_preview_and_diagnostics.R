# =========================================================
# Script 02 - Parametric baseline and dependence diagnostics
# =========================================================
# AEMET annual precipitation, mainland Spain, 2024.
# Coordinates: UTM 30N (EPSG:25830) in kilometres.
#
# A simple linear trend on the log scale is fitted only as an
# exploratory step, to remove the large-scale mean before inspecting
# the spatial dependence in the residuals through the empirical
# variogram, directional variograms and an independence test. The
# linear trend is not the model used later; it only motivates the
# nonparametric, dependence-corrected approach.
#
# Reads the frozen dataset from script 00. Paths are relative to the
# project root.
# =========================================================


# =========================================================
# 1. Load libraries
# =========================================================

library(dplyr)
library(sf)
library(sp)
library(gstat)
library(npsp)
library(sm)


# =========================================================
# 2. Read frozen dataset and prepare data
# =========================================================

aemet_annual_precip_utm <- readRDS(
  "data/aemet_2024_annual_precip_points_sf.rds"
)

model_data <- st_drop_geometry(aemet_annual_precip_utm)

coords <- as.matrix(model_data[, c("x", "y")])


# =========================================================
# 3. Exploratory linear trend on the log scale
# =========================================================
# Fitted only to remove the dominant large-scale mean before looking at
# the dependence in the residuals.

modelo_lin_log <- lm(
  log_annual_precip ~ x + y,
  data = model_data,
  na.action = na.exclude
)

summary(modelo_lin_log)

model_data$mu_hat_lin_log <- fitted(modelo_lin_log)
model_data$res_ols_log <- residuals(modelo_lin_log)

# Diagnostic plots (not included in the thesis): the fitted linear trend
# surface and the spatial pattern of the residuals. The residual map
# shows that spatial structure remains after removing the linear trend,
# which is what the variogram below quantifies. Residuals are clipped at
# their 98th percentile so that a few extreme values do not dominate the
# colour scale.

spoints(
  coords[, 1],
  coords[, 2],
  model_data$mu_hat_lin_log,
  main = "",
  xlab = "UTM X coordinate (km)",
  ylab = "UTM Y coordinate (km)",
  asp = 1,
  cex = 0.5
)

q_lim_log <- quantile(abs(model_data$res_ols_log), 0.98, na.rm = TRUE)

res_ols_log_plot <- pmax(
  pmin(model_data$res_ols_log, q_lim_log),
  -q_lim_log
)

spoints(
  coords[, 1],
  coords[, 2],
  res_ols_log_plot,
  main = "",
  xlab = "UTM X coordinate (km)",
  ylab = "UTM Y coordinate (km)",
  asp = 1,
  cex = 0.5
)


# =========================================================
# 4. Empirical variogram of the log-OLS residuals
# =========================================================
# Both the classical estimator and the robust Cressie-Hawkins estimator
# are computed; the robust one is less sensitive to atypical residuals.

datos_vario_log <- model_data %>%
  filter(
    !is.na(res_ols_log),
    !is.na(x),
    !is.na(y)
  )

coordinates(datos_vario_log) <- ~ x + y

bb_log <- bbox(datos_vario_log)

maxlag_log <- sqrt(
  (bb_log["x", "max"] - bb_log["x", "min"])^2 +
    (bb_log["y", "max"] - bb_log["y", "min"])^2
)

vario_res_log <- variogram(
  res_ols_log ~ 1,
  data = datos_vario_log,
  cutoff = maxlag_log / 2,
  width = (maxlag_log / 2) / 20
)

plot(
  vario_res_log,
  main = "",
  xlab = "Distance (km)",
  ylab = "Semivariance",
  pch = 16
)

vario_res_log_rob <- variogram(
  res_ols_log ~ 1,
  data = datos_vario_log,
  cutoff = maxlag_log / 2,
  width = (maxlag_log / 2) / 20,
  cressie = TRUE
)

plot(
  vario_res_log_rob$dist,
  vario_res_log_rob$gamma,
  pch = 1,
  ylim = c(0, 0.30),
  xlab = "Distance (km)",
  ylab = "Semivariance"
)

# =========================================================
# 5. Directional variograms (anisotropy check)
# =========================================================

vario_dir_log <- variogram(
  res_ols_log ~ 1,
  data = datos_vario_log,
  cutoff = maxlag_log / 3,
  width = (maxlag_log / 3) / 15,
  alpha = c(0, 45, 90, 135),
  cressie = TRUE
  
)

plot(
  vario_dir_log,
  main = ""
)

print(vario_dir_log)

# Approximate sill per direction (mean of the last lags), as a quick
# numeric read of the directional plot.
vario_dir_log %>%
  group_by(dir.hor) %>%
  arrange(dist) %>%
  summarise(
    sill_aprox   = mean(tail(gamma, 4)),
    gamma_max    = max(gamma),
    n_pares_tot  = sum(np),
    n_pares_cola = sum(tail(np, 4)),
    .groups = "drop"
  )


# =========================================================
# 6. Spatial independence test
# =========================================================
# Tests the residuals from the linear trend against the independent
# reference model; rejecting independence motivates the dependence-
# corrected methodology.

sm.variogram(
  coords,
  model_data$res_ols_log,
  model = "independent"
)

