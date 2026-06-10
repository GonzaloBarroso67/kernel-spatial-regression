# =========================================================
# Script 01 - Exploratory analysis
# =========================================================
# AEMET annual precipitation, mainland Spain, 2024.
# Coordinates: UTM 30N (EPSG:25830) in kilometres.
#
# This script reads the frozen dataset produced by script 00. The
# mainland restriction and the collapsing of duplicate locations are
# already done there, so here the data are only described, not
# transformed. Paths are relative to the project root.
# =========================================================


# =========================================================
# 1. Load libraries
# =========================================================

library(dplyr)
library(sf)
library(npsp)
library(e1071)
library(ggplot2)


# =========================================================
# 2. Read frozen dataset and prepare modelling variables
# =========================================================
# The frozen dataset already has one row per unique location, so the
# duplicate check below is only a sanity check and should return FALSE.

aemet_annual_precip_utm <- readRDS(
  "data/aemet_2024_annual_precip_points_sf.rds"
)

model_data <- st_drop_geometry(aemet_annual_precip_utm)

any(duplicated(model_data[, c("x", "y")]))

x <- model_data$x
y <- model_data$y

coords <- as.matrix(model_data[, c("x", "y")])

z <- model_data$annual_precip
log_z <- model_data$log_annual_precip


# =========================================================
# 3. Basic dataset information
# =========================================================

cat("Number of stations:", nrow(model_data), "\n")
cat("Number of provinces:", n_distinct(model_data$provincia_inventory), "\n")

summary(model_data$n_valid_prec)


# =========================================================
# 4. Exploratory analysis of the response: original scale
# =========================================================

summary(z)
sd(z, na.rm = TRUE)
var(z, na.rm = TRUE)

skewness(z, na.rm = TRUE)
kurtosis(z, na.rm = TRUE)


# =========================================================
# 5. Exploratory analysis of the response: log scale
# =========================================================

summary(log_z)
sd(log_z, na.rm = TRUE)
var(log_z, na.rm = TRUE)
skewness(log_z, na.rm = TRUE)
kurtosis(log_z, na.rm = TRUE)

summary(model_data[, c("x", "y")])


# =========================================================
# 6. Spatial plots with spoints
# =========================================================

spoints(
  coords[, 1],
  coords[, 2],
  z,
  main = "Annual precipitation at AEMET stations\nMainland Spain, 2024",
  xlab = "UTM X coordinate (km)",
  ylab = "UTM Y coordinate (km)",
  asp = 1,
  cex = 0.5
)

spoints(
  coords[, 1],
  coords[, 2],
  log_z,
  main = "",
  xlab = "UTM X coordinate (km)",
  ylab = "UTM Y coordinate (km)",
  asp = 1,
  cex = 0.5
)


# =========================================================
# 7. Univariate plots: original scale
# =========================================================

dens <- density(z, na.rm = TRUE)

hist(
  z,
  freq = FALSE,
  main = "",
  xlab = "Annual precipitation (mm)",
  col = "grey80",
  border = "white",
  ylim = c(0, max(dens$y))
)
lines(dens, lwd = 2)

boxplot(
  z,
  main = "",
  ylab = "Annual precipitation (mm)"
)


# =========================================================
# 8. Univariate plots: log scale
# =========================================================

hist(
  log_z,
  freq = FALSE,
  main = "",
  xlab = "log(annual precipitation)",
  col = "grey80",
  border = "white"
)
lines(density(log_z, na.rm = TRUE), lwd = 2)

boxplot(
  log_z,
  main = "",
  ylab = "log(annual precipitation)"
)

# =========================================================
# 9. Exploratory analysis versus projected coordinates
# Log scale
# =========================================================

old.par <- par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))

plot(
  x, log_z,
  pch = 16,
  cex = 0.6,
  col = rgb(0, 0, 0, 0.35),
  xlab = "UTM X coordinate (km)",
  ylab = "log(annual precipitation)",
  main = ""
)

lines(
  lowess(x, log_z),
  lwd = 2,
  col = "blue"
)

plot(
  y, log_z,
  pch = 16,
  cex = 0.6,
  col = rgb(0, 0, 0, 0.35),
  xlab = "UTM Y coordinate (km)",
  ylab = "log(annual precipitation)",
  main = ""
)

lines(
  lowess(y, log_z),
  lwd = 2,
  col = "blue"
)

par(old.par)


# =========================================================
# 10. Descriptive summary by province
# =========================================================

province_summary <- model_data %>%
  group_by(provincia_inventory) %>%
  summarise(
    n_stations = n(),
    mean_precip = mean(annual_precip, na.rm = TRUE),
    median_precip = median(annual_precip, na.rm = TRUE),
    sd_precip = sd(annual_precip, na.rm = TRUE),
    min_precip = min(annual_precip, na.rm = TRUE),
    max_precip = max(annual_precip, na.rm = TRUE),
    mean_log_precip = mean(log_annual_precip, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_precip))

province_summary


# =========================================================
# 11. Wettest and driest stations
# =========================================================

wettest_stations <- model_data %>%
  arrange(desc(annual_precip)) %>%
  select(
    indicativo,
    provincia_inventory,
    annual_precip,
    log_annual_precip,
    n_valid_prec
  ) %>%
  head(10)

wettest_stations


driest_stations <- model_data %>%
  arrange(annual_precip) %>%
  select(
    indicativo,
    provincia_inventory,
    annual_precip,
    log_annual_precip,
    n_valid_prec
  ) %>%
  head(10)

driest_stations


# =========================================================
# 12. Basic distance analysis between stations
# =========================================================
# Since coords contains x and y in kilometres, the distances
# computed here are also in kilometres.

dist_matrix <- as.matrix(dist(coords))

nearest_dist <- apply(
  dist_matrix + diag(Inf, nrow(dist_matrix)),
  1,
  min
)

summary(nearest_dist)

hist(
  nearest_dist,
  main = "Nearest-neighbour distance between stations",
  xlab = "Distance (km)",
  col = "grey80",
  border = "white"
)


# =========================================================
# 13. Missing-data diagnostic
# =========================================================
# Check how missing days are distributed over the year, to assess the
# rescaling-to-366 assumption (script 00). If gaps concentrated in wet
# or dry months, the rescaling could bias the annual estimate. Missing
# days here include the few values coded as "Acum", which are treated as
# missing in script 00. Counts cover the stations that pass the quality
# filter.

missing_by_station_month <- readRDS(
  "data/aemet_2024_missing_by_station_month.rds"
)

# Total number of "Acum" values across all final stations.
sum(missing_by_station_month$n_acum)

# Overall distribution of missing days by month (summed over stations).
missing_by_month <- missing_by_station_month %>%
  group_by(month) %>%
  summarise(n_missing = sum(n_missing), .groups = "drop")

ggplot(missing_by_month, aes(x = month, y = n_missing)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks = 1:12) +
  xlab("Month") +
  ylab("Total missing days") +
  ggtitle("Distribution of missing days by month")

# Station-by-station pattern for the stations with the most gaps
# (between 20 and 66 missing days; 66 is the maximum allowed by the
# quality filter). This shows whether their gaps concentrate in
# particular months or are spread across the year.
stations_high_missing <- missing_by_station_month %>%
  group_by(indicativo) %>%
  summarise(total_missing = sum(n_missing), .groups = "drop") %>%
  filter(total_missing >= 20, total_missing <= 66)

missing_pattern <- missing_by_station_month %>%
  filter(indicativo %in% stations_high_missing$indicativo)

ggplot(missing_pattern, aes(x = month, y = n_missing)) +
  geom_bar(stat = "identity") +
  facet_wrap(~indicativo) +
  scale_x_continuous(breaks = 1:12) +
  xlab("Month") +
  ylab("Missing days") +
  ggtitle("Missing days by month for high-missing stations")
