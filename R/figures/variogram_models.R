# =========================================================
# variogram_models.R
# Figure: parametric variogram models (Gaussian, Exponential, Spherical)
# Packages: gstat
# =========================================================
library(gstat)

if (!dir.exists("figures")) dir.create("figures")

# Parameters
nugget <- 0
psill  <- 1
sill   <- nugget + psill
practical_range <- 1000
maxdist <- 1500

# One entry per model, in plotting order: gstat code, range, colour
models <- list(
  Gaussian    = list(code = "Gau", range = practical_range / sqrt(3), col = "darkgreen"),
  Exponential = list(code = "Exp", range = practical_range / 3,       col = "blue"),
  Spherical   = list(code = "Sph", range = practical_range,           col = "red")
)

# Theoretical semivariogram curves
curves <- lapply(models, function(p)
  variogramLine(vgm(psill, p$code, p$range, nugget), maxdist = maxdist, n = 1000))

pdf("figures/variogram_models.pdf", width = 8.5, height = 4.5)
par(mar = c(5, 5, 2, 10), xpd = FALSE)

plot(NA, xlim = c(0, maxdist), ylim = c(0, 1.15), las = 1, xaxs = "i", yaxs = "i",
     xlab = "Separation distance", ylab = "Semivariance")
for (i in seq_along(models))
  lines(curves[[i]]$dist, curves[[i]]$gamma, lwd = 2, col = models[[i]]$col)

segments(0, sill, maxdist, sill, lty = 2, col = "red")
segments(practical_range, 0, practical_range, 1.15, lty = 2, col = "red")

text(90, sill + 0.075, "Sill", cex = 1)
text(practical_range + 35, 0.08, "Practical range", pos = 4, cex = 1)

par(xpd = NA)
legend(maxdist + 60, 0.70, bty = "n", lwd = 2,
       legend = names(models), col = sapply(models, `[[`, "col"))

dev.off()

