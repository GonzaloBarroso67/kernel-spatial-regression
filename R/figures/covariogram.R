# =========================================================
# covariogram.R
# Figure: exponential variogram and covariogram with nugget
# Packages: base R only
# =========================================================

if (!dir.exists("figures")) dir.create("figures")

# Parameters
phi    <- 0.2            # scale parameter -> practical range 3*phi = 0.6
sigma2 <- 1             # sill = variance = C(0)
c0     <- 0.2            # nugget
ps     <- sigma2 - c0    # partial sill = 0.8
pr     <- 3 * phi        # practical range = 0.6

h     <- seq(0, 0.8, length.out = 500)
h_pos <- h[h > 0]

# Exponential variogram with nugget (h > 0)
gamma_h <- c0 + ps * (1 - exp(-h_pos / phi))

# Covariogram C(h) = sigma2 - gamma(h)
C_h <- sigma2 - gamma_h

pdf("figures/covariogram.pdf", width = 6, height = 4.5)
par(mar = c(4.5, 4.5, 1, 1))

plot(h_pos, gamma_h, type = "l", lwd = 2,
     xlim = c(0, 0.8), ylim = c(0, sigma2 + 0.05),
     xlab = "Distance", ylab = expression(gamma(h)), las = 1)
lines(h_pos, C_h, lty = 2, lwd = 2)

abline(h = sigma2, col = "gray50", lty = 3)
abline(h = 0,      col = "gray50", lty = 3)
abline(v = pr,     col = "gray50", lty = 3)
abline(v = 0,      col = "gray50", lty = 3)

points(0, sigma2, pch = 16, cex = 1.2)
text(0, sigma2, expression(sigma^2), pos = 4, cex = 0.9)
points(0, c0, pch = 1, cex = 1.2, lwd = 2)
text(0, c0, expression(c[0]), pos = 4, cex = 0.9)
text(pr, 0.00, expression(a[p]), pos = 2, cex = 0.9)

dev.off()

