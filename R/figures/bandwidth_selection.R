# =========================================================
# bandwidth_selection.R
# Figure: local linear estimator for three bandwidths
# Packages: base R only
# =========================================================

if (!dir.exists("figures")) dir.create("figures")

set.seed(6)
n <- 100

# sparse design, ORDERED by x (so dependence is by distance, not by index)
x <- sort(runif(n, -5, 5))

# true function with clear curvature
m <- function(x) 2 * sin(x)

# distance-based dependent errors (1D analogue of the spatial model)
sigma     <- 0.8            # error SD
cor_range <- 0.25           # correlation range
D   <- as.matrix(dist(x))
Sig <- sigma^2 * exp(-D / cor_range)
e   <- as.vector(t(chol(Sig)) %*% rnorm(n))
y   <- m(x) + e

# local linear estimator at a point
ll <- function(x0, h) {
  w  <- dnorm((x - x0) / h)
  X  <- cbind(1, x - x0)
  Xw <- X * w
  beta <- tryCatch(solve(crossprod(Xw, X), crossprod(Xw, y)),
                   error = function(e) c(NA, NA))
  beta[1]
}

# three bandwidths set by hand (undersmooth / middle / oversmooth)
h_small <- 0.1
h_med   <- 0.4
h_large <- 1.5

xx  <- seq(min(x), max(x), length = 400)
fit <- function(h) sapply(xx, ll, h = h)

pdf("figures/bandwidth_selection.pdf", width = 7, height = 5)
plot(x, y, pch = 1, xlab = "", ylab = "")
curve(m, add = TRUE, lwd = 2)               # true function (black)
lines(xx, fit(h_small), col = "green")      # undersmoothing
lines(xx, fit(h_med),   col = "red")        # middle
lines(xx, fit(h_large), col = "blue")       # oversmoothing
dev.off()
  
