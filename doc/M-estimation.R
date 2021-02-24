## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(shapeNA)
library(mvtnorm)

## -----------------------------------------------------------------------------
set.seed(1)
mu <- c(0, 0)
S <- matrix(c(1, 0.5, 0.5, 2), ncol = 2)
X <- rmvt(500, sigma = S, df = 8)

## -----------------------------------------------------------------------------
muHat <- colMeans(X)
SigmaHat <- cov(X)

## -----------------------------------------------------------------------------
plot(X)
lines(shapeNA:::ellipse(mu, S), col = 2, lwd = 2)
points(t(mu), col = 2, pch = 3, lwd = 2)
lines(shapeNA:::ellipse(muHat, SigmaHat), col = 3, lty = 3, lwd = 2)
points(t(muHat), col = 3, pch = 4, lwd = 2)

## -----------------------------------------------------------------------------
plot(X)
for (i in 4:1) {
  res <- powerShape(X, alpha = (i-1)/4)
  lines(shapeNA:::ellipseShape(res), col = i+1, lwd = 2)
  points(t(res$mu), col = i+1, pch = i+1)
}
legend(2.5, -2, paste('alpha', (3:0/4)), col = 2:5, lty = 1, lwd = 2)

## -----------------------------------------------------------------------------
summary(res)

## -----------------------------------------------------------------------------
all.equal(res$mu, muHat)
all.equal(shapeNA:::toCov(res), SigmaHat)

