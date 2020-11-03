#' Compute M-Estimators of Shape for Data Without Missing Values
#'
#' `powerShape`, `tylerShape` and `classicShape` compute
#' power M-estimators of shape. Using `powerShape` and `classicShape`
#' it is even possible to compute M-estimators of covariance matrices.
#' These functions also compute estimates of location if no `center` is supplied.
#'
#' For multivariate normally distributed data, `classicShape` is an
#' ML-estimator. This is a special case of the power M-estimator with tail index
#' `alpha` = 0 and returns the empirical covariance matrix and the empirical mean
#' vector.
#'
#' `tylerShape` maximizes the likelihood function after projecting the observed
#' data of each individual onto the unit hypersphere, in which case we obtain
#' an angular central Gaussian distribution. This is a special case of the power
#' M-estimator with tail index `alpha` = 1 and returns Tyler's M-estimator of
#' scatter and an affine equivariant multivariate median.
# TODO: Hettmansperger and randles
#'
#' `powerShape` requires an additional parameter, the so-called tail index `alpha`.
#' For asymptotic normality, this index should be chosen taking into consideration
#' the data. For heavy tailed data, the index should be closer to 1, for light
#' tailed data the index should be chosen closer to 0.
#'
#' @aliases powerShape
#' @aliases tylerShape
#' @aliases classicShape
#'
#' @usage powerShape(x, alpha, center = NULL,
#'     normalization = c("det", "trace", "one"), maxiter = 1e4, eps = 1e-6)
#' @usage tylerShape(x, center = NULL,
#'     normalization = c("det", "trace", "one"), maxiter = 1e4, eps = 1e-6)
#' @usage classicShape(x, center = NULL,
#'     normalization = c("det", "trace", "one"), maxiter = 1e4, eps = 1e-6)
#'
#' @param x A numeric data matrix or data.frame without missing data.
#' @param alpha Tail index, a numeric value from the interval `[0, 1]`. Determines
#'   the power function. For more information see 'Details'.
#' @param center An optional vector of the data's center, if NULL the center
#'   will be estimated simultaneously to the shape estimate.
#' @param normalization A string, determines scale of returned shape estimate.
#'   The possible values are \itemize{
#'     \item{'det'} s.t. the returned shape estimate has determinant 1.
#'     \item{'trace'} s.t. the returned shape estimate has trace `p`.
#'     \item{'one'} s.t. the returned shape estimate's first entry is 1.
#'   }
#' @param maxiter A positive integer, restricting the maximum number of iterations.
#' @param eps A numeric, specifying tolerance level of when the iteration stops.
#'
#' @return A list with class 'shapeNA' containing the following components:
#' \describe{
#'   \item{S}{the estimated shape matrix.}
#'   \item{scale}{the scale with which the shape matrix may be scaled to obtain a scatter estimate. If `alpha` == 1, then this value is meaningless.}
#'   \item{mu}{the location parameter, either provided by the user or estimated.}
#'   \item{alpha}{the tail index with which the Power M-estimator has been called.}
#'   \item{naBlocks}{`NULL`, since `powerShape` operates only on complete data.}
#'   \item{iterations}{number of computed iterations before convergence.}
#'   \item{call}{the matched call.}
#' }
#'
#' @export
#'
#' @references Tyler, D.E. (1987). A Distribution-Free M-Estimator of Multivariate Scatter. The Annals of Statistics, 15, 234.251. <doi:10.1214/aos/1176350263>.
#' @references Frahm, G., & Jaekel, U. (2010). A generalization of Tyler's M-estimators to the case of incomplete data. Computational Statistics & Data Analysis, 54, 374-393. <doi:10.1016/j.csda.2009.08.019>.
#' @references Frahm, G., Nordhausen, K., & Oja, H. (2020). M-estimation with incomplete and dependent multivariate data. Journal of Multivariate Analysis, 176, 104569. <doi:10.1016/j.jmva.2019.104569>.
#'
#' @seealso [powerShapeNA], [tylerShapeNA] and [classicShapeNA] for the
#'     corresponding functions for data with missing values.
#'
#' @examples
#'     ## Generate example data
#'     S <- toeplitz(c(1, 0.1))
#'     x <- mvtnorm::rmvt(100, S)
#'     ## Compute some M-estimators
#'     res0 <- classicShape(x, center = c(0, 0))
#'     res1 <- powerShape(x, alpha = 0.67, normalization = 'one')
#'     res2 <- tylerShape(x, normalization = 'trace')
#'     ## Get location estimates
#'     res1$mu
#'     res2$mu
#'     ## Get shape estimates
#'     res0$S
#'     res1$S
#'     res2$S
#'     ## Print summary
#'     summary(res0)
powerShape <- function(x, alpha, center = NULL, normalization = c("det", "trace", "one"), maxiter = 1e4, eps = 1e-6) {
  if (any(is.na(x))) {
    stop("Missing values found. Use powerShapeNA()")
  }
  fctCall <- match.call()
  powerfct <- powerFunction(alpha)
  scatterNormFct <- normalizationFunction(normalization)

  try(
    {
    if (is.null(center)) {
      res <- mestimator_mean_cov(x, powerfct, scatterNormFct, maxiter, eps)
    } else {
      xCentered <- sweep(x, 2, center)
      zeroEntry <- rowSums(xCentered) == 0
      if (any(zeroEntry)) {
        xCentered <- xCentered[!zeroEntry, ]
        z <- sum(zeroEntry)
        if (z==1) {
          message("Found ", z, " observation coinciding with given center.")
        } else {
          message("Found ", z, " observations coinciding with given center.")
        }
      }
      res <- mestimator_cov(xCentered, powerfct, scatterNormFct, maxiter, eps)
      res$mu <- center
    }
    res$alpha <- alpha
    res$call <- fctCall
    return(res)
    }
  )
}
