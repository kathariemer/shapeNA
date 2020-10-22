#' Computing M-Estimators of Shape for Data Without Missing Values
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
#' @usage powerShape(x, alpha, center = NULL, normalization = c("det", "trace", "one"), maxiter = 1e4, eps = 1e-6)
#' @usage tylerShape(x, center = NULL, normalization = c("det", "trace", "one"), maxiter = 1e4, eps = 1e-6)
#' @usage classicShape(x, center = NULL, normalization = c("det", "trace", "one"), maxiter = 1e4, eps = 1e-6)
#'
#' @param x A numeric data matrix or data.frame without missing data.
#' @param alpha Tail index, a numeric value from the interval [0, 1]. Determines
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
#' @return A `shapeNA` object with a shape estimate `S` and either a center `mu`,
#'   which was either the supplied `center` vector or has been estimated.
#'
#' @return a shapeNA object, which contains a shape and center estimate
#' @export
#'
#' @references Frahm, G., & Jaekel, U. (2010). A generalization of Tyler’s M-estimators to the case of incomplete data. Computational Statistics & Data Analysis, 54(2), 374-393.
#' @references Frahm, G., Nordhausen, K., & Oja, H. (2020). M-estimation with incomplete and dependent multivariate data. Journal of Multivariate Analysis, 176, 104569.
#'
#' @seealso \link{powerShapeNA}
#' @seealso \link{tylerShapeNA}
#' @seealso \link{classicShapeNA}
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
