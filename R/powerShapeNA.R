#' Computing M-Estimators of Shape for Data With Missing Values
#'
#' `powerShapeNA`, `tylerShapeNA` and `classicShapeNA` compute
#' power M-estimators of shape for data with missing values. The underlying
#' missingness mechanism should be MCAR. These functions also compute estimates
#' of location if no `center` is supplied.
#'
#' For multivariate normally distributed data, `classicShapeNA` is an
#' ML-estimator. This is a special case of the power M-estimator with tail index
#' `alpha` = 0 and returns the empirical covariance matrix and the empirical mean
#' vector.
#'
#' `tylerShapeNA` maximizes the likelihood function after projecting the observed
#' data of each individual onto the unit hypersphere, in which case we obtain
#' an angular central Gaussian distribution. This is a special case of the power
#' M-estimator with tail index `alpha` = 1 and returns Tyler's M-estimator of
#' scatter and an affine equivariant multivariate median.
# TODO: Hettmansperger and randles
#'
#' `powerShapeNA` requires an additional parameter, the so-called tail index `alpha`.
#' For asymptotic normality, this index should be chosen taking into consideration
#' the data. For heavy tailed data, the index should be closer to 1, for light
#' tailed data the index should be chosen closer to 0.
#'
#' @aliases powerShapeNA
#' @aliases classicShapeNA
#' @aliases tylerShapeNA
#'
#' @usage powerShapeNA(x, alpha, center = NULL, normalization = c("det", "trace", "one"),
#'          maxiter = 1e4, eps = 1e-6)
#' @usage tylerShapeNA(x, center = NULL, normalization = c("det", "trace", "one"),
#'           maxiter = 1e4, eps = 1e-6)
#' @usage classicShapeNA(x, center = NULL, normalization = c("det", "trace", "one"),
#'          maxiter = 1e4, eps = 1e-6)
#'
#' @param x A data matrix or data.frame with missing data and `p > 2` columns.
#'   Representing sample from generalized elliptical distribution and MCAR missingness
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
#' @export
#'
#' @references Frahm, G., & Jaekel, U. (2010). A generalization of Tyler’s M-estimators to the case of incomplete data. Computational Statistics & Data Analysis, 54, 374-393. <doi:10.1016/j.csda.2009.08.019>.
#' @references Frahm, G., Nordhausen, K., & Oja, H. (2020). M-estimation with incomplete and dependent multivariate data. Journal of Multivariate Analysis, 176, 104569. <doi:10.1016/j.jmva.2019.104569>.
#'
#' @examples
#'     ## Generate a data set with missing values
#'     x <- mvtnorm::rmvt(100, toeplitz(seq(1, 0.1, length.out = 3)), df = 5)
#'     y <- mice::ampute(x, mech='MCAR')$amp
#'
#'     ## Compute some M-estimators
#'     res0 <- classicShapeNA(y, center = c(0, 0, 0))
#'     res1 <- powerShapeNA(y, alpha = 0.67, normalization = 'one')
#'     res2 <- tylerShapeNA(y, normalization = 'trace')
#'
#'     ## Get location estimates
#'     res1$mu
#'     res2$mu
#'     ## Get shape estimates
#'     res0$S
#'     res1$S
#'     res2$S
#'
#'     ## Print summary
#'     summary(res0)
#'     ## Inspect missingness pattern
#'     plot(res0$naBlocks)
#'     barplotMissProp(res0$naBlocks)
powerShapeNA <- function(x, alpha, center = NULL, normalization = c("det", "trace", "one"), maxiter = 1e4, eps = 1e-6) {
  if (!any(is.na(x))) {
    stop("No missing values found. Use powerShape().")
  }
  fctCall <- match.call()
  powerfct <- powerFunction(alpha)
  scatterNormFct <- normalizationFunction(normalization)

  try(
    {
    if (is.null(center)) {
      res <- mestimator_mean_cov(x, powerfct, scatterNormFct, maxiter,eps)
    } else {
      xCentered <- sweep(x, 2, center)
      zeroEntry <- rowSums(xCentered, na.rm = TRUE) == 0
      if (any(zeroEntry)) {
        xCentered <- xCentered[!zeroEntry, ]
        z <- sum(zeroEntry)
        if (z==1) {
          message("Found ", z, " observation coinciding with given center.")
        } else {
          message("Found ", z, " observations coinciding with given center.")
        }
      }
      res <- mestimator_cov(xCentered, powerfct, scatterNormFct, maxiter,eps)
      res$mu <- center
    }
    res$alpha <- alpha
    res$call <- fctCall
    return(res)
    }
  )
}
