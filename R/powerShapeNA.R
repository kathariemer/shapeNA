#' Computing M-Estimators of Shape for Data With Missing Values
#'
#' todo: add description here
#'
#' @aliases powerShapeNA
#' @aliases classicShapeNA
#' @aliases tylerShapeNA
#'
#' @usage powerShapeNA(x, alpha, center = NULL, normalization = c("det", "trace", "one"), 
#'          maxiter = 1e4, tol = 1e-6)
#' @usage tylerShapeNA(x, center = NULL, normalization = c("det", "trace", "one"), 
#'           maxiter = 1e4, tol = 1e-6)
#' @usage classicShapeNA(x, center = NULL, normalization = c("det", "trace", "one"),
#'          maxiter = 1e4, tol = 1e-6)
#'
#' @param x data matrix or data.frame with missing data and more than 2 columns. Representing sample from continuous distribution and MCAR missingness
#' @param alpha numeric, determines power function
#' @param center optional vector of center, if NULL the center will be estimated simultaneously to the shape estimate
#' @param normalization string, determines scale of returned shape estimate
#' @param maxiter integer, maximum number of iterations
#' @param tol numeric, tolerance level
#'
#' @return shape and center estimate
#' @export
#'
#' @references Frahm, G., & Jaekel, U. (2010). A generalization of Tyler’s M-estimators to the case of incomplete data. Computational Statistics & Data Analysis, 54, 374-393. <doi:10.1016/j.csda.2009.08.019>.
#' @references Frahm, G., Nordhausen, K., & Oja, H. (2020). M-estimation with incomplete and dependent multivariate data. Journal of Multivariate Analysis, 176, 104569. <10.1016/j.jmva.2019.104569>.
#'
#' @examples
#'     ## generate data set with missing values
#'     x <- mvtnorm::rmvt(100, toeplitz(seq(1, 0.1, length.out = 3)), df = 5)
#'     y <- mice::ampute(x, mech='MCAR')$amp
#'     ## compute M-estimate
#'     res <- powerShapeNA(y, alpha = 0.5)
#'     summary(res)
powerShapeNA <- function(x, alpha, center = NULL, normalization = c("det", "trace", "one"), maxiter = 1e4, tol = 1e-6) {
  if (!any(is.na(x))) {
    stop("No missing values found. Use powerShape().")
  }
  fctCall <- match.call()
  powerfct <- powerFunction(alpha)
  scatterNormFct <- normalizationFunction(normalization)

  try(
    {
    if (is.null(center)) {
      res <- mestimator_mean_cov(x, powerfct, scatterNormFct, maxiter, tol)
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
      res <- mestimator_cov(xCentered, powerfct, scatterNormFct, maxiter, tol)
      res$mu <- center
    }
    res$alpha <- alpha
    res$call <- fctCall
    return(res)
    }
  )
}
