#' @export
classicShape <- function(x, center = NULL,  normalization = c("det", "trace", "one"), maxiter = 1e4, eps = 1e-6) {
  if (any(is.na(x))) {
    stop("Missing values found. Use classicShapeNA().")
  }
  try(
    res <- powerShape(x, 0, center, normalization, maxiter, eps)
  )
  res$call <- match.call()
  return(res)
}

