mestimator_cov <- function(x, ...) {
  UseMethod("mestimator_cov")
}

mestimator_cov.default <- function(x, ...) {
  stop("No implementation for object of provided class. Please supply a data matrix of observations.")
}

mestimator_cov.data.frame <- function(x, ...) {
  mestimator_cov.matrix(as.matrix(x), ...)
}

mestimator_cov.matrix <- function(x, powerfct, normalization, maxiter=1e4, tol=1e-6, ...) {
  if (tol <= 0 || maxiter <= 0) {
    stop("Nonpositive arguments maxiter or tol.")
  }
  if (any(is.na(x))) {
    return(mestimator_cov.naBlocks(naBlocks(x), powerfct, normalization, maxiter, tol, ...))
  }

  i <- 0
  dist <- 2*tol
  n <- nrow(x)
  p <- ncol(x)
  S0 <- diag(p)
  while (i < maxiter && dist > tol) {
    # squared mahalanobis distances xi
    xi <- rowSums(x * t(solve(S0, t(x)))) # define weights
    try(w <- powerfct(xi, p = p)$w)
    S <- (crossprod(x, w*x))/n

    dist <- norm(S-S0, "M")
    S0 <- S
    i <- i+1
  }

  if (i >= maxiter) {
    warning(paste("No convergence in", i, "steps."))
  }

  shape <- normalization(S)
  res <- list(S=shape$S, scale=shape$scale, mu=NULL, alpha=NULL, iterations=i, naBlocks=NULL)
  class(res) <- "shapeNA"
  return(res)
}

mestimator_cov.naBlocks <- function(x, powerfct, normalization, maxiter, tol, ...) {
  if (tol <= 0 || maxiter <= 0) {
    stop("Nonpositive arguments maxiter or tol.")
  }
  i <- 0
  dist <- 2*tol
  y <- x$data
  n <- nrow(y)
  p <- ncol(y)

  # use closure
  covOfSubset <- function(x, S, varCount, n) {
      Sinv <- solve(S)
      if (!isSymmetric(Sinv)) {
        Sinv <- (Sinv + t(Sinv))/2
      }
      xi <- rowSums((x%*%Sinv)* x)
      try(
        w <- powerfct(xi, p, varCount)$w
      )
      a <- (crossprod(x, w*x))
      return(Sinv %*% a %*% Sinv - n * Sinv)
  }

  S0 <- diag(p)
  blockIdx <- x$N
  blockPattern <- x$P
  while (i < maxiter && dist > tol) {
    # split data into blocks of same missingness
    rows <- 1:blockIdx[1]
    tryCatch(
    a_Sigma <- covOfSubset(
      y[rows, , drop = FALSE], S0, n=length(rows), varCount = p)
      , error = function(e) {
        message(paste("Iteration", i))
        message(e)
        stop("Matrix not positive definite")
      }
    )

    for (j in 2:length(blockIdx)) {
      rows <- (blockIdx[j-1]+1):blockIdx[j]
      cols <- asBinaryVector(blockPattern[j], p)
      tryCatch(
      aux <- covOfSubset(
        y[rows, cols, drop = FALSE], S0[cols, cols],
        n = length(rows), varCount = sum(cols))
      , error = function(e) {
        message(paste("Iteration", i, "block", j))
        message(e)
        stop("Matrix not positive definite")
      }
      )
      a_Sigma[cols, cols] <- a_Sigma[cols, cols] + aux
    }
    # adjust scatter estimate
    a_Sigma <- S0 %*% (a_Sigma/n) %*% S0
    dist <- norm(a_Sigma, "M")
    S0 <- S0 + a_Sigma
    i <- i+1
  }

  if (i >= maxiter) {
    warning(paste("No convergence in", i, "steps."))
  }

  ogOrder <- order(x$permutation)
  S0 <- S0[ogOrder, ogOrder]

  shape <- normalization(S0)
  res <- list(S=shape$S, scale=shape$scale, mu=NULL, alpha=NULL, iterations=i, naBlocks=x)
  class(res) <- "shapeNA"
  return(res)
}

