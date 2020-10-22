#' Reorder Data With Missing Values
#'
#' Reorder a data set with `NA` entries to form blocks of missing values. The resulting
#' data will have increasing missingness along the rows and along the columns. The rows
#' are ordered s.t. the first block consists of complete observations, and
#' the following blocks are ordered from most frequent missingness pattern to least
#' frequent missingness pattern.
#' This method may fail, as it has been designed as a preprocessing step for
#' shape estimations.
#'
#' @param data A matrix with `NA` values
#' @param cleanup A logical flag. If `TRUE`, observations with less than 2 responses
#'    are discarded
#' @param plot A logical flag. If true, a plot of the missingness pattern is produced.
#'
#' @return a `naBlocks` object, which is a list with \itemize{
#'     \item `data` the reordered data matrix
#'     \item `permutation` the permutation of the columns, which was applied to reorder
#'         the columns according to the number of `NA`s
#'     \item `rowPermutation` the permutation of the rows, which generates the blocks
#'     \item `N` a vector of all row indices. Each row number points to the beginning
#'         of a new missingness pattern
#'     \item `D` a vector specifying the missingness pattern for each block.
#'     \item `P` a vector specifying the number of observed variables per block.
#'     \item `kn` a vector specifying the percentage of observed responses per variable.
#' }
naBlocks <- function(data, cleanup=TRUE, plot=FALSE) {
  ## remove observations without responses or only 1 response ##
  x <- data
  if (cleanup) {
    hasValue <- !is.na(x) # R matrix
    #rmIdx <- apply(hasValue, 1, function(x) sum(x) <=1)
    rmIdx <- rowSums(hasValue)<=1
    if (any(rmIdx)) {
      x <- data[!rmIdx, ]
      message(paste("Removed", sum(rmIdx), "observations with less than 2 responses each."))
    }
  }
  ## compute dimensions of cleaned up data ##
  n <- nrow(x)
  p <- ncol(x)
  if (n < 2) {
    stop(paste("Only", n, "observations."))
  } else if (n <= p) {
    stop('Dimension exceeds sample size.')
  }
  colNames <- colnames(x)
  if (is.null(colNames)) {
    colNames <- paste("V", 1:p, sep = '')
  }
  ## reorder data to form contiguous blocks with the same pattern ##
  hasValue <- !is.na(x) # R matrix
  respPerVar <- colSums(hasValue)
  xIdx <- order(respPerVar, decreasing = TRUE)
  x <- x[, xIdx]
  hasValue <- hasValue[, xIdx]

  ## get indices where patterns change and the number of observations/the corresponding pattern ##
  blockInfo <- getMissingnessBlocks(hasValue)
  rle <- blockInfo$rle
  N <- cumsum(rle$lengths)
  if (length(N) == 1) {
    warning("No missing data.\n")
    return(x)
  }
  D <- rowSums(hasValue[blockInfo$order[N],])
  P <- rle$values  # save pattern info

  x <- x[blockInfo$order, ]

  if (N[1] < D[1]) {
    # first block has fewer responses than rows
    cat("N: ", N, "\nD: ", D, "\n")
     stop('Missingness is too strong! Estimator cannot converge.')
  }
  res <- list(data=as.matrix(x), permutation=xIdx, rowPermutation=blockInfo$order,
              N=N, D=D, P=P, kn=respPerVar/n, colNames = colNames)
  class(res) <- "naBlocks"
  if (plot) {
    plot(res)
  }
  return(res)
}


#' plot missingness pattern of data
#'
#' @param x A `naBlocks` object
#' @param ... additional parameters
#'
#' @export
#' @examples
#'     x <- mvtnorm::rmvt(100, toeplitz(seq(1, 0.1, length.out = 3)), df = 5)
#'     y <- mice::ampute(x, mech='MCAR')$amp
#'     res <- classicShapeNA(y)
#'     plot(res$naBlocks)
plot.naBlocks <- function(x, orderProp=TRUE, ...) {
  idx <- x$N
  bprop <- (idx - c(0, idx[-length(idx)]))/nrow(x$data)
  if (orderProp){
    bIdx <- order(bprop, decreasing = TRUE)
    bprop <- bprop[bIdx]
    plotAsBoolMat(x$data[idx[bIdx],], round(bprop,3))
  } else {
    bIdx <- 1:length(idx)
    plotAsBoolMat(x$data[idx,], round(bprop,3))
  }
  invisible(list(blockProps=bprop, blockIdx=bIdx))
}

# given a number n, return its binary representation as a vector
#
# @param n decimal number
# @param d length of vector
# @param logical whether to return a numeric or logical vector
asBinaryVector <- function(n, d, logical = TRUE) {
  b <- numeric(d)
  for (i in d:1) {
    b[i] <- n %% 2
    n <- n %/% 2
  }
  if (n > 0) {
    warning(paste("Number was larger than 2^",d,"-1. Allow for ~",floor(log2(n))+1," more digits.\n", sep=''))
  }
  if (logical) {
    b <- as.logical(b)
  }
  return(b)
}

# Given a integer, which represents a binary pattern, return a vector, which defines the corresponding indices
#
# @param x integer representing a binary pattern
# @param d maximal exponent, i.e. 2^d > x
#
# @return subset of 1:d, which selects the entries corresponding to coefficient 1 of binary representation of x
binaryToIndex <- function(x, d) {
  b <- numeric(d)
  for (i in d:1) {
    b[i] <- x %% 2
    x <- x %/% 2
  }
  if (x > 0) {
    warning(paste("Number was larger than 2^",d,"-1. Allow for ~",floor(log2(x))+1," more digits.\n", sep=''))
  }
  b <- (1:d)[b==1]
  return(b)
}

# plot boolean matrix as a grid of colored rectangles
# blue fields are TRUE, white fields are FALSE, yellow fields are NA
# @param R boolean matrix
plotAsBoolMat <- function(R, blockNames, ...) {
  # R <- obj$pattern
  color2 <- colors["white"]
  if (any(is.na(R))) {
    color2 <- colors["yellow"]
    B <- is.na(R)
  } else {
    B <- R == 0
  }
  r <- nrow(R)
  c <- ncol(R)
  plot(1, type='n', xlim=c(0,c+1), ylim = c(r+1,0), axes=FALSE, xlab=NA, ylab=NA)
  xleft <- rep(1, r)
  ybottom <- 0:(r-1)
  for (i in 1:c) {
    rowColor <- ifelse(B[,i], color2, colors["blue"])
    rect(xleft, ybottom, xleft+1, ybottom+1, col = rowColor, ...)
    xleft <- xleft+1
  }
  text(0.5+(1:c),r+1, colnames(R))
  if (!missing(blockNames)){
    text(0, (1:r)-0.5, blockNames)
  }
}

# transform numeric input to have values between 0 and 1
# the largest value gets mapped to 0 (black), the smallest to 1 (white)
grayscale <- function(A) {
  m <- min(A)
  M <- max(A)
  return((A-m)/(m-M)+1)
}

#' print method for elements of class shapeNA
#'
#' Only print M-estimates and alpha level
#'
#' @param x A `shapeNA` object
#' @param ... Additional parameters.
#'
#' @export
#' @examples
#'     x <- mvtnorm::rmvt(100, toeplitz(seq(1, 0.1, length.out = 3)), df = 5)
#'     res <- tylerShape(x)
#'     res ## equivalent to call print(res)
print.shapeNA <- function(x, ...) {
  if (is.null(x$mu)) {
    print(list(alpha=x$alpha, S=x$S), ...)
  } else {
    print(list(alpha=x$alpha, S=x$S, mu=x$mu), ...)
  }
}

#' Crude visualization of shape estimate
#'
#' Plot each matrix entry as a cell, with dark cells indicating high values and
#' light values indicate small values.
#'
#' @param x A `shapeNA` oopbject
#' @param message A logical, If `TRUE`, a similar summary is printed in the console.
#' @param ... Additional parameters.
#'
#' @export
#' @examples
#'     x <- mvtnorm::rmvt(100, toeplitz(seq(1, 0.1, length.out = 3)), df = 5)
#'     y <- mice::ampute(x, mech='MCAR')$amp
#'     res <- tylerShapeNA(y)
#'     plot(res)
plot.shapeNA <- function(x, message=TRUE, ...) {
  imageS <- x$S
  p <- ncol(imageS)
  vars <- x$naBlocks$colNames
  ## "flip" image to have 1st entry in lower left corner
  imageS <- t(imageS[p:1,])
  graphics::image(imageS, xaxt = 'n', yaxt = 'n', ...)
  axis(side = 1, at = seq(0, 1, length.out = p),
       labels = vars, lwd = 0)
  if (!is.null(x$naBlocks) && message) {
    df <- data.frame(
      Variables = vars,
      Observed = round(x$naBlocks$kn, 2),
      row.names = NULL)
    format(df)
  }
}

# trace of matrix, i.e. sum of diagonal entries
tr <- function(X) {
  return(sum(diag(X)))
}

# return a scaled scale estimator S, which satisfies trace(S) = p
scaleTr <- function(Sigma) {
  #s <- nrow(Sigma)/sum(diag(Sigma))
  s <- tr(Sigma)/nrow(Sigma)
  if (s == 0) {
    stop("Cannot scale matrix, trace equals 0")
  }
  return(list(S = Sigma/s, scale=s))
}

scaleDet <- function(Sigma) {
  d <- det(Sigma)^(1/ncol(Sigma))
  if (d == 0) {
    stop("Cannot scale matrix, determinant equals 0")
  }
  return(list(S = Sigma/d, scale=d))
}

scaleOne <- function(Sigma) {
  s <- Sigma[1,1]
  if (s == 0) {
    stop("Cannot scale matrix, first entry equals 0")
  }
  return(list(S = Sigma/Sigma[1,1], scale = s))
}

colors <- c(red = "#d55e00",
            purple = "#cc79a7",
            blue = "#0072b2",
            yellow = "#FFCD70",
            green = "#009e73",
            black = "#000000",
            white = "#ffffff")

# use eigen decomposition to compute square root of matrix
# @param Sigma symmetric, pos-definite matrix
matroot <- function(Sigma) {
  A <- eigen(Sigma)
  Gamma <- A$vectors
  D <- diag(sqrt(A$values))
  return(tcrossprod(Gamma%*%D, Gamma))
}

# runlength encoding for missingness matrix, where
# the rows get transformed to an integer corresponding
# to it's missingness pattern (big endian)
# Example: c(1,1,0) <-> 4\*1 + 2\*1 + 1\*0
# @param R reordered missingness matrix
# @return permutation for reordering of data and runlength encoded integer pattern,
#   i.e a missingness of (0,1,0,1) translates to 8*0 + 4*1 + 2*0 * 1*1 = 5
getMissingnessBlocks <- function(R) {
  p2 <- 2^((ncol(R)-1):0)
  pattern <- rowSums(sweep(R, 2, p2, "*"))
  lev <- summary(as.factor(pattern))
  patternF <- factor(pattern, levels = names(sort(lev, decreasing = TRUE)))
  #rowOrder <- order(rowSums(R), pattern, decreasing = c(TRUE, FALSE))
  rowOrder <- order(patternF)
  return(list(order=rowOrder, rle = rle(pattern[rowOrder])))
}

#' Print Missingness Pattern
#'
#' Print the pattern of missingness in the supplied data, as a matrix for 1s,
#' representing a column vector of responses and 0s, representing a column vector
#' of `NA`s.
#'
#' The first row shows the column names. The leftmost column, without column names,
#'  shows the number of rows per block and the rightmost column, titled with `#`
#'  shows the number of observed variables in the block.
#'
#' @param x A `naBlocks` object
#' @param ... Additional parameters.
#'
#' @export
#' @examples
#'     x <- mvtnorm::rmvnorm(200, mean = c(0, 0))
#'     classicShape(x)
print.naBlocks <- function(x, ...) {
  p <- length(x$permutation)
  pattern <- sapply(x$P, asBinaryVector, d=p, logical=FALSE)
  res <- matrix(c(t(pattern), D=p-x$D), nrow=length(x$N))
  rCount <- c(x$N, 0) - c(0, x$N)
  row.names(res) <- rCount[1:length(x$N)]
  colNames <- c(paste("V", x$permutation, sep=""), " #")
  colnames(res) <- colNames
  print(res, ...)
}

#' Barplot Showcasing Missingness Proportion of the Original Data
#'
#' Visualize the proportion of missingness per variable in a barplot.
#'
#' @param obj A `naBlocks` object
#' @param sortNA A logical. If `FALSE` the original variable order is kept,
#'     else the variables are ordered from least to most missingness
#' @param ... Additional graphical arguments for \code{\link[graphics]{barplot}}
#'
#' @export
#' @examples
#'     S <- toeplitz(seq(1, 0.1, length.out = 3))
#'     x <- mvtnorm::rmvt(100, S, df = 5)
#'     y <- mice::ampute(x, mech='MCAR')$amp
#'     res <- classicShapeNA(y)
#'     barplotMissProp(res$naBlocks)
barplotMissProp <- function(obj, sortNA = FALSE, ...) {
  if (is.null(obj$data)) {
    stop('No missing values in the data.')
  }
  mprop <- colSums(is.na(obj$data))/nrow(obj$data)
  if (!sortNA) {
    mprop <- mprop[order(obj$permutation)]
  }
  graphics::barplot(mprop, ...)
  invisible(mprop)
}

#' summary method for class shapeNA
#'
#' @param object an object of class shapeNA, usually from a call to powerShape or similar functions
#' @param ... further arguments
#'
#' @return object of class shapeNA
#' @export
#'
#' @examples
#'     obj <- tylerShape(mvtnorm::rmvt(100, diag(3)))
#'     summary(obj)
summary.shapeNA <- function(object, ...) {
  class(object) <- 'summary.shapeNA'
  return(object)
}

#' print method for class summary.shapeNA
#'
#' @param x object returned from summary.shapeNA
#' @param ... further arguments
#'
#' @return invisibly return NULL
#' @export
#'
#' @examples
#'     obj <- tylerShape(mvtnorm::rmvt(100, diag(3)))
#'     print(summary(obj))
print.summary.shapeNA <- function(x, ...) {
  cat('Call: ')
  print(x$call)
  cat('---\nalpha:', x$alpha)
  cat('\nCenter:\n')
  print(x$mu, digits=3)
  cat('\nShape:\n')
  print(x$S, digits=3)
  if (!is.null(x$naBlocks)) {
    cat('\nMissingness pattern:\n')
    print(x$naBlocks)
  }
  cat('---\nComputed', x$iterations, 'iterations.\n')
  invisible(NULL)
}

# get power function to pass to powerShape/powerShapeNA
powerFunction <- function(alpha) {
  if (alpha > 1 || alpha < 0) {
    stop("Ensure 0 <= alpha <= 1.")
  }
  else if (alpha == 1) {
    powerfct <- power_weight_tyler
  } else if (alpha==0) {
    powerfct <- power_weight_gauss
  } else {
    powerfct <-
      function(xi, p, m)
        power_weight(xi, p, m, alpha = alpha)
  }
  return(powerfct)
}

# get normalization function to pass to powerShape/powerShapeNA
normalizationFunction <- function(normalization) {
  ## determine normalization of scatter estimate ##
  if (length(normalization) >= 1) {
    normalization <- normalization[1]
  }
  normalization <- pmatch(tolower(normalization), c("det", "trace", "one"))
  if (normalization == 1) {
    scatterNormFct <- scaleDet
  } else if (normalization == 2) {
    scatterNormFct <- scaleTr
  } else if (normalization == 3) {
    scatterNormFct <- scaleOne
  } else {
    stop("Unknown method of normalization")
  }
  return(scatterNormFct)
}

# for power M-estimates with alpha < 1: restore covariance estimate
toCov <- function(obj) {
  if (obj$alpha == 1) {
    stop("Covariance matrix for Tyler's M-estimate undefined")
  }
  return(obj$scale * obj$S)
}

#' Ellipse For Covariance Matrix From M-Estimator
#'
#' Given a `shapeNA` object, get (x,y)-coordinates of points, which for the
#' estimated scatter (shape, if `alpha` == 1) and location have mahalanobis
#' distance 1.
#'
#' @param obj A `shapeNA` object
#' @param idx A vector of two distinct indices which specify the submatrix
#' @param n number of points on the ellipse
#'
#' @examples
#'     S <- toeplitz(c(1, 0.3, 0.7))
#'     set.seed(123)
#'     x <- mvtnorm::rmvt(100, S, df = 3)
#'     obj <- powerShape(x, alpha = 0.85)
#'     ## Plot variables 1 and 3
#'     idx <- c(1,3)
#'     plot(x[, idx])
#'
#'     ## Plot projection of true covariance matrix
#'     lines(shapeNA:::ellipse(S[idx, idx], c(0, 0)), col = 2)
#'     ## Plot base-R estimate
#'     lines(shapeNA:::ellipse(cov(x[, idx]), c(0, 0)), col = 3, lty = 2)
#'     ## Plot M-estimate
#'     lines(shapeNA:::ellipseShape(obj, idx=idx), col = 4, lty = 2)
ellipseShape <- function(obj, idx=1:2, n = 250) {
  S <- if (obj$alpha == 1) {
    obj$S
  } else {
    toCov(obj)
  }
  S <- S[idx, idx]
  mu <- (obj$mu)[idx]
  return(ellipse(S, mu, n))
}

#' Ellipse For Covariance Matrix
#'
#' Given a center and covariance, get (x,y)-coordinates of points, which have
#' mahalanobis distance 1.
#'
#' @param S A 2x2 covariance matrix
#' @param mu A 2-dimensional vector for the ellipse's center
#' @param n number of points on the ellipse
#'
#' @examples
#'     S <- toeplitz(c(1, 0.3, 0.7))
#'     x <- mvtnorm::rmvt(100, S, df = 3)
#'
#'     ## Plot variables 1 and 3
#'     idx <- c(1,3)
#'     plot(x[, idx])
#'     ## Plot projection of true covariance matrix
#'     lines(shapeNA:::ellipse(S[idx, idx], c(0, 0)), col = 2)
ellipse <- function(S, mu, n = 250) {
  t <- seq(0, 2*pi, length.out = n)
  circle <- matrix(c(cos(t), sin(t)), ncol = 2)
  m <- sqrt(mahalanobis(circle, mu, S))
  return(sweep(circle, 1, m, FUN = '/'))
}
