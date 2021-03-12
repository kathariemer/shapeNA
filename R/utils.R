#' Reorder Data with Missing Values
#'
#' Reorder a data set with `NA` entries to form blocks of missing values. The
#' resulting data will have increasing missingness along the rows and along the
#' columns. The rows are ordered such that the first block consists of complete
#' observations, and the following blocks are ordered from most frequent
#' missingness pattern to least frequent missingness pattern.
#'
#' In case of ties, that is if two patterns occur with the same frequency, the block
#' whose pattern occurs first will be ordered in front of the other block.
#'
#' This method may fail if the missingness is too strong or if the number of
#' observations is too low (the number of observations has to exceed the number
#' of variables), as it has been designed as a preprocessing step for shape
#' estimations.
#'
#' @param x A matrix with missing values.
#' @param cleanup A logical flag. If `TRUE`, observations with less than 2 responses
#'    are discarded.
#' @param plot A logical flag. If `TRUE`, a plot of the missingness pattern is produced.
#'
#' @return A list of class `naBlocks` with components:
#'     \item{x}{The reordered data matrix.}
#'     \item{permutation}{The permutation of the columns that was applied to reorder the columns according to the number of `NA`s.}
#'     \item{rowPermutation}{The permutation of the rows that generates the blocks.}
#'     \item{N}{A vector of all row indices. Each row number points to the beginning of a new missingness pattern.}
#'     \item{D}{A vector specifying the missingness pattern for each block.}
#'     \item{P}{A vector specifying the number of observed variables per block.}
#'     \item{kn}{A vector specifying the percentage of observed responses per variable.}
naBlocks <- function(x, cleanup=TRUE, plot=FALSE) {
  ## remove observations without responses or only 1 response ##
  if (cleanup) {
    hasValue <- !is.na(x) # R matrix
    #rmIdx <- apply(hasValue, 1, function(x) sum(x) <=1)
    rmIdx <- rowSums(hasValue)<=1
    if (any(rmIdx)) {
      x <- x[!rmIdx, ]
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
  # columns reordered from most to least observations
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

  # rows reordered to form blocks
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


#' Plot Missingness Pattern of Data
#'
#' Function to visualize the missingness patterns for objects of class `naBlocks`.
#'
#' @param x A `naBlocks` object.
#' @param ... Additional parameters passed on to \code{\link[graphics]{rect}}.
#'
#' @return No return value.
#'
#' @export
#' @examples
#'     x <- mvtnorm::rmvt(100, toeplitz(seq(1, 0.1, length.out = 3)), df = 5)
#'     y <- mice::ampute(x, mech='MCAR')$amp
#'     res <- classicShapeNA(y)
#'     plot(res$naBlocks)
plot.naBlocks <- function(x, ...) {
  idx <- x$N
  bprop <- (idx - c(0, idx[-length(idx)]))/nrow(x$data)
  bIdx <- 1:length(idx)
  plotAsBoolMat(x$data[idx,], round(bprop,3), ...)
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

#' Matrix Plot for Binary Matrices
#'
#' If the matrix contains missing values, instead of the matrix itself, a boolean
#' matrix representing the missingness is plotted, with blue cells representing
#' observed and yellow cells representing missing values.
#' Otherwise, cells with entries 0 are colored white and all other cells are colored
#' blue.
#'
#' @param R Matrix with missing values or with entries \code{{0, 1}}.
#' @param blockNames An optional vector with length \code{nrow(R)}, which is
#'     displayed to the left of the plot.
#' @param ... Additional arguments passed on to \code{\link[graphics]{rect}}.
#'
#' @noRd
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
  graphics::plot(1, type='n', xlim=c(0,c+1), ylim = c(r+1,0), axes=FALSE, xlab=NA, ylab=NA)
  xleft <- rep(1, r)
  ybottom <- 0:(r-1)
  for (i in 1:c) {
    rowColor <- ifelse(B[,i], color2, colors["blue"])
    graphics::rect(xleft, ybottom, xleft+1, ybottom+1, col = rowColor, ...)
    xleft <- xleft+1
  }
  graphics::text(0.5+(1:c),r+1, colnames(R))
  if (!missing(blockNames)){
    graphics::text(0, (1:r)-0.5, blockNames)
  }
}

# transform numeric input to have values between 0 and 1
# the largest value gets mapped to 0 (black), the smallest to 1 (white)
grayscale <- function(A) {
  m <- min(A)
  M <- max(A)
  return((A-m)/(m-M)+1)
}

#' Print Method for Objects of Class `shapeNA`
#'
#' Prints the chosen value of `alpha` as well as the estimated shape and
#' location for objects of class `shapeNA`.
#'
#' @param x A `shapeNA` object
#' @param ... Additional parameters passed to lower level \code{\link[base]{print}}.
#'
#' @return No return value.
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
  invisible(NULL)
}

#' Visualization of Shape Estimate
#'
#' Function to visualize the shape matrix from objects of class `shapeNA` by
#' plotting a heatmap where light colored cells indicate small values and dark
#' colored cells indicate high values.
#'
#' @param x A `shapeNA` oopbject
#' @param message A logical, If `TRUE`, the percentage of observed values per
#'     variable is printed in the console.
#' @param ... Additional parameters passed to \code{\link[graphics]{image}}.
#'
#' @return A matrix with the proportion of observed values for each variable.
#'
#' @export
#' @examples
#'     x <- mvtnorm::rmvt(100, toeplitz(seq(1, 0.1, length.out = 3)), df = 5)
#'     y <- mice::ampute(x, mech='MCAR')$amp
#'     res <- tylerShapeNA(y)
#'     ## default plot
#'     plot(res)
#'     ## plot result in gray scale - reverse order to get a palette starting
#'     ## with the lightest instead of the darkest color
#'     plot(res, col = gray.colors(9, rev = TRUE))
plot.shapeNA <- function(x, message = TRUE, ...) {
  imageS <- x$S
  p <- ncol(imageS)
  vars <- x$naBlocks$colNames
  ## "flip" image to have 1st entry in lower left corner
  imageS <- t(imageS[p:1,])
  graphics::image(imageS, xaxt = 'n', yaxt = 'n', ...)
  graphics::axis(side = 1, at = seq(0, 1, length.out = p),
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

#' Get Row Indices which Produce a Block Structure of Missingness
#'
#' Use a mapping of binary vectors to integers, that is
#' \deqn{(x_1, x_2, \ldots, x_p) \mapsto \sum_{i=1}^p 2^{p-i}x_i}
#' to handle the missingness for each record in a more compact way. Return a
#' row order, such that the missingness forms blocks in the matrix `R` and also
#' return a run length encoding of the integer vector, which shows for each
#' pattern the number of consecutive rows, which have the same missingness
#' structure.
#'
#' @param R A boolean matrix.
#'
#' @return A list with
#'     \item{order}{A vector, which subsets the input matrix R in such a way,
#'     that the boolean values form blocks.}
#'     \item{rle}{An `rle` object, whose `values` give the blocks pattern as
#'     integers and whose `lengths` give the size of each block.}
#' @noRd
getMissingnessBlocks <- function(R) {
  p2 <- 2^((ncol(R)-1):0)
  # encode missingness as integer
  pattern <- rowSums(sweep(R, 2, p2, "*"))
  # named vector with pattern named with corresponding integer and it's
  # frequency as value
  lev <- summary(as.factor(pattern))
  patternF <- factor(pattern, levels = names(sort(lev, decreasing = TRUE)))
  rowOrder <- order(patternF)
  return(list(order=rowOrder, rle = rle(pattern[rowOrder])))
}

#' Print Missingness Pattern
#'
#' Print the pattern of missingness in the supplied data, as a block matrix.
#' Observed data are represented by 1, missing values by 0.
#for 1s, representing a column vector of responses and 0s, representing a column vector of `NA`s.
#'
#' The first row shows the column names. The leftmost column, without column
#' name, shows the number of rows per block and the rightmost column with name
#' `#` shows the number of observed variables in the block.
#'
#' @param x An `naBlocks` object.
#' @param ... Additional parameters passed to \code{\link[base]{print}}.
#'
#' @return A named matrix representing the missingness pattern of the data.
#'
#' @export
#' @examples
#'     x <- mvtnorm::rmvt(100, toeplitz(seq(1, 0.1, length.out = 3)), df = 5)
#'     y <- mice::ampute(x, mech='MCAR')$amp
#'     res <- classicShapeNA(y)
#'     print(res$naBlocks)
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
#' @param height A `shapeNA` object.
#' @param sortNA A logical. If `FALSE`, the original variable order is kept.
#'     Otherwise the variables are ordered from least to most missingness.
#' @param ... Additional graphical arguments passed to
#'     \code{\link[graphics]{barplot}}.
#'
#' @return Invisibly returns a named vector holding the proportion of
#' missingness per variable.
#'
#' @importFrom graphics barplot
#'
#' @seealso \code{\link[graphics]{barplot}}
#'
#' @export
#' @examples
#'     S <- toeplitz(seq(1, 0.1, length.out = 3))
#'     x <- mvtnorm::rmvt(100, S, df = 5)
#'     y <- mice::ampute(x, mech='MCAR')$amp
#'     res <- classicShapeNA(y)
#'     barplot(res)
barplot.shapeNA <- function(height, sortNA = FALSE, ...) {
  if (is.null(height$naBlocks)) {
    stop('No missing values in the data.')
  }
  blocks <- height$naBlocks
  mprop <- colSums(is.na(blocks$data))/nrow(blocks$data)
  if (!sortNA) {
    mprop <- mprop[order(blocks$permutation)]
  }
  barplot(mprop, ..., main = 'Proportion of missing values')
  invisible(mprop)
}

#' Summary Method for Class `shapeNA`
#'
#' Summary methods for objects from class `shapeNA`.
#'
#' @param object An object of class `shapeNA`, usually from a call to
#'     \code{\link{powerShape}} or similar functions.
#' @param ... Further arguments to be passed to or from methods.
#'
#' @return A `summary.shapeNA` object. For objects of this class, the `print`
#' method tries to format the location and shape estimate in a readable format
#' and also shows the number of iterations, before the algorithm converged.
#'
#' @export
#'
#' @examples
#'     obj <- tylerShape(mvtnorm::rmvt(100, diag(3)))
#'     summary(obj)
summary.shapeNA <- function(object, ...) {
  class(object) <- 'summary.shapeNA'
  return(object)
}

#' Print Method for Class `summary.shapeNA`
#'
#' @param x Object returned from \code{\link{summary.shapeNA}}.
#' @param ... Further arguments to be passed to or from methods.
# @param ... Further arguments, which will be ignored.
#'
#' @return No return value.
#'
#' @export
#'
#' @examples
#'     obj <- tylerShape(mvtnorm::rmvt(100, diag(3)))
#'     print(summary(obj))
print.summary.shapeNA <- function(x, ...) {
  cat('Call: ')
  print(x$call)
  cat('---\nAlpha:', x$alpha)
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

#' Scatter Estimates from `shapeNA` Objects
#'
#' For Power M-estimates with tail index `alpha < 1`, the resulting estimate
#' has a scale. For these shape estimates, scatter matrices can be computed.
# cases, a scatter estimate can be computed.
#' Results from
#' \code{\link{tylerShape}} and \code{\link{tylerShapeNA}} give no scatter
#' estimates. In these cases the function returns `NA`.
#'
#' @param obj `shapeNA` object, resulting from a call to
#'  \code{\link{powerShape}} and other functions from the same family.
#'
#' @return Scatter matrix estimate, or only `NA` if `alpha` = 1.
#'
#' @export
#'
#' @examples
#'     S <- toeplitz(c(1, 0.3, 0.7))
#'     set.seed(123)
#'     x <- mvtnorm::rmvt(100, S, df = 3)
#'     obj_det <- powerShape(x, alpha = 0.85, normalization = 'det')
#'     shape2scatter(obj_det)
#'     obj_tr <- powerShape(x, alpha = 0.85, normalization = 'trace')
#'     shape2scatter(obj_tr)
#'     obj_one <- powerShape(x, alpha = 0.85, normalization = 'one')
#'     shape2scatter(obj_one)
shape2scatter <- function(obj) {
  if (obj$alpha == 1) {
    message("Covariance matrix for Tyler's M-estimate undefined")
    return(NA)
  }
  return(obj$scale * obj$S)
}
