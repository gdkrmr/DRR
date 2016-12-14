#' Fast implementation for Kernel Ridge Regression.
#'
#' Constructs a learner for the divide and conquer version of KRR.
#'
#' This function is to be used with the CVST package as a drop in
#' replacement for \code{\link[CVST]{constructKRRLearner}}. The
#' implementation approximates the inversion of the kernel Matrix
#' using the divide an conquer scheme, lowering computational and
#' memory complexity from \eqn{O(n^3)} and \eqn{O(n^2)} to
#' \eqn{O(n^3/m^2)} and \eqn{O(n^2/m^2)} respectively, where m are the
#' number of blocks to be used (parameter nblocks). Theoretically safe
#' values for \eqn{m} are \eqn{< n^{1/3}}, but practically \eqn{m} may
#' be a little bit larger. The function will issue a warning, if the
#' value for \eqn{m} is too large.
#'
#'
#'
#' @return Returns a learner similar to \code{\link[CVST]{constructKRRLearner}}
#'     suitable for the use with \code{\link[CVST]{CV}} and
#'     \code{\link[CVST]{fastCV}}.
#'
#' @seealso \code{\link[CVST]{constructLearner}}
#'
#' @references
#' Zhang, Y., Duchi, J.C., Wainwright, M.J., 2013. Divide and Conquer
#'     Kernel Ridge Regression: A Distributed Algorithm with Minimax
#'     Optimal Rates. arXiv:1305.5029 [cs, math, stat].
#'
#' @examples
#' ns <- noisySinc(1000)
#' nsTest <- noisySinc(1000)
#'
#' fast.krr <- constructFastKRRLearner()
#' fast.p <- list(kernel="rbfdot", sigma=100, lambda=.1/getN(ns), nblocks = 4)
#' system.time(fast.m <- fast.krr$learn(ns, fast.p))
#' fast.pred <- fast.krr$predict(fast.m, nsTest)
#' sum((fast.pred - nsTest$y)^2) / getN(nsTest)
#'
#' \dontrun{
#' krr <- CVST::constructKRRLearner()
#' p <- list(kernel="rbfdot", sigma=100, lambda=.1/getN(ns))
#' system.time(m <- krr$learn(ns, p))
#' pred <- krr$predict(m, nsTest)
#' sum((pred - nsTest$y)^2) / getN(nsTest)
#'
#' plot(ns, col = '#00000030', pch = 19)
#' lines(sort(nsTest$x), fast.pred[order(nsTest$x)], col = '#00C000', lty = 2)
#' lines(sort(nsTest$x), pred[order(nsTest$x)], col = '#0000C0', lty = 2)
#' legend('topleft', legend = c('fast KRR', 'KRR'),
#'        col = c('#00C000', '#0000C0'), lty = 2)
#' }
#'
#' @import Matrix
#' @import CVST
#' @import kernlab
#' @export
constructFastKRRLearner <- function() { # nolint
    if (!requireNamespace("CVST")) stop("require the 'CVST' package")
    if (!requireNamespace("kernlab")) stop("require 'kernlab' package")

    learn.krr <- function(data, params) {
        stopifnot(CVST::isRegression(data))
        nblocks <- params$nblocks
        nobs <- nrow(data$x)
        if (log(nblocks) / log(nobs) > 1 / 3) {
            warning(
                "Number of blocks too large wrt. number of observations,",
                " log(m)/log(N) = ",
                sprintf("%.2f", log(nblocks)), "/", sprintf("%.2f", log(nobs)),
                " = ", sprintf("%.2f", log(nblocks) / log(nobs)),
                ", should be < 1/3, ",
                "you results may suffer numerical inaccurracy. ",
                "For detail see Zhang et. al. (2013)"
            )
        }

        ## make kernel function
        kpar <- params[setdiff(names(params), c("kernel", "lambda", "nblocks"))]
        kernel <- get_kernel_fun(params$kernel, kpar)

        ## make blocks for samples
        shuff <- sample(1:nobs)
        blocksizes <- make_blocks(nobs, nblocks)
        bends <- cumsum(blocksizes)
        bstarts <- c(1, bends[-nblocks] + 1)

        ## we make nblock models for the subsamples
        ## this can be parallelized:
        models <- list()
        for (i in 1:nblocks) {
            i_indices <- shuff[ bstarts[i]:bends[i] ]
            models[[i]] <- krr(data$x[i_indices, ],
                               kernel,
                               data$y[i_indices],
                               params$lambda)
        }
        return(models)
    } # end learn.krr

    ## newData cannot be renamed because it is the convention of the
    ## CVST package
    predict.krr <- function(models, newData) {
        stopifnot(CVST::isRegression(newData))

        n_models <- length(models)
        pred <- rep(0, nrow(newData$x))
        for (i in 1:n_models) {
            pred <- pred + krr.predict(newData$x, models[[i]])
        }
        pred <- pred / n_models
        return(as.matrix(pred))
    }


  return(CVST::constructLearner(learn.krr, predict.krr))
}

# dividing nobs into nblocks parts that have approximately the same size
#
# @param nobs total number of observations
# @param nblocks number of blocks
#
# @return vector of integers of length \code{nblocks} that sums up to
#   \code{nobs}
#
make_blocks <- function(nobs, nblocks) {

  maxbs <- nobs %/% nblocks
  rest <-  nobs %% nblocks

  res <- rep(maxbs, nblocks)
  if (rest > 0) res[1:rest] <- maxbs + 1
  return(res)
}


## get the kernel function out of the kernlab namespace:
get_kernel_fun <- function (kernel, pars) {
    if (!methods::is(kernel, "kernel")) {
        if (methods::is(kernel, "function")) {
            kernel <- deparse(substitute(kernel))
        } else {
            kernel <- get(kernel, asNamespace("kernlab"))
        }
        kernel <- do.call(kernel, pars)
    }
    return(kernel)
}


## internal functions from cvst package, have to be here, because CVST
## does not export them and CRAN does not allow the use of unexported
## functions.

## CVST:::.krr
krr <- function (data, kernel, y, lambda) {
    K <- kernlab::kernelMatrix(kernel, data)
    N <- nrow(K)
    alpha <- solve(Matrix(K + diag(lambda, N))) %*% y
    return(list(data = data, kernel = kernel, alpha = alpha))
}

## CVST:::.krr.predict
krr.predict <- function (new_data, krr) {
    k <- kernlab::kernelMatrix(krr$kernel, new_data, krr$data)
    return(k %*% krr$alpha)
}
