#' Fast implementation for Kernel Ridge Regression.
#' 
#' Constructs a learner for the divide and conquer version of KRR.
#'
#' This function is to be used with the CVST package as a dropin
#' replacement for \code{\link[CVST]{constructKRRLearner}}. The
#' implementation approximates the inversion of the kernel Matrix
#' using the divide an conquer scheme, lowering computational and
#' memory complexity from O(n^3) and O(n^2) to O(n^3/m^2) and
#' O(n^2/m^2) respectively, where m are the number of blocks to be
#' used (parameter nblocks). In theory safe values for m are <
#' n^(1/3), the function will issue a warning, if the value for m is
#' too large.
#'
#' 
#'
#' @return Returns a learner of type \code{\link[CVST]{CVST.learner}}
#'     suitable for \code{\link[CVST]{CV}} and
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
#' @export
constructFastKRRLearner <- function () {
    if(!requireNamespace("CVST")) stop("require the 'CVST' package")
    if(!requireNamespace("kernlab")) stop("require 'kernlab' package")

    learn.krr <- function(data, params) {
        stopifnot(CVST::isRegression(data))
        nblocks <- params$nblocks
        nobs <- nrow(data$x)
        if (log(nblocks) / log(nobs) > 1/3) {
            warning(
                "Number of blocks too large wrt. number of observations, log(m)/log(N) ",
                sprintf("%.2f", log(nblocks)), "/", sprintf("%.2f", log(nobs)),
                " = ", sprintf("%.2f", log(nblocks)/log(nobs)),
                ", should be < 1/3, you results may suffer numerical inaccurracy. ",
                "For detail see Zhang et. al. (2013)"
            )
        }

        ## make kernel function
        kpar <- params[setdiff(names(params), c("kernel", "lambda", "nblocks"))]
        kernel <- get_kernel_fun(params$kernel, kpar)

        ## make blocks for samples
        shuff <- sample(1:nobs)
        blocksizes <- makeBlocks(nobs, nblocks)
        bends <- cumsum(blocksizes)
        bstarts <- c(1, bends[-nblocks] + 1)

        ## we make nblock models for ththe subsamples
        ## this can be parallelized:
        models <- list()
        for(i in 1:nblocks) {
            iIndices <- shuff[ bstarts[i]:bends[i] ]
            models[[i]] <- CVST:::.krr(data$x[iIndices,], kernel, data$y[iIndices], params$lambda)
        }

        return(models)
    } # end learn.krr
    
    predict.krr <- function(models, newData) {
        stopifnot(CVST::isRegression(newData))

        nModels <- length(models)
        pred <- rep(0, nrow(newData$x))
        for(i in 1:nModels) {
            pred <- pred + CVST:::.krr.predict(newData$x, models[[i]])
        }
        pred <- pred / nModels
        return(pred)
    }
        
               
  return(CVST::constructLearner(learn.krr, predict.krr))
}

#' Fast implementation of Kernel Ridge Regression
#'
#' To be used as a dropin replacement of .krr to be used with constructFastKRRLearner
#'
#' implements a divide and conquer scheme that is easier on large datasets
#'  
#' @param data output of CVST::constructData
#' @param kernel kernel function or string, see
#'   \code{\link[kernlab]{kernel-klass}}
## .fastkrr <- function (data, kernel, y, lambda, nblocks) {
##     nobs <- nrow(data)

##     if(!requireNamespace("kernlab"))
##         stop("require package 'kernlab'")

    
##     shuff <- sample(1:nobs)

##     blocksizes <- makeBlocks(nobs, nblocks)
##     bends <- cumsum(blocksizes)
##     bstarts <- c(1, bends[-nblocks] + 1)

##     ## this can be parallelized:
##     alpha <- numeric(nobs)
##     for(i in 1:nblocks) {
##         iIndices <- shuff[ bstarts[i]:bends[i] ]
##         iK <- kernlab::kernelMatrix(kernel, data[iIndices,])
##         iN <- length(iIndices)

##         alpha[iIndices] <- solve(iK + diag(lambda, iN), y[iIndices])
##     }

##     return(list(data = data, kernel = kernel, alpha = alpha))
## }

#' dividing nobs into nblocks parts that have approximately the same size
#'
#' @param nobs total number of observations
#' @param nblocks number of blocks
#'
#' @return vector of integers of length \code{nblocks} that sums up to
#'   \code{nobs}
#' 
makeBlocks <- function (nobs, nblocks) {
 
  maxbs <- nobs %/% nblocks
  rest <-  nobs %% nblocks

  res <- rep(maxbs, nblocks)
  if(rest > 0) res[1:rest] <- maxbs + 1
  return(res)
}


## get the kernel function out of the kernlab namespace:
get_kernel_fun <- function (kernel, pars) {
    if (!is(kernel,"kernel")) {
        if (is(kernel,"function")) {
            kernel <- deparse(substitute(kernel))
        } else {
            kernel <- get(kernel, asNamespace('kernlab'))
        }
        kernel <- do.call(kernel, pars)
    }
    return(kernel)
}



