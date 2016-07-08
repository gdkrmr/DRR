#' Dimensionality Reduction via Regression
#'
#' \code{DRR} implements the Dimensionality Reduction via Regression using
#' Kernel Ridge Regression.
#'
#' Parameter combination will be formed and crossvalidation used to
#' select the best combination. Crossvalidation uses
#' \code{\link[CVST]{CV}} or \code{\link[CVST]{fastCV}}.
#'
#' Pre-treatment of the data using a PCA and scaling is made
#' \eqn{\alpha = Vx}.  the representation in reduced dimensions is
#' \deqn{y_i = \alpha - f_i(\alpha_1, \ldots, \alpha_{i-1})} then the
#' final DRR representation is: \deqn{r = (\alpha_1, y_2, y_3, \ldots,
#' y_d)}.  DRR is invertible by \deqn{\alpha_i = y_i + f_i(\alpha_1,
#' \alpha_2, \ldots, alpha_{i-1})} If less dimensions are estimated,
#' there will be less inverse functions and calculating the inverse
#' will be inaccurate.
#' 
#' 
#' @references
#' Laparra, V., Malo, J., Camps-Valls, G., 2015. Dimensionality
#'     Reduction via Regression in Hyperspectral Imagery. IEEE Journal
#'     of Selected Topics in Signal Processing 9,
#'     1026–1036. doi:10.1109/JSTSP.2015.2417833
#'
#' @param X input data, a matrix.
#' @param ndim the number of output dimensions and regression
#'     functions to be estimated, see details for inversion.
#' @param lambda the penalty term for the Kernel Ridge Regression.
#' @param kernel a kernel function or string, see
#'     \code{\link[kernlab]{kernel}} for details.
#' @param kernel.pars a list with parameters for the kernel. each
#'     parameter can be a vector, crossvalidation will choose the best
#'     combination.
#' @param pca logical, do a preprocessing using pca.
#' @param pca.center logical, center data before applying pca.
#' @param pca.scale logical, scale data before applying pca.
#' @param fastcv if \code{TRUE} uses \code{\link[CVTS]{fastCV}}, if
#'     \code{FALSE} uses \code{\link[CVST]{CV}} for crossvalidation.
#' @param cv.folds if using normal crossvalidation, the number of
#'     folds to be used.
#' @param fastcv.test an optional separate test data set to be used
#'     for \code{\link[CVTS]{fastCV}}, handed over as option
#'     \code{test} to \code{\link[CVTS]{fastCV}}.
#' @param fastkrr.nblocks the number of blocks used for fast KRR,
#'     higher numbers are faster to compute but may introduce
#'     numerical inaccurracies, see
#'     \code{\link{constructFastKRRLearner}} for details.
#' @param verbose logical, should the crossvalidation report back.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item{fitted.data}{The data in reduced dimensions.}
#' \item{pca.means}{The means used to center the original data.}
#' \item{pca.scale}{The standard deviations used to scale the original data.}
#' \item{pca.rotation}{The rotation matrix of the PCA.}
#' \item{models}{A list of models used to estimate each dimension.}
#' \item{apply}{A function to fit new data to the estimated model.}
#' \item(inverse}{A function to untransform data.}
#' }
#' 
#'
#' tt <- seq(0,4*pi, length.out = 2000)
#' helix <- cbind(
#'   x = 3 * cos(tt) + rnorm(length(tt), sd = seq(0.1, 1.4, length.out = length(tt))),
#'   y = 3 * sin(tt) + rnorm(length(tt), sd = seq(0.1, 1.4, length.out = length(tt))),
#'   z = 2 * tt      + rnorm(length(tt), sd = seq(0.1, 1.4, length.out = length(tt)))
#' )
#'
#' library(profvis)
#' p <- profvis({
#'   drr <- DRR(helix, ndim = 3, cv.folds = 5,
#'              fastkrr.nblocks = 4, verbose = TRUE,
#'              fastcv = FALSE)
#' })
#' print(p)
#'
#' 
#' \dontrun{
#' library(rgl)
#' plot3d()
#' plot3d(helix)
#' points3d(drr$inverse(drr$fitted.data[,1,drop = FALSE]), col = 'blue')
#'
#' 
#' plot3d(drr$fitted.data)
#' xx <- seq(-0, 250, length.out = 25)
#' yy <- seq(-40, 40, length.out = 5)
#' zz <- seq(-40, 40, length.out = 5)
#'
#' dd <- as.matrix(expand.grid(xx, yy, zz))
#' plot3d(helix)
#' for(y in yy) for(x in xx)
#'                  rgl.linestrips(drr$inverse(cbind(x, y, zz)), col = 'blue')
#' for(y in yy) for(z in zz)
#'                  rgl.linestrips(drr$inverse(cbind(xx, y, z)), col = 'blue')
#' for(x in xx) for(z in zz)
#'                  rgl.linestrips(drr$inverse(cbind(x, yy, z)), col = 'blue')
#' 
#' 
#' @export
DRR <- function (X, ndim = ncol(X),
                 lambda = c(0, 10^(-2:1)),
                 kernel      = 'rbfdot',
                 kernel.pars = list(sigma = 10^(-2:3)),
                 pca        = TRUE,
                 pca.center = TRUE,
                 pca.scale  = TRUE,
                 fastcv      = FALSE,
                 cv.folds    = 5,
                 fastcv.test = NULL,
                 fastkrr.nblocks = 4,
                 verbose  = TRUE)  {
    if((!fastcv) && (cv.folds <= 1)) stop("need more than one fold for crossvalidation")
    if(cv.folds %% 1 != 0) stop("cv.folds must be a whole number")
    if(fastkrr.nblocks < 1) stop("fastkrr.nblocks must be at least 1")
    if(fastkrr.nblocks %% 1 != 0) stop('fastkrr.nblocks must be a whole number')
    if(!requireNamespace("CVST")) stop("require the 'CVST' package")
    if(!requireNamespace("kernlab")) stop("require 'kernlab' package")
    if(ndim < ncol(X)) warning('ndim < data dimensionality, the inverse functions will be incomplete!')
    if(ndim > ncol(X)) ndim <- ncol(X)

    if(pca) {
        pca <- prcomp(X, center = pca.center, scale. = pca.scale)
        if (!pca.center) pca$center <- rep(0, ncol(X))
        if (!pca.scale) pca$scale   <- rep(1, ncol(X))
    } else {
        pca <- list()
        pca$x <- X
        pca$rotation <- diag(1, ncol(X), ncol(X))
        pca$center <- rep(0, ncol(X))
        pca$scale <- rep(1, ncol(X))
    }

    
    alpha <- pca$x
    d <- ndim

    kpars <- kernel.pars
    kpars$kernel <- kernel
    kpars$lambda <- lambda
    kpars$nblocks <- fastkrr.nblocks

    krrl <- constructFastKRRLearner()
       
    p <- do.call(CVST::constructParams, kpars)
    
    Y <- matrix(NA_real_, nrow = nrow(X), ncol = d)
    models <- list()
    if(d > 1) for (i in d:2) {
        message(Sys.time(), ": Constructing Axis ", d-i+1, "/", d)
        data <- CVST::constructData(
            x = alpha[,1:(i-1), drop = FALSE],
            y = alpha[,i]
        )

        cat("predictors: ", colnames(alpha)[1:(i-1)], "dependent: ", colnames(alpha)[i], '\n')

     
        res <- if (fastcv) {
                   CVST::fastCV(
                       data, krrl, p,
                       CVST::constructCVSTModel(),
                       test = fastcv.test,
                       verbose = verbose
                   )
               } else {
                   CVST::CV(
                       data, krrl, p,
                       fold = cv.folds,
                       verbose = verbose
                   )
               }
        model <- krrl$learn(data, res[[1]])

        models[[i]] <- model
        Y[,i] <- alpha[,i] - krrl$predict(model, data)
    }
    ## we don't need to construct the very last dimension
    message(Sys.time(), ": Constructing Axis ", d, "/", d)
    Y[,1] <- alpha[,1]
    models[[1]] <- list()

    
    appl <- function(x) {
        ## apply PCA
        dat <- scale(x, pca$center, pca$scale)
        dat <-  dat %*% pca$rotation
        
        ## apply KRR
        outdat <- matrix(NA_real_, ncol = d, nrow = nrow(x))
        if(d > 1) for (i in d:2)
            outdat[,i] <-
                dat[,i] - krrl$predict(
                              models[[i]],
                              CVST:::constructData(x = dat[,1:(i-1), drop = FALSE],
                                                   y = NA)
                          )
        
        outdat[,1] <- dat[,1]
    
        return(outdat)                
    }

    inv <- function(x){
        #browser()
        dat <- cbind(x, matrix(0, nrow(x), ncol(X)-ncol(x)))

        outdat <- dat #matrix(NA_real_, nrow(x), ncol(X))

        ## krr
        #outdat[,1] <- dat[,1]
        if(d > 1) for (i in 2:d)
            outdat[,i] <- dat[,i] + krrl$predict(
                              models[[i]],
                              CVST:::constructData(x = outdat[,1:(i-1), drop = FALSE],
                                                   y = NA)
                          )
        
        ## inverse pca
        outdat <- outdat %*% t(pca$rotation)
        outdat <- sweep(outdat, 2L, pca$scale, "*")
        outdat <- sweep(outdat, 2L, pca$center, "+")
        
        return(outdat)        
    }
                                   
    return(list(
        fitted.data = Y,
        pca.means = pca$center,
        pca.scale = pca$scale,
        pca.rotation = pca$rotation,
        models = models,
        apply = appl,
        inverse = inv
    ))
}

