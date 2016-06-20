
#' DRR
#'
#' Dimensionality reduction via regression
#'
#' Currently uses linear regression, so output should be a PCA!!!
#'
#' @param X input data, a matrix.
#'
#' @return The data in reduced dimensions.
#' 
DRR <- function (X)  {
  require(CVST)
  alpha <- prcomp(X)$x
  d <- ncol(X)

  Y <- matrix(NA_real_, nrow = nrow(X), ncol = d)
  for (i in d:2) {
    data <- constructData(
      x = alpha[,1:(i-1), drop = FALSE],
      y = alpha[,i]
    )
    krr <- constructKRRLearner()
    p   <- constructParams(kernel = "rbfdot",
                           sigma=10^(-1:5),
                           lambda=10^(-2:3))
    invisible(res <- fastCV(data, krr, p, constructCVSTModel()))
    alpha_hat_i <- krr$learn(data, res[[1]])
    ## alpha_hat_i <- lm.fit(alpha[,1:(i-1), drop = FALSE],
    ##                       alpha[,i, drop = FALSE])$fitted.values
    Y[,d-i+1] <- alpha[,i] - alpha_hat_i$alpha[,1]
  }
  Y[,d] <- alpha[,d]
  return(Y)
}
