# Dimensionality Reduction via Regression

`drr` Implements Dimensionality Reduction via Regression using Kernel
Ridge Regression.

## Usage

``` r
drr(
  X,
  ndim = ncol(X),
  lambda = c(0, 10^(-3:2)),
  kernel = "rbfdot",
  kernel.pars = list(sigma = 10^(-3:4)),
  pca = TRUE,
  pca.center = TRUE,
  pca.scale = FALSE,
  fastcv = FALSE,
  cv.folds = 5,
  fastcv.test = NULL,
  fastkrr.nblocks = 4,
  verbose = TRUE
)
```

## Arguments

- X:

  input data, a matrix.

- ndim:

  the number of output dimensions and regression functions to be
  estimated, see details for inversion.

- lambda:

  the penalty term for the Kernel Ridge Regression.

- kernel:

  a kernel function or string, see
  [`kernel-class`](https://rdrr.io/pkg/kernlab/man/kernel-class.html)
  for details.

- kernel.pars:

  a list with parameters for the kernel. each parameter can be a vector,
  crossvalidation will choose the best combination.

- pca:

  logical, do a preprocessing using pca.

- pca.center:

  logical, center data before applying pca.

- pca.scale:

  logical, scale data before applying pca.

- fastcv:

  if `TRUE` uses [`fastCV`](https://rdrr.io/pkg/CVST/man/fastCV.html),
  if `FALSE` uses [`CV`](https://rdrr.io/pkg/CVST/man/CV.html) for
  crossvalidation.

- cv.folds:

  if using normal crossvalidation, the number of folds to be used.

- fastcv.test:

  an optional separate test data set to be used for
  [`fastCV`](https://rdrr.io/pkg/CVST/man/fastCV.html), handed over as
  option `test` to [`fastCV`](https://rdrr.io/pkg/CVST/man/fastCV.html).

- fastkrr.nblocks:

  the number of blocks used for fast KRR, higher numbers are faster to
  compute but may introduce numerical inaccurracies, see
  [`constructFastKRRLearner`](constructFastKRRLearner.md) for details.

- verbose:

  logical, should the crossvalidation report back.

## Value

A list the following items:

- "fitted.data" The data in reduced dimensions.

- "pca.means" The means used to center the original data.

- "pca.scale" The standard deviations used to scale the original data.

- "pca.rotation" The rotation matrix of the PCA.

- "models" A list of models used to estimate each dimension.

- "apply" A function to fit new data to the estimated model.

- "inverse" A function to untransform data.

## Details

Parameter combination will be formed and cross-validation used to select
the best combination. Cross-validation uses
[`CV`](https://rdrr.io/pkg/CVST/man/CV.html) or
[`fastCV`](https://rdrr.io/pkg/CVST/man/fastCV.html).

Pre-treatment of the data using a PCA and scaling is made \\\alpha =
Vx\\. the representation in reduced dimensions is

\$\$y_i = \alpha - f_i(\alpha_1, \ldots, \alpha\_{i-1})\$\$

then the final DRR representation is:

\$\$r = (\alpha_1, y_2, y_3, \ldots,y_d)\$\$

DRR is invertible by

\$\$\alpha_i = y_i + f_i(\alpha_1,\alpha_2, \ldots, alpha\_{i-1})\$\$

If less dimensions are estimated, there will be less inverse functions
and calculating the inverse will be inaccurate.

## References

Laparra, V., Malo, J., Camps-Valls, G., 2015. Dimensionality Reduction
via Regression in Hyperspectral Imagery. IEEE Journal of Selected Topics
in Signal Processing 9, 1026-1036. doi:10.1109/JSTSP.2015.2417833

## Examples

``` r
tt <- seq(0,4*pi, length.out = 200)
helix <- cbind(
  x = 3 * cos(tt) + rnorm(length(tt), sd = seq(0.1, 1.4, length.out = length(tt))),
  y = 3 * sin(tt) + rnorm(length(tt), sd = seq(0.1, 1.4, length.out = length(tt))),
  z = 2 * tt      + rnorm(length(tt), sd = seq(0.1, 1.4, length.out = length(tt)))
)
helix <- helix[sample(nrow(helix)),] # shuffling data is important!!
system.time(
drr.fit  <- drr(helix, ndim = 3, cv.folds = 4,
                lambda = 10^(-2:1),
                kernel.pars = list(sigma = 10^(0:3)),
                fastkrr.nblocks = 2, verbose = TRUE,
                fastcv = FALSE)
)
#> 2026-05-15 15:15:10.289662: Constructing Axis 1/3
#> predictors:  PC1 PC2 dependent:  PC3 
#> sigma=1 kernel=rbfdot lambda=0.01 nblocks=2 ( 2.186397 )
#> sigma=10 kernel=rbfdot lambda=0.01 nblocks=2 ( 3.243588 )
#> sigma=100 kernel=rbfdot lambda=0.01 nblocks=2 ( 4.017735 )
#> sigma=1000 kernel=rbfdot lambda=0.01 nblocks=2 ( 4.176362 )
#> sigma=1 kernel=rbfdot lambda=0.1 nblocks=2 ( 2.025058 )
#> sigma=10 kernel=rbfdot lambda=0.1 nblocks=2 ( 3.308414 )
#> sigma=100 kernel=rbfdot lambda=0.1 nblocks=2 ( 4.029985 )
#> sigma=1000 kernel=rbfdot lambda=0.1 nblocks=2 ( 4.177476 )
#> sigma=1 kernel=rbfdot lambda=1 nblocks=2 ( 2.146495 )
#> sigma=10 kernel=rbfdot lambda=1 nblocks=2 ( 3.608597 )
#> sigma=100 kernel=rbfdot lambda=1 nblocks=2 ( 4.096955 )
#> sigma=1000 kernel=rbfdot lambda=1 nblocks=2 ( 4.183313 )
#> sigma=1 kernel=rbfdot lambda=10 nblocks=2 ( 3.460475 )
#> sigma=10 kernel=rbfdot lambda=10 nblocks=2 ( 4.062405 )
#> sigma=100 kernel=rbfdot lambda=10 nblocks=2 ( 4.172848 )
#> sigma=1000 kernel=rbfdot lambda=10 nblocks=2 ( 4.189501 )
#> 2026-05-15 15:15:11.013802: Constructing Axis 2/3
#> predictors:  PC1 dependent:  PC2 
#> sigma=1 kernel=rbfdot lambda=0.01 nblocks=2 ( 2.372495 )
#> sigma=10 kernel=rbfdot lambda=0.01 nblocks=2 ( 2.838575 )
#> sigma=100 kernel=rbfdot lambda=0.01 nblocks=2 ( 3.695045 )
#> sigma=1000 kernel=rbfdot lambda=0.01 nblocks=2 ( 4.726495 )
#> sigma=1 kernel=rbfdot lambda=0.1 nblocks=2 ( 1.911542 )
#> sigma=10 kernel=rbfdot lambda=0.1 nblocks=2 ( 2.451259 )
#> sigma=100 kernel=rbfdot lambda=0.1 nblocks=2 ( 3.786396 )
#> sigma=1000 kernel=rbfdot lambda=0.1 nblocks=2 ( 4.781136 )
#> sigma=1 kernel=rbfdot lambda=1 nblocks=2 ( 1.827136 )
#> sigma=10 kernel=rbfdot lambda=1 nblocks=2 ( 2.815618 )
#> sigma=100 kernel=rbfdot lambda=1 nblocks=2 ( 4.222793 )
#> sigma=1000 kernel=rbfdot lambda=1 nblocks=2 ( 4.987592 )
#> sigma=1 kernel=rbfdot lambda=10 nblocks=2 ( 3.51036 )
#> sigma=10 kernel=rbfdot lambda=10 nblocks=2 ( 4.515441 )
#> sigma=100 kernel=rbfdot lambda=10 nblocks=2 ( 5.086705 )
#> sigma=1000 kernel=rbfdot lambda=10 nblocks=2 ( 5.279499 )
#> 2026-05-15 15:15:11.619287: Constructing Axis 3/3
#>    user  system elapsed 
#>   1.662   2.011   1.334 

if (FALSE) { # \dontrun{
library(rgl)
plot3d(helix)
points3d(drr.fit$inverse(drr.fit$fitted.data[,1,drop = FALSE]), col = 'blue')
points3d(drr.fit$inverse(drr.fit$fitted.data[,1:2]),             col = 'red')

plot3d(drr.fit$fitted.data)
pad <- -3
fd <- drr.fit$fitted.data
xx <- seq(min(fd[,1]),       max(fd[,1]),       length.out = 25)
yy <- seq(min(fd[,2]) - pad, max(fd[,2]) + pad, length.out = 5)
zz <- seq(min(fd[,3]) - pad, max(fd[,3]) + pad, length.out = 5)

dd <- as.matrix(expand.grid(xx, yy, zz))
plot3d(helix)
for(y in yy) for(x in xx)
  rgl.linestrips(drr.fit$inverse(cbind(x, y, zz)), col = 'blue')
for(y in yy) for(z in zz)
  rgl.linestrips(drr.fit$inverse(cbind(xx, y, z)), col = 'blue')
for(x in xx) for(z in zz)
  rgl.linestrips(drr.fit$inverse(cbind(x, yy, z)), col = 'blue')
} # }
```
