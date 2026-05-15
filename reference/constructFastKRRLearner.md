# Fast implementation for Kernel Ridge Regression.

Constructs a learner for the divide and conquer version of KRR.

## Usage

``` r
constructFastKRRLearner()
```

## Value

Returns a learner similar to
[`constructKRRLearner`](https://rdrr.io/pkg/CVST/man/constructLearner.html)
suitable for the use with [`CV`](https://rdrr.io/pkg/CVST/man/CV.html)
and [`fastCV`](https://rdrr.io/pkg/CVST/man/fastCV.html).

## Details

This function is to be used with the CVST package as a drop in
replacement for
[`constructKRRLearner`](https://rdrr.io/pkg/CVST/man/constructLearner.html).
The implementation approximates the inversion of the kernel Matrix using
the divide an conquer scheme, lowering computational and memory
complexity from \\O(n^3)\\ and \\O(n^2)\\ to \\O(n^3/m^2)\\ and
\\O(n^2/m^2)\\ respectively, where m are the number of blocks to be used
(parameter nblocks). Theoretically safe values for \\m\\ are \\\<
n^{1/3}\\, but practically \\m\\ may be a little bit larger. The
function will issue a warning, if the value for \\m\\ is too large.

## References

Zhang, Y., Duchi, J.C., Wainwright, M.J., 2013. Divide and Conquer
Kernel Ridge Regression: A Distributed Algorithm with Minimax Optimal
Rates. arXiv:1305.5029 \[cs, math, stat\].

## See also

[`constructLearner`](https://rdrr.io/pkg/CVST/man/constructLearner.html)

## Examples

``` r
ns <- noisySinc(1000)
nsTest <- noisySinc(1000)

fast.krr <- constructFastKRRLearner()
fast.p <- list(kernel="rbfdot", sigma=100, lambda=.1/getN(ns), nblocks = 4)
system.time(fast.m <- fast.krr$learn(ns, fast.p))
#>    user  system elapsed 
#>   0.078   0.122   0.065 
fast.pred <- fast.krr$predict(fast.m, nsTest)
sum((fast.pred - nsTest$y)^2) / getN(nsTest)
#> [1] 0.01886754

if (FALSE) { # \dontrun{
krr <- CVST::constructKRRLearner()
p <- list(kernel="rbfdot", sigma=100, lambda=.1/getN(ns))
system.time(m <- krr$learn(ns, p))
pred <- krr$predict(m, nsTest)
sum((pred - nsTest$y)^2) / getN(nsTest)

plot(ns, col = '#00000030', pch = 19)
lines(sort(nsTest$x), fast.pred[order(nsTest$x)], col = '#00C000', lty = 2)
lines(sort(nsTest$x), pred[order(nsTest$x)], col = '#0000C0', lty = 2)
legend('topleft', legend = c('fast KRR', 'KRR'),
       col = c('#00C000', '#0000C0'), lty = 2)
} # }
```
