getwd()
library(resper)

set.seed(seed=1, kind="default")

mklein <- toeplitz(c(1,1,0,0,0,0,0))
mgross <- toeplitz(c(rep(1,4), rep(0,36)))

mzufall <- WithinDeltaMat(cumsum(runif(20)), 2)

resper:::HuberGNonRecursive(5)
resper:::HuberGNonRecursive(280)
resper:::HuberProbs(mklein, 1)
resper:::HuberProbs(mklein, 4)
resper:::HuberProbs(mgross, 1)
resper:::HuberProbs(mgross, 14)
