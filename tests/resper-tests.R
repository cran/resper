getwd()
source('../../resper/R/helper.R')
source('../../resper/R/byrow-sample.R')
source('../../resper/R/clot-sample.R')

set.seed(seed=1, kind="default")

mklein <- toeplitz(c(1,1,0,0,0,0,0))
mgross <- toeplitz(c(rep(1,4), rep(0,36)))

mzufall <- WithinDeltaMat(cumsum(runif(20)), 2)

HuberGNonRecursive(5)
HuberGNonRecursive(280)
HuberProbs(mklein, 1)
HuberProbs(mklein, 4)
HuberProbs(mgross, 1)
HuberProbs(mgross, 14)
