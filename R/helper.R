# \subsection{Auxiliary functions}
# \label{sec:helper}
# 
# \subsection{The recursive function}
# \label{sec:recfun}
# The novelty introduced by Huber is the variation of a
# bound of the permanent from above \cite{bregman1973} that works as a 
# probability: The
# value of a function of a matrix is always at least as big as the sum
# of the values of the function of the minors, developed around a column
# of the matrix. It is the recursive function on the dimension $n$ of
# the matrix:
# \begin{eqnarray}
#   \label{eq:huberg}
#   G(n) :=
#   \begin{cases}
#      \text{e} & \text{for }n=1 \\
#       G(n-1)+1+0.5/G(n-1)+0.6/G(n-1)^2 & \text{for }n>1
#   \end{cases}
# \end{eqnarray}

HuberGInner <- function(n) {
    if (n == 1) 
        exp(1)
    else {
        gnm <- HubergGInner(n - 1)
        gnm + 1 + 0.5/gnm + 0.6/gnm/gnm
    }
}
 
# R does not like recursion too much, so a non-recursive version is
# employed here:
HuberGNonRecursive <- function(n) {
    for (i in 1:n) {
        if (i == 1) 
            gnm <- exp(1)
        else gnm <- gnm + 1 + 0.5/gnm + 0.6/gnm/gnm
    }
    gnm
}
# 
# A wrapper for the function is used that checks for valid entries:
# 
HuberGE <- function(n) {
    if (n > floor(n)) {
        stop("n must be of integer value")
    }
    else if (n < 0) {
        stop("n must be at least 1")
    }
    else if (n == 0) 
        0
    else {
        HuberGNonRecursive(n)/exp(1)
    }
}
# \subsection{The selection probability}
# \label{sec:selectprob}
# 
# The actual probability is given by
# \begin{eqnarray}
#   \label{eq:permboundary}
#   M(A) & := & \prod_i G(c_i) ,
# \end{eqnarray}
# where $c_i$ is
# the sum of the $i$th column (i.e. the number of ones in it).
# 
###  Formula (3)
PermBound <- function(mat) {
        prod(sapply(apply(mat, 2, sum), HuberGE))
}
# 
# As said, this is a probability because
# \begin{eqnarray}
#   \label{eq:probability}
#   M(A) \ge \sum_j a_{ij} M(A(\breve{\imath}, \breve{\jmath})) 
# \quad\forall\quad i.
# \end{eqnarray}
#  The \verb+drop=FALSE+ argument retains the array structure even if the 
# dimension of the matrix is 1 (see \cite{R-FAQ}). The vector of the above sum 
# elements is needed to sample a row for a given column. The function pastes 
# the difference of the sum of these elements from one to the end.
# 
HuberProbs <- function(at, i) {
    n <- nrow(at)
    Mrt <- sapply(1:n, function(j) {
        if (at[i, j] > 0) {
            PermBound(at[-i, -j, drop=FALSE])
        }
        else {
            0
        }
    })
    c(Mrt, PermBound(at) - sum(Mrt))
}
# The most important application of a tube matrix
# in this context is the one that determines
# permutability between elements that are close
# to each other. Therefore, a function is desirable that
# returns such a matrix from an ordered sequence
# $(t_i)_{i\in \{1,\ldots , n\}}$ and a lag $\Delta$.
# The rows and columns of the matrix
# correspond to the elements in the ordered seqence,
# and an element $a_{ij}$ is set to \verb+TRUE+
# exactly when $|t_i-t_j| < \Delta$. The Boolean entries
# in the matrix are converted to 0s and 1s when arithmetic
# functions are applied to them.
WithinDeltaMat <- function(seq, delta) {
  outer(seq, seq+delta, '<') & outer(seq+delta, seq, '>')
}

