# \subsection{Original accept and reject algorithm}
# \label{sec:byrow-sample}
# 
# The following function runs the algorithm described in
# \cite{huber2003}. It picks the first row of the matrix and randomly
# selects a column, where the sample elements are weighted by the
# probability formula given above. If the last element of the sampling
# vector is selected (the one that has fills up the probability sum to
# one), the whole sample is rejected and restarted. This is accomplished
# by a \verb+repeat+ loop encapsulated by another \verb+repeat+ loop,
# and the Boolean variable \verb+accept+. The outermost loop is broken
# out of iff \verb+accept==TRUE+.
###  Algorithm from Figure 2    
ResperByrow <- function(mat) {
  reject <- 0
  repeat { # the outer loop; exited after proper sample
    # The accept flag is initialized to \verb+TRUE+, the row sums and 
    # dimension
    # are calculated, and the indices are initialized as a reference. These 
    # are
    # needed because the matrix is reduced to its minors later, to maintain 
    # the
    # reference from the column of a minor to the column of the matrix. Then
    # the inner loop is entered. This \verb+repeat+ loop can be broken out of
    # successfully (if a permutation is complete) or with a failure (if a
    # permutation failed to be sampled). If it is exited successfully, exit
    # this loop also, if not, increment the variable \verb+reject+, 
    # re-initialize and re-enter the inner loop.
    # The
    # variable \verb+reject+ records the number of failures. It is returned 
    # with the permutation in a
    # common data structure.
    accept <- TRUE
    at <- mat
    rt <- apply(at, 1, sum)
    n <- nrow(at)
    rownames(at) <- colnames(at) <- 1:n
    indizes <- 1:n
    perm <- NULL
    repeat { # the inner loop; exited on failed sample for retry
      # Failures within this loop can occur if the matrix includes
      # a zero row, or if the random sampling process selects the
      # ``column'' beyond the matrix.
      # 
      # If a column of the matrix is sampled, the vector of permutation is
      # extended by the selected column. Note that the column number has to
      # denote the column number of the original matrix, not of the current,
      # reduced, matrix. Therefore, the column numbers of the original matrix
      # are passed as colnames and referenced when the permanent is extended.
      ##  1:n doesn't change, even if n is changed
      ##  within the loop 
      for (i in 1:n) {
        ##  reject if any row with only zeros
        if (prod(rt)==0) {
          accept <- FALSE
          break
        ##  trap special case of 1 by 1 matrix  
        } else if (n==1) {
          elt <- 1
        } else {
          # First, create the vector of sample weights (permanent bounds of
          # all the minors), including
          # the extra element that completes the sum so that it
          # equals the permanent bound of the current matrix.
          # Then, sample and reject if this extra element is selected.
          Mrt <- HuberProbs(at, 1)
          elt <- sample(n+1, 1, prob=Mrt)
          if (elt==n+1) {
            accept <- FALSE
            break
          }
                            }
        ##  update permutation vector
        perm <- c(perm, indizes[elt])
        if (n > 1) {
          at <- at[-1, -elt]
          indizes <- indizes[-elt]
          if (is.matrix(at)) {
            ##  re-calculate row sums and dimension
            rt <- apply(at, 1, sum)
            n <- nrow(at)
            ##  trap special case of 1 by 1 matrix  
          } else {
            ##  don't reduce further, just
            ##  re-calculate row sums and dimension
            rt <- at
            n <- 1
          }   
        }
      }  
      break  #  proper permutation sampled
    }   
    if (accept) break else {
      reject <- reject+1
    }
  }
  list(perm=perm, reject=reject)
}
