# \subsection{The separating accept and reject algorithm}
# \label{sec:clot-sample}
# This function uses an accept-reject algorithm that tries
# to break up the matrix into a block-diagonal structure.
# It selects a row whose elimination will bring the matrix
# closer to a block-diagonal structure, samples an appropriate
# cell from the 1s in the row (denoting the position to
# which the element corresponding to the row will be shifted), eliminates the 
# cell row and column from the matrix and
# reiterates the algorithm on the minor. The sampling
# weights are constructed by an envelope probability
# discovered by \cite{huber2003}, which guarantees that the
# sum of weights over all minors is less than or equal
# to the weight of the matrix itself. If the algorithm
# samples the ``less than'' part, it will reject and
# restart the current attempt.
# 
# The modified algorithm does not take the rows in the given order,
# but picks a row by certain criteria.
# 
# The first criterion addresses possible separators,
# that is, sets of few columns that, when removed,
# leave the matrix with a block diagonal structure.
# If one views the matrix as an adjacency matrix of
# a graph, the task is now to look for waists of the
# graph.
# 
# As one can rely on the matrix having a tube structure,
# a primitive algorithm is sufficient to look for
# waists in the corresponding graph.
# These columns are detected by counting the number
# of ones below the main diagonal.
# 
###  selection criteria for separator
firstcrit <- function(mat) {
  n <- nrow(mat)
  lowerdiag <- outer(1:n, 1:n, ">=")
  ##  the most desirable is the column _after_ the one
  ##  with the least ones from the main diagonal
  apply(cbind(rep(1,n), mat[,1:n-1])*
        ##  do not choose columns with 1s to the bottom
        ##  therefore heavy weight to bottom row
        array(rep(c(rep(1, n-1), n), n), dim=c(n,n)) *
        lowerdiag, 2, sum) 
}
# 
# If one aims to split a matrix into separate blocks,
# one would like the blocks to have the same size.
# Therefore, the second criterion favors
# the middle rows and therefore is a convex,
# symmetric function of the row index.
# 
secondcrit <- function(mat) {
    n <- nrow(mat)
    (1:n) * (n:1)/(n + 1)/(n + 1) * 8
}
# 
# The row with the minimal sum of first and second criterion is selected.
# 
selectrow <- function(mat) {
  which.min(firstcrit(mat)-secondcrit(mat))
}
# The main idea of speeding up the process is to
# select the rows to reduce the matrix such as to
# obtain minors with a block-diagonal structure.
# If a block-diagonal structure is detected, the
# algorithm calls itself recursively on the blocks
# and pastes the results together to get the whole
# permutation.
# 
# Again, this is not a function that works on
# general symmetric 0-1 matrices.
# Instead of globally searching for a block structure,
# it looks if subsequent rows have at least one zero
# in either column.
seps <- function(mat) {
    if (dim(mat)[1] < 3) 
        NULL
    else {
        n <- nrow(mat)
        sumoff <- c(sum(mat[2:n, 1] + mat[1, 2:n]), sapply(2:(n - 
            1), function(i) {
            sum(mat[1:i, (i + 1):n]) + sum(mat[(i + 1):n, 1:i])
        }))
        which(sumoff == 0)
    }
}
# 
# 
# The next row is selected according to the two criteria, and the column
# is sampled according to the weights obtained by equation \ref{eq:huber-g}. 
# If the element beyond the matrix rows is
# selected, the current sample is rejected and restarted.
ResperClotInner <- function(mat) {
  reject <- 0
  repeat {
    # The permutation structure is initialized as a two-column
    # matrix, the first column denoting the row indices and the
    # second the column indices. The outer wrapper function
    # converts these to a permutation vector.
    # 
    # The acceptance flag is initialized to \verb+TRUE+, the
    # row sums and dimension are calculated. The rownames and
    # the colnames are initialized to the indices if they are
    # not present. If there are rownames and colnames already,
    # do not overwrite them as the function can be called
    # recursively.
    # 
    # Then, sample, reject, and try again until a proper sample
    # is selected.
    perm <- c(NULL, NULL)
    accept <- TRUE
    n <- nrow(mat)
    if (is.null(colnames(mat))) {
      colnames(mat) <- 1:n
      rownames(mat) <- 1:n
    }
    at <- mat
    rt <- apply(at, 1, sum)
    repeat {
      # First, trap the special case where the matrix has only
      # one dimension. In this case, fill the permanent structure
      # with the last row and column index, and return successfully.
      if (n==1) {
        perm <- rbind(perm, as.numeric(c(rownames(at), colnames(at))))
        break
      }
      # If the matrix contains only ones, one can sample from
      # the unrestricted set of permutations.
      else if (prod(at)==1) {
        # In the unrestricted case, one can simply use R's \verb+sample()+
        # algorithm. which is used for the column index column of the
        # permutation structure, which is then updated row-wise by the
        # permutation of the current matrix.
        unrestricted <- array(as.numeric(c(rownames(at),
                                           sample(colnames(at),
                                                  size=ncol(at)))),
                              dim=c(nrow(at),2))
        perm <- rbind(perm, unrestricted) 
        ## exit successfully 
        break
      } 
      # If there is a block-diagonal structure, call
      # function recursively on the blocks.
      else if (length(seps <- seps(at))>0) {
        bstart <- c(1, seps+1)
        bend <- c(seps, n)
        for (i in (1:length(bstart))) {
          ## recursively call the function on the blocks
          a <- ResperClotInner(at[bstart[i]:bend[i],
                                    bstart[i]:bend[i], drop=FALSE])
          perm <- rbind(perm, a$perm)                          
          reject <- reject+a$reject
        }
        ##  exit successfully
        break 
      }  
      # The next test is on the main diagonal containing a 0. By virtue of the
      # tube structure of the original matrix, if any of its minors has a 0 in
      # the main diagonal, this minor necessarily contains a rectangular
      # submatrix of only 0s that includes either both the first row and the
      # last column or the last row and the first column. This is a sufficient 
      # condition for the permanent of this minor being 0 (a result cited in 
      # \cite{minc1978}), which in turn is reason enough to reject and 
      # restart.
      else if (prod(diag(at))==0) {
        accept <- FALSE
        break
      }
      #   Now that we've handled the special cases,
      #   let's treat the normal case.
      else {
        # The row is sampled according to the optimality
        # criteria. The column is sampled at random.
        # If the element beyond the matrix columns is
        # selected, the accept-reject algorithm rejects.
        i <- selectrow(at) 
        Mrt <- HuberProbs(at, i)
        j <- sample(n+1, 1, prob=Mrt)
        if (j==n+1) {
          accept <- FALSE
          break
        }        ##  update permutation structure
        perm <- rbind(perm, as.numeric(c(rownames(at)[i],
                                         colnames(at)[j])))
        at <- at[-i,-j, drop=FALSE]
        n <- n-1
        if (dim(at)[1]==0) break
      }    }
    if (accept) {
      break
      } 
    else {
      reject <- reject+1
    }
    
          }
  list(perm=perm, reject=reject)
}
# 
# Finally, a wrapper function is written that reduces the permanent data
# structure to a single permanent vector, as in the function for the
# Huber algorithm.
# 
ResperClot <- function(mat) {
  a <- ResperClotInner(mat)
  a$perm <- a$perm[order(a$perm[,1]), 2]
  a
}
