%@+leo-ver=4
%@+node:@file ../../man/ResperByrow.Rd
%@@language latex

\name{ResperByrow}
\alias{ResperByrow}
%- Also NEED an '\alias' for EACH other topic documented here.
\encoding{
%@nonl
%@<<Encoding>>
%@+node:<< Encoding >>
utf-8
%@nonl
%@-node:<< Encoding >>
%@afterref
}
\title{
%@nonl
%@<<Title>>
%@+node:<< Title >>
Sample from a set of restricted permutations
%@nonl
%@-node:<< Title >>
%@afterref
}
\description{
  %@  <<Description>>
  %@+node:<< Description >>
  Sample from a set of permutations restricted by a Boolean 
  or 0-1 matrix \eqn{A}{A}, where the rows correspond to the original
  elements and the columns to the positions in the permutation.
  \eqn{A_{i,j}=1}{A_{i,j}=1} means that the \eqn{i}{i}th element may occur at the
  \eqn{j}{j}th position in the permuted sequence.    
  %@nonl
  %@-node:<< Description >>
  %@nl
}
\usage{
ResperByrow(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  %@  <<Arguments>>
  %@+node:<< Arguments >>
  \item{mat}{
  %@nonl
  %@<<mat>>
  %@+node:<< mat >>
  A Boolean or 0-1 matrix. A permutation can only be
  gemerated if at least one element in each row and
  each column is equal to 1. 
  %@nonl
  %@-node:<< mat >>
  %@afterref
}
  %@nonl
  %@-node:<< Arguments >>
  %@nl
}
\details{
  %@  <<Details>>
  %@+node:<< Details >>
  This function uses an accept-reject algorithm.
  It selects the first row, samples an appropriate 
  cell from the 1s in the row (denoting the position to 
  which the element corresponding to the row will be shifted), eliminates the cell row and column from the matrix and 
  reiterates the algorithm on the minor. The sampling 
  weights are constructed by an envelope probability 
  discovered by \cite{huber2003}, which guarantees that the 
  sum of weights over all minors is less than or equal
  to the weight of the matrix itself. If the algorithm
  samples the ``less than'' part, it will reject and
  restart the current attempt.  
  %@nonl
  %@-node:<< Details >>
  %@nl
  }
\value{
  %@  <<Value>>
  %@+node:<< Value >>
  A sequence of values corresponding to the row indices 
  of \verb+t+.
  %@nonl
  %@-node:<< Value >>
  %@nl
}
\references{ }
\author{ 
%@nonl
%@<<Author>>
%@+node:<< Author >>
Johannes Hüsing
%@-node:<< Author >>
%@afterref
 }
\note{ 
%@nonl
%@<<Notes>>
%@+node:<< Notes >>
The function does not check if the input is a proper 0-1 matrix.
%@nonl
%@-node:<< Notes >>
%@afterref
 }

\seealso{ 
%@nonl
%@<<See Also>>
%@+node:<< See Also >>
\code{link{WithinDeltaMat}}, \code{link{ResperClot}}  
%@nonl
%@-node:<< See Also >>
%@afterref
 }
\examples{
%@<<Examples>>
%@+node:<< Examples >>
##  No element may be shifted by more than two 
##  positions

mat <- toeplitz(c(rep(1,3), rep(0,3)))
perm <- ResperByrow(mat)
perm

  
%@-node:<< Examples >>
%@nl
}
\keyword{datagen}
%@nonl
%@-node:@file ../../man/ResperByrow.Rd
%@-leo
