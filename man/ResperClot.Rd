%@+leo-ver=4
%@+node:@file ../../man/ResperClot.Rd
%@@language latex

\name{ResperClot}
\alias{ResperClot}
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
Sample from a set of time window restricted permutations
%@nonl
%@-node:<< Title >>
%@afterref
}
\description{
  %@  <<Description>>
  %@+node:<< Description >>
  Sample from a set of permutations of points in time, where two elements cannot be juxtaposed if the lag between them is longer than a predefined time window
  %@nonl
  %@-node:<< Description >>
  %@nl
}
\usage{
ResperClot(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  %@  <<Arguments>>
  %@+node:<< Arguments >>
  \item{mat}{
  %@nonl
  %@<<mat>>
  %@+node:<< mat >>
  A Boolean or 0-1 matrix with a ```tube'' structure 
  following the rule: 
  \deqn{i<j \Leftrightarrow A_{i,j+1}\le A_{i,j}\le A_{i+1,j}}{i<j ==> A_{i,j+1} <= A_{i,j} <= A_{i+1,j}}
  (for the Boolean matrix, replace \eqn{\le}{<=} with \eqn{\Rightarrow}{==>}). 
  The rows and columns of the matrix can be interpreted as
  time points, and  corresponds to
  the situations where the elements are time points, and 
  a permutation is not valid if an element \eqn{i}{i} occurs
  before another one \eqn{j}{j} when \eqn{t(i)>t(j)+\Delta}{t(i)>t(j)+delta}}. The matrix
  is given by the sequence and \eqn{\Delta}{delta}.
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
  %@+at
  % This function uses an accept-reject algorithm that tries
  % to break up the matrix into a block-diagonal structure.
  % It selects a row whose elimination will bring the matrix
  % closer to a block-diagonal structure, samples an appropriate
  % cell from the 1s in the row (denoting the position to
  % which the element corresponding to the row will be shifted), eliminates 
  % the cell row and column from the matrix and
  % reiterates the algorithm on the minor. The sampling
  % weights are constructed by an envelope probability
  % discovered by \cite{huber2003}, which guarantees that the
  % sum of weights over all minors is less than or equal
  % to the weight of the matrix itself. If the algorithm
  % samples the ``less than'' part, it will reject and
  % restart the current attempt.
  %@-at
  %@@c  
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
The matrix is not checked for a tube structure
when calling the function.
%@nonl
%@-node:<< Notes >>
%@afterref
 }

\seealso{ 
%@nonl
%@<<See Also>>
%@+node:<< See Also >>
\code{\link{WithinDeltaMat}}, \code{\link{ResperByrow}}  
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
perm <- ResperClot(mat)
perm

  
%@-node:<< Examples >>
%@nl
}
\keyword{datagen}
%@nonl
%@-node:@file ../../man/ResperClot.Rd
%@-leo
