%@+leo-ver=4
%@+node:@file ../../man/WithinDeltaMat.Rd
%@@language latex

\name{WithinDeltaMat}
\alias{WithinDeltaMat}
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
Generating an adjacency matrix for |x-y| < delta
%@nonl
%@-node:<< Title >>
%@afterref
}
\description{
  %@  <<Description>>
  %@+node:<< Description >>
  Given a sequence of numbers and a "window" delta,
  the function generates an adjacency matrix for the
  elements of the sequence and the relation |x-y| < delta,
  where x and y are two elements of the sequence.
  %@-node:<< Description >>
  %@nl
}
\usage{
WithinDeltaMat(seq, delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  %@  <<Arguments>>
  %@+node:<< Arguments >>
  \item{seq}{
  %@nonl
  %@<<seq>>
  %@+node:<< seq >>
  an ordered sequence of numbers
  %@nonl
  %@-node:<< seq >>
  %@afterref
}
  \item{delta}{
  %@nonl
  %@<<delta>>
  %@+node:<< delta >>
  a nonnegative number
  %@nonl
  %@-node:<< delta >>
  %@afterref
}
  %@nonl
  %@-node:<< Arguments >>
  %@nl
}
\details{
  %@  <<Details>>
  %@+node:<< Details >>
  %@-node:<< Details >>
  %@nl
  }
\value{
  %@  <<Value>>
  %@+node:<< Value >>
  A matrix of Boolean values
  (\code{TRUE} for a pair of elements whose distance
  is less than \code{delta}).
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
%@-node:<< Notes >>
%@afterref
 }

\seealso{ 
%@nonl
%@<<See Also>>
%@+node:<< See Also >>
%@-node:<< See Also >>
%@afterref
 }
\examples{
%@<<Examples>>
%@+node:<< Examples >>
##  this returns a Toeplitz matrix

WithinDeltaMat(1:10, 4)

##  this generates a matrix from a Poisson sequence

run <- cumsum(rgamma(100, 1))
WithinDeltaMat(run, 3)
%@nonl
%@-node:<< Examples >>
%@nl
}
\keyword{datagen}% at least one, from doc/KEYWORDS
%@-node:@file ../../man/WithinDeltaMat.Rd
%@-leo
