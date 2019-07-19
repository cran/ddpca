\name{ProjSDD}
\alias{ProjSDD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Projection onto the Symmetric Diagonally Dominant Cone
}
\description{
Given a matrix C, this function outputs the projection of C onto the cones of symmetric diagonally domimant matrices using Dykstra's projection algorithm. 
}
\usage{
ProjSDD(A, max_iter_SDD = 20, eps = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Input matrix of size \eqn{n\times n}
}
  \item{max_iter_SDD}{
Maximal number of iterations of the Dykstra's projection algorithm 
}
  \item{eps}{
The iterations will stop either when the Frobenious norm of difference matrix between two updates is less than \code{eps} or after \code{max_iter_SDD} steps. If set to \code{NA}, then no check will be done during iterations and the iteration will stop after \code{max_iter_SDD} steps. Default is \code{NA}.
}
}
\details{
This function projects the input matrix \eqn{C} of size \eqn{n\times n} onto the cones of symmetric diagonally domimant matrices defined as 
\deqn{ \{A = (a_{ij})_{1\le i\le n, 1\le j\le n} : a_{ij} = a_{ji}, a_{jj} \ge \sum_{k\not=j} |a_{jk}| \quad \textrm{for all} \quad 1\le j\le n, 1\le i\le n \}}

It makes use of Dykstra's algorithm, which is a variation of iterative projection algorithm. The two key steps are projection onto the diagonally domimant cone by calling function \code{ProjDD} and projection onto the symmetric matrix cone by simple symmetrization. 

More details can be found in Mendoza, M., Raydan, M. and Tarazaga, P., 1998. Computing the nearest diagonally dominant matrix.
}
\value{
A \eqn{n\times n} symmetric diagonally dominant matrix
}
\references{
Mendoza, M., Raydan, M. and Tarazaga, P., 1998. Computing the nearest diagonally dominant matrix. Numerical linear algebra with applications, 5(6), pp.461-474.

Ke, Z., Xue, L. and Yang, F., 2019. Diagonally Dominant Principal Component Analysis. Journal of Computational and Graphic Statistics, under review.
}
\author{
Fan Yang <fyang1@uchicago.edu>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ProjDD}}
}
\examples{
ProjSDD(matrix(runif(100),nrow=10))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{projection}
