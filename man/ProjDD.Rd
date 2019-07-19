\name{ProjDD}
\alias{ProjDD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Projection onto the Diagonally Dominant Cone
}
\description{
Given a matrix C, this function outputs the projection of C onto the cones of diagonally domimant matrices.
}
\usage{
ProjDD(C)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{C}{
A \eqn{n \times n} matrix
}
}
\details{
This function projects the input matrix \eqn{C} of size \eqn{n\times n} onto the cones of diagonally domimant matrices defined as 
\deqn{ \{A = (a_{ij})_{1\le i\le n, 1\le j\le n} : a_{jj} \ge \sum_{k\not=j} |a_{jk}| \quad  \textrm{for all} \quad 1\le j\le n \}}
The algorithm is described in Mendoza, M., Raydan, M. and Tarazaga, P., 1998. Computing the nearest diagonally dominant matrix.
}
\value{
A \eqn{n\times n} diagonally dominant matrix
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
\code{\link{ProjSDD}}
}
\examples{
ProjDD(matrix(runif(100),nrow=10))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{projection}