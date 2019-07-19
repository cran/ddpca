\name{DDPCA_convex}
\alias{DDPCA_convex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagonally Dominant Principal Component Analysis using Convex approach
}
\description{
This function decomposes a positive semidefinite matrix into a low rank component, and a diagonally dominant component by solving a convex relaxation using ADMM. 
}
\usage{
DDPCA_convex(Sigma, lambda, rho = 20, max_iter_convex = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sigma}{
Input matrix of size \eqn{n\times n}
}
  \item{lambda}{
The parameter in the convex program that controls the rank of the low rank component 
}
  \item{rho}{
The parameter used in each ADMM update.
}
  \item{max_iter_convex}{
Maximal number of iterations of ADMM update. 
}
}
\details{
This function decomposes a positive semidefinite matrix \code{Sigma} into a low rank component \code{L} and a symmetric diagonally dominant component \code{A}, by solving the following convex program 
\deqn{\textrm{minimize} \quad 0.5*\|\Sigma - L - A\|^2 + \lambda \|L\|_{*}}
\deqn{\textrm{subject to} \quad A\in SDD}
where \eqn{\|L\|_{*}}{||L||_{*}} is the nuclear norm of \code{L} (sum of singular values) and \code{SDD} is the symmetric diagonally dominant cone. 

}
\value{
A list containing the following items
\item{L}{The low rank component}
\item{A}{The diagonally dominant component}
}
\references{
Ke, Z., Xue, L. and Yang, F., 2019. Diagonally Dominant Principal Component Analysis. Journal of Computational and Graphic Statistics, under review.
}
\author{
Fan Yang <fyang1@uchicago.edu>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DDPCA_nonconvex}}
}
\examples{
library(MASS)
p = 30
n = 30
k = 3
rho = 0.5
a = 0:(p-1)
Sigma_mu = rho^abs(outer(a,a,'-'))
Sigma_mu = (diag(p) + Sigma_mu)/2 # Now Sigma_mu is a symmetric diagonally dominant matrix
mu = mvrnorm(n,rep(0,p),Sigma_mu)
B = matrix(rnorm(p*k),nrow = p)
F = matrix(rnorm(k*n),nrow = k)
Y = mu + t(B \%*\% F) 
Sigma_sample = cov(Y)
result = DDPCA_convex(Sigma_sample,lambda=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ddpca}
\keyword{convex}