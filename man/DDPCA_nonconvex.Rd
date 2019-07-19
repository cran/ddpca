\name{DDPCA_nonconvex}
\alias{DDPCA_nonconvex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagonally Dominant Principal Component Analysis using Nonconvex approach
}
\description{
This function decomposes a positive semidefinite matrix into a low rank component, and a diagonally dominant component using an iterative projection algorithm. 
}
\usage{
DDPCA_nonconvex(Sigma, K, max_iter_nonconvex = 15, SDD_approx = TRUE, 
max_iter_SDD = 20, eps = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sigma}{
Input matrix of size \eqn{n\times n}
}
  \item{K}{
A positive integer indicating the rank of the low rank component.
}
  \item{max_iter_nonconvex}{
Maximal number of iterations of the iterative projection algorithm. 
}
  \item{SDD_approx}{
If set to TRUE, then the projection onto SDD cone step in each iteration will be replaced by projection onto DD cone followed by symmetrization. This approximation will reduce the computational cost, but the output matrix \code{A} may only be approximately diagonally dominant. 
}
  \item{max_iter_SDD,eps}{
Arguments in function \code{ProjSDD}. Matters only when \code{SDD_approx = False}.
}

}
\details{
This function performs iterative projection algorithm to decompose a positive semidefinite matrix \code{Sigma} into a low rank component \code{L} and a symmetric diagonally dominant component \code{A}. The projection onto the set of low rank matrices is done via eigenvalue decomposition, while the projection onto the symmetric diagonally dominant (SDD) cone is done via function \code{ProjSDD}, unless \code{SDD_approx = TRUE} where an approximation is used to speed up the algorithm. 
}
\value{
A list containing the following items
\item{L}{The low rank component}
\item{A}{The diagonally dominant component}
%% ...
}
\references{
Ke, Z., Xue, L. and Yang, F., 2019. Diagonally Dominant Principal Component Analysis. Journal of Computational and Graphic Statistics, under review.
}
\author{
Fan Yang <fyang1@uchicago.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DDPCA_convex}}
}
\examples{
library(MASS)
p = 200
n = 200
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
result = DDPCA_nonconvex(Sigma_sample,K=k)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ddpca}
\keyword{nonconvex}