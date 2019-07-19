\name{DDHC}
\alias{DDHC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DD-HC test
}
\description{
Combining DDPCA with orthodox Higher Criticism for detecting sparse mean effect.
}
\usage{
DDHC(X, known_Sigma = NA, method = "nonconvex", K = 1, lambda = 3, 
max_iter_nonconvex = 15 ,SDD_approx = TRUE, max_iter_SDD = 20, eps = NA, 
rho = 20, max_iter_convex = 50, alpha = 0.5, pvalcut = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A \eqn{n\times p} data matrix, where each row is drawn i.i.d from \eqn{\mathcal{N}(\mu,\Sigma)} 
}
  \item{known_Sigma}{
The true covariance matrix of data. Default NA. If NA, then \eqn{\Sigma} will be estimated from data matrix X. 
}
  \item{method}{
Either "convex" or "noncovex", indicating which method to use for DDPCA.
}
  \item{K}{
Argument in function \code{DDPCA_nonconvex}. Need to be specified when \code{method = "nonconvex"}
}
  \item{lambda}{
Argument in function \code{DDPCA_convex}. Need to be specified when \code{method = "convex"}
}
  \item{max_iter_nonconvex}{
Argument in function \code{DDPCA_nonconvex}. 
}
 \item{SDD_approx}{
Argument in function \code{DDPCA_nonconvex}. 
 }
 \item{max_iter_SDD}{
Argument in function \code{DDPCA_nonconvex}. 
 }
 \item{eps}{
 Argument in function \code{DDPCA_nonconvex}. 
 }
 \item{rho}{
 Argument in function \code{DDPCA_convex}. 
 }
 \item{max_iter_convex}{
 Argument in function \code{DDPCA_convex}. 
 }
 \item{alpha}{
 Argument in function \code{HCdetection}.
 }
 \item{pvalcut}{
 Argument in function \code{HCdetection}.
 }
}
\details{
See reference paper for more details. 
}
\value{
Returns a list containing the following items
\item{H}{0 or 1 scalar indicating whether \eqn{H_0} the global null is rejected (1) or not rejected (0). The use of \code{H} is not recommended as it's approximately valid only when \code{p} is sufficiently large and mean effect in alternative is really sparse.}
\item{HCT}{DD-HC Test statistic}
}
\references{
Ke, Z., Xue, L. and Yang, F., 2019. Diagonally Dominant Principal Component Analysis. Journal of Computational and Graphic Statistics, under review.
}
\author{
Fan Yang <fyang1@uchicago.edu>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{IHCDD}}, \code{\link{HCdetection}}, \code{\link{DDPCA_convex}}, \code{\link{DDPCA_nonconvex}}
}
\examples{
library(MASS)
n = 200
p = 200
k = 3
rho = 0.5
a = 0:(p-1)
Sigma_mu = rho^abs(outer(a,a,'-'))
Sigma_mu = (diag(p) + Sigma_mu)/2 # Now Sigma_mu is a symmetric diagonally dominant matrix
B = matrix(rnorm(p*k),nrow = p)
Sigma = Sigma_mu + B \%*\% t(B)
X = mvrnorm(n,rep(0,p),Sigma)
results = DDHC(X,K = k)
print(results$H)
print(results$HCT)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ddpca }
\keyword{ Higher Criticism}
