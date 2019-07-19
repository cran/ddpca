\name{HCdetection}
\alias{HCdetection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Higher Criticism for detecting rare and weak signals
}
\description{
This function takes a bunch of p-values as input and ouput the Higher Criticism statistics as well as the decision (rejection or not).
}
\usage{
HCdetection(p, alpha = 0.5, pvalcut = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
A vector of size \code{n} containing p-values from data
}
  \item{alpha}{
A number between 0 and 1. The smallest alpha*n p-values will be
used to calculate the HC statistic. Default is 0.5.
}
  \item{pvalcut}{
A number between 0 and 1. Those small p-values (smaller than
pvalcut) will be taken away to avoid heavy tails of test
statistic. Set it to \code{NA} is equivalent to setting it to \eqn{1/n}.
}
}
\details{
This function is an adaptation of the Matlab code here \url{http://www.stat.cmu.edu/~jiashun/Research/software/HC/}
}
\value{
Returns a list containing the following items
\item{H}{0 or 1 scalar indicating whether \eqn{H_0} the global null is rejected (1) or not rejected (0)}
\item{HCT}{Higher Criticism test statistic}
%% ...
}
\references{
Donoho, D. and Jin, J., Higher criticism for detecting sparse heterogeneous mixtures. Ann. Statist. 32 (2004), no. 3, 962--994. 

Ke, Z., Xue, L. and Yang, F., 2019. Diagonally Dominant Principal Component Analysis. Journal of Computational and Graphic Statistics, under review.
}
\author{
Fan Yang <fyang1@uchicago.edu>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
n = 1e5
data = rnorm(n)
p = 2*(1 - pnorm(abs(data)))
result = HCdetection(p)
print(result$H)
print(result$HCT)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Higher Criticism}

