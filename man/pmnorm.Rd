\name{pmnorm}
\alias{pmnorm}
\title{
MVN Rectangle Probabilities
}

\description{
Multivariate normal rectangle probabilities using Schervish's method
}

\usage{
pmnorm(lb, ub, mu, sigma, eps = 1.e-05)
}
\arguments{
\item{lb}{
   vector of lower limits of integral/probability
}
\item{ub}{
   vector of upper limits of integral/probability 
}
\item{mu}{
   mean vector of the multivariate normal density 
}
\item{sigma}{
  covariance matrix, it is assumed to be positive-definite
}
\item{eps}{
  tolerance for integration
}
}
\value{
\item{pr}{
  probability of the multivariate normal rectangle region
}
\item{perr}{
  estimated accuracy
}
\item{ifault}{
  return codes from the referenced paper

  = 0 if no problems
 
  = 1 or 2 if eps too small
 
  = 3 if dimension is not between 1 and 6 inclusive

  = 4 if covariance matrix is not positive-definite
}
}

\references{
  Schervish, M.J. (1984). Multivariate normal probabilities with error bound.
  Appl. Statist., 33, 81-94.
}
\seealso{
  \code{\link{mvnapp}}.
}
\author{H. Joe, Statistics Department, UBC}
\examples{

rh<-0.3
m<-2
a<-c(-1,-1)
b<-c(1,1)
mu<-rep(0,m)
s<-matrix(c(1,rh,rh,1),2,2)
print(pmnorm(a,b,mu,s))

m<-3
a<-c(-1,-1,-2)
b<-c(1,1,.5)
mu<-rep(0,m)
s<-matrix(c(1,rh,rh,rh,1,rh,rh,rh,1),3,3)
print(pmnorm(a,b,mu,s))

m<-4
a<-c(-1,-2.5,-2,-1.5)
b<-c(1.68,1.11,.5,.25)
mu<-rep(0,m)
s<-matrix(c(1,rh,rh,rh,rh,1,rh,rh,rh,rh,1,rh,rh,rh,rh,1),4,4)
print(pmnorm(a,b,mu,s))
}
\keyword{models}
