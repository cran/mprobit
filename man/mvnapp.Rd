\name{mvnapp}
\alias{mvnapp}
\title{
MVN Rectangle Probabilities
}
\description{
Approximation to multivariate normal rectangle probabilities using methods in Joe
(1995), JASA
}
\usage{
mvnapp(lb,ub,mu,sigma,type=1,eps=1.e-05,nsim=0)
}
\arguments{
\item{lb}{
  vector of lower limits of integral/probability
}
\item{ub}{
  vector of upper limits of integral/probability
}
\item{mu}{
  mean vector
}
\item{sigma}{
  covariance matrix, it is assumed to be positive-definite
}
\item{type}{
  indicator, type=1 refers to the first order approximation, 
  type=2 is the second order approximation.
}
\item{eps}{
  accuracy/tolerance for bivariate marginal rectangle probabilities
} 
\item{nsim}{
  an optional integer if random permutations are used in the approximation
  for dimension >=6; nsim=2000 recommended for dim>=6
} 
}
\value{
\item{prob}{
  rectangle probability with approximation
}
\item{esterr}{
  indicator of accuracy in the approximation
}
\item{ifail}{
  = 0 if no problems
 
  >= 1 if problems from using Schervish's code in dimensions 2 to 4.
}
}

\references{
  Joe, H (1995). Approximations to multivariate  normal rectangle
  probabilities based on conditional expectations.  J. Amer. Statist.
  Assoc., 90, 957-964.

  Schervish, M.J. (1984). Multivariate normal probabilities with error bound.
  Appl. Statist., 33, 81-94.
}
\seealso{
  \code{\link{pmnorm}}.
}
\author{H. Joe, Statistics Department, UBC}
\examples{
m<-2
rh<-0.5
a<-c(-1,-1)
b<-c(1,1)
mu<-rep(0,m)
s<-matrix(c(1,rh,rh,1),2,2)
print(pmnorm(a,b,mu,s))
print(mvnapp(a,b,mu,s))
print(mvnapp(a,b,mu,s,type=2))
print(mvnapp(a,b,mu,s,type=2,nsim=3))

m<-3
rh<-0.3
a<-c(-1,-1,-2)
b<-c(1,1,.5)
mu<-rep(0,m)
s<-matrix(c(1,.5,.6,.5,1,.7,.6,.7,1),3,3)
print(pmnorm(a,b,mu,s))
print(mvnapp(a,b,mu,s))
print(mvnapp(a,b,mu,s,type=2))
print(mvnapp(a,b,mu,s,type=2,nsim=3))

m<-4
rh<- -0.1
a<-c(-1,-2.5,-2,-1.5)
b<-c(1.68,1.11,.5,.25)
mu<-rep(0,m)
s<-matrix(c(1,.5,.3,.4,.5,1,.5,.4,.3,.5,1,.4,.4,.4,.4,1),m,m)
print(pmnorm(a,b,mu,s))
print(mvnapp(a,b,mu,s))
print(mvnapp(a,b,mu,s,type=2))
print(mvnapp(a,b,mu,s,type=2,nsim=3))

m<-5
rh<-.4
a<-rep(-1,m)
b<-rep(2,m)
mu<-rep(0,m)
s<-matrix(c(1,rh,rh,rh,rh,rh,1,rh,rh,rh,rh,rh,1,rh,rh,rh,rh,rh,1,
       rh,rh,rh,rh,rh,1),m,m)
print(mvnapp(a,b,mu,s))
print(mvnapp(a,b,mu,s,type=2))
print(mvnapp(a,b,mu,s,type=2,nsim=3))

m<-6
a<-c(-1,-1,-1,-1.5,-1,-2)
b<-rep(7,m)
mu<-rep(0,m)
s<-matrix(c(1,rh,rh,rh,rh,rh,rh,1,rh,rh,rh,rh,rh,rh,1,rh,rh,rh,rh,rh,rh,1,
       rh,rh,rh,rh,rh,rh,1,rh,rh,rh,rh,rh,rh,1),m,m)
print(mvnapp(a,b,mu,s))
print(mvnapp(a,b,mu,s,type=2))
print(mvnapp(a,b,mu,s,type=2,nsim=3))
}
\keyword{models}
