\name{mvn.deriv}
\alias{mvn.deriv}
\alias{mvn.deriv.margin}
\alias{mvn.deriv.rho}
\title{
Derivatives of Multivariate Normal Rectangle Probabilities
}
\description{
Derivatives of Multivariate Normal Rectangle Probabilities based on
Approximations
}
\usage{
mvn.deriv.margin(lb,ub,mu,sigma,k,ksign,type=1,eps=1.e-05,nsim=0)
mvn.deriv.rho(lb,ub,mu,sigma,j1,k1,type=1,eps=1.e-05,nsim=0)
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
\item{k}{
  margin for which derivative is to be taken, that is,
  deriv of mvnapp(lb,ub,mu,sigma) with respect to lb[k] or ub[k];
}
\item{ksign}{
  =-1 for deriv of mvnapp(lb,ub,mu,sigma) with respect to lb[k]
  =+1 for deriv of mvnapp(lb,ub,mu,sigma) with respect to ub[k]
}
\item{j1}{
  correlation for which derivative is to be taken, that is,
  deriv of mvnapp(lb,ub,mu,sigma) with respect to rho[j1,k1],
  where rho is a correlation corresponding to sigma
}
\item{k1}{
  See above explanation with j1
}
}
\value{
  derivative with respect to margin lb[k], ub[k], or correlation rho[j1][k1]
  corresponding to sigmamatrix
}

\references{
  Joe, H. (1995) JASA  TODO, complete
}
\seealso{
  \code{\link{mvnapp}}
}
\author{H. Joe, Statistics Department, UBC}
\examples{
# step size for numerical derivatives (accuracy of mvnapp etc may be about 1.e-4 to 1.e-5)
heps = 1.e-3

cat("compare numerical and analytical derivatives based on mvnapp\n")

cat("\ncase 1: dim=3\n");
m=3
mu=rep(0,m)
a=c(0,0,0)
b=c(1,1.5,2)
rr=matrix(c(1,.3,.3,.3,1,.4,.3,.4,1),m,m)
pr=mvnapp(a,b,mu,rr)$pr
# not checking ifail returned from mvnapp
cat("pr=mvnapp(avec,bvec,mu=0,sigma=corrmat)=",pr,"\n")

cat("derivative wrt a_k,b_k, k=1,...,",m,"\n")
for(k in 1:m)
{ cat("  k=", k, " lower\n")
  a2=a
  a2[k]=a[k]+heps
  pr2=mvnapp(a2,b,mu,rr)$pr
  da.numerical = (pr2-pr)/heps 
  da.analytic= mvn.deriv.margin(a,b,mu,rr,k,-1)$deriv
  cat("   numerical: ", da.numerical, ", analytic: ", da.analytic,"\n")
  cat("  k=", k, " upper\n")
  b2=b
  b2[k]=b[k]+heps
  pr2=mvnapp(a,b2,mu,rr)$pr
  db.numerical = (pr2-pr)/heps 
  db.analytic= mvn.deriv.margin(a,b,mu,rr,k,1)$deriv
  cat("   numerical: ", db.numerical, ", analytic: ", db.analytic,"\n")
}

cat("derivative wrt rho(j,k)\n")
for(j in 1:(m-1))
{ for(k in (j+1):m)
  { cat("  (j,k)=", j,k,"\n")
    rr2=rr
    rr2[j,k]=rr[j,k]+heps
    rr2[k,j]=rr[k,j]+heps
    pr2=mvnapp(a,b,mu,rr2)$pr
    drh.numerical = (pr2-pr)/heps 
    drh.analytic= mvn.deriv.rho(a,b,mu,rr2,j,k)$deriv
    cat("   numerical: ", drh.numerical, ", analytic: ", drh.analytic,"\n")
  }
}

#=============================================

cat("\ncase 2: dim=5\n");
m=5
mu=rep(0,m)
a=c(0,0,0,-1,-1)
b=c(1,1.5,2,2,2)
rr=matrix(c(1,.3,.3,.3,.4, .3,1,.4,.4,.4,  .3,.4,1,.4,.4,
  .3,.4,.4,1,.4, .4,.4,.4,.4,1),m,m)
pr=mvnapp(a,b,mu,rr)$pr
# not checking ifail returned from mvnapp
cat("pr=mvnapp(avec,bvec,mu=0,sigma=corrmat)=",pr,"\n")

cat("derivative wrt a_k,b_k, k=1,...,",m,"\n")
for(k in 1:m)
{ cat("  k=", k, " lower\n")
  a2=a
  a2[k]=a[k]+heps
  pr2=mvnapp(a2,b,mu,rr)$pr
  da.numerical = (pr2-pr)/heps 
  da.analytic= mvn.deriv.margin(a,b,mu,rr,k,-1)$deriv
  cat("   numerical: ", da.numerical, ", analytic: ", da.analytic,"\n")
  cat("  k=", k, " upper\n")
  b2=b
  b2[k]=b[k]+heps
  pr2=mvnapp(a,b2,mu,rr)$pr
  db.numerical = (pr2-pr)/heps 
  db.analytic= mvn.deriv.margin(a,b,mu,rr,k,1)$deriv
  cat("   numerical: ", db.numerical, ", analytic: ", db.analytic,"\n")
}


cat("derivative wrt rho(j,k): first order approx\n")
for(j in 1:(m-1))
{ for(k in (j+1):m)
  { cat("  (j,k)=", j,k,"\n")
    rr2=rr
    rr2[j,k]=rr[j,k]+heps
    rr2[k,j]=rr[k,j]+heps
    pr2=mvnapp(a,b,mu,rr2)$pr
    drh.numerical = (pr2-pr)/heps 
    drh.analytic= mvn.deriv.rho(a,b,mu,rr2,j,k)$deriv
    cat("   numerical: ", drh.numerical, ", analytic: ", drh.analytic,"\n")
  }
}

cat("\nsecond order approx\n")
pr=mvnapp(a,b,mu,rr,type=2)$pr
cat("pr=mvnapp(avec,bvec,mu=0,sigma=corrmat,type=2)=",pr,"\n")

cat("derivative wrt a_k,b_k, k=1,...,",m,"\n")
for(k in 1:m)
{ cat("  k=", k, " lower\n")
  a2=a
  a2[k]=a[k]+heps
  pr2=mvnapp(a2,b,mu,rr,type=2)$pr
  da.numerical = (pr2-pr)/heps 
  da.analytic= mvn.deriv.margin(a,b,mu,rr,k,-1,type=2)$deriv
  cat("   numerical: ", da.numerical, ", analytic: ", da.analytic,"\n")
  cat("  k=", k, " upper\n")
  b2=b
  b2[k]=b[k]+heps
  pr2=mvnapp(a,b2,mu,rr,type=2)$pr
  db.numerical = (pr2-pr)/heps 
  db.analytic= mvn.deriv.margin(a,b,mu,rr,k,1,type=2)$deriv
  cat("   numerical: ", db.numerical, ", analytic: ", db.analytic,"\n")
}

cat("derivative wrt rho(j,k): second order approx\n")
for(j in 1:(m-1))
{ for(k in (j+1):m)
  { cat("  (j,k)=", j,k,"\n")
    rr2=rr
    rr2[j,k]=rr[j,k]+heps
    rr2[k,j]=rr[k,j]+heps
    pr2=mvnapp(a,b,mu,rr2,type=2)$pr
    drh.numerical = (pr2-pr)/heps 
    drh.analytic= mvn.deriv.rho(a,b,mu,rr2,j,k,type=2)$deriv
    cat("   numerical: ", drh.numerical, ", analytic: ", drh.analytic,"\n")
  }
}
}
\keyword{models}
