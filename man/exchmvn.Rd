\name{exchmvn}
\alias{exchmvn}
\alias{exchmvn.deriv.margin}
\alias{exchmvn.deriv.rho}
\title{
Exchangeable (positive) multivariate normal
}
\description{
Rectangle probability and derivatives of positive exchangeable multivariate
normal}
\usage{
exchmvn(lb, ub, rh, mu=0, scale=1, eps = 1.e-06)
exchmvn.deriv.margin(lb, ub, rh, k, ksign, eps = 1.e-06)
exchmvn.deriv.rho(lb, ub, rh, eps = 1.e-06)
}
\arguments{
\item{lb}{
   vector of lower limits of integral/probability  
}
\item{ub}{
   vector of upper limits of integral/probability 
}
\item{rh}{
   correlation, rho
}
\item{mu}{
   mean vector
}
\item{scale}{
   standard deviation
}
\item{eps}{
  tolerance for numerical integration
}
\item{k}{
  margin for which derivative is to be taken, that is,
  deriv of exchmvn(lb,ub,rh) with respect to lb[k] or ub[k];
  use exchmvn.deriv.rh for deriv of exchmvn(lb,ub,rh) with respect to rho
}
\item{ksign}{
  =-1 for deriv of exchmvn(lb,ub,rh) with respect to lb[k]
  =+1 for deriv of exchmvn(lb,ub,rh) with respect to ub[k]
}
}
\value{
  rectangle probability or a derivative
}
\references{
  Kotz S and Johnson NL (1972). Continuous Multivariate Distributions. Wiley, New York,
  page 48.
}
\seealso{
}
\author{H. Joe, Statistics Department, UBC}
\examples{
# The tests here show clearly what the function parameters are.
# step size for numerical derivatives (accuracy of exchmvn etc about 1.e-6)
heps = 1.e-4

cat("case 1: m=3\n")
m=3
a=c(-1,-1,-1)
b=c(2,1.5,1)
rh=.6
pr=exchmvn(a,b,rh)
cat("pr=exchmvn(avec,bvec,rh)=",pr,"\n")
cat("derivative wrt rho\n")
rh2=rh+heps
pr2=exchmvn(a,b,rh2)
drh.numerical= (pr2-pr)/heps
drh.analytic= exchmvn.deriv.rho(a,b,rh)
cat("   numerical: ", drh.numerical, ", analytic: ", drh.analytic,"\n")

cat("derivative wrt a_k,b_k, k=1,...,",m,"\n")
for(k in 1:m)
{ cat("  k=", k, " lower\n")
  a2=a
  a2[k]=a[k]+heps
  pr2=exchmvn(a2,b,rh)
  da.numerical = (pr2-pr)/heps 
  da.analytic= exchmvn.deriv.margin(a,b,rh,k,-1)
  cat("   numerical: ", da.numerical, ", analytic: ", da.analytic,"\n")
  cat("  k=", k, " upper\n")
  b2=b
  b2[k]=b[k]+heps
  pr2=exchmvn(a,b2,rh)
  db.numerical = (pr2-pr)/heps 
  db.analytic= exchmvn.deriv.margin(a,b,rh,k,1)
  cat("   numerical: ", db.numerical, ", analytic: ", db.analytic,"\n")
}


cat("\ncase 2: m=5\n")
m=5
a=rep(-1,m)
b=c(2,1.5,1,1.5,2)
rh=.6
pr=exchmvn(a,b,rh)
cat("pr=exchmvn(avec,bvec,rh)=",pr,"\n")
cat("derivative wrt rho\n")
rh2=rh+heps
pr2=exchmvn(a,b,rh2)
drh.numerical= (pr2-pr)/heps
drh.analytic= exchmvn.deriv.rho(a,b,rh)
cat("   numerical: ", drh.numerical, ", analytic: ", drh.analytic,"\n")

cat("derivative wrt a_k,b_k, k=1,...,",m,"\n")
for(k in 1:m)
{ cat("  k=", k, " lower\n")
  a2=a
  a2[k]=a[k]+heps
  pr2=exchmvn(a2,b,rh)
  da.numerical = (pr2-pr)/heps 
  da.analytic= exchmvn.deriv.margin(a,b,rh,k,-1)
  cat("   numerical: ", da.numerical, ", analytic: ", da.analytic,"\n")
  cat("  k=", k, " upper\n")
  b2=b
  b2[k]=b[k]+heps
  pr2=exchmvn(a,b2,rh)
  db.numerical = (pr2-pr)/heps 
  db.analytic= exchmvn.deriv.margin(a,b,rh,k,1)
  cat("   numerical: ", db.numerical, ", analytic: ", db.analytic,"\n")
}

}
\keyword{models}
