\name{ordprobit.univar}
\alias{ordprobit.univar}
\title{
Maximum Likelihood for Ordinal Probit Model
}
\description{
Maximum Likelihood for Ordinal Probit: 
Newton-Raphson minimization of negative log-likelihood. 
}
\usage{
ordprobit.univar(x,y,iprint=0,maxiter=20,toler=1.e-6)
}
\arguments{
\item{x}{
   vector or matrix of explanatory variables. Each row corresponds
   to an observation and each column to a variable.  The number of rows of x
   should equal the number of data values in y, and there should be fewer 
   columns than rows.  Missing values are not allowed.
}
\item{y}{
   numeric vector containing the ordinal response. 
   The values must be in the range 1,2,..., number of categories.
   Missing values are not allowed. 
}
\item{iprint}{
  logical indicator, default is FALSE, for whether the iterations for 
  numerical maximum likelihood should be printed.
}
\item{maxiter}{
  maximum number of Newton-Raphson iterations, default = 20.
}
\item{toler}{
  tolerance for convergence in Newton-Raphson iterations, default = 1.e-6.
}
}
\value{
  list of MLE of parameters and their associated standard 
  errors, in the order cutpt1,...,cutpt(number of categ-1),b1,...b(number of covariates).
\item{negloglik}{
  value of negative log-likelihood, evaluated at MLE
}
\item{cutpts}{
  MLE of ordered cutpoint parameters
}
\item{beta}{
  MLE of regression parameters
}
\item{cov}{
  estimated covariance matrix of the parameters
}
}
\details{
  If ordprobit for repeated measures ordinal probit fails to converge
  from the simple starting point in that function, this function
  ordprobit.univar should provide a better starting point.
  It is also equivalent to ordprobit with an identity latent
  correlation matrix.

  The ordinal probit model is similar to the ordinal logit model
  (proportion odds logistic regression : polr in library MASS),
  The parameter estimate of ordinal logit are roughly 1.8 to 2 times
  those of ordinal probit (the signs of the parameters in polr
  may be different, as this function may be using a different
  orientation for the latent variable. 
}
\references{
  Anderson, J.A. and Pemberton, J.D. (1985). The grouped continuous model for
    multivariate ordered categorical variables and covariate adjustment.  
    Biometrics,  41, 875-885.

}
\examples{
data(ordinalex)
x=as.vector(ordinalex$x)
y=ordinalex$y
ord.univar = ordprobit.univar(x,y)
print(ord.univar)
startp=c(ord.univar$cutpts,ord.univar$beta,0.5)
ord.exch <- ordprobit.exch(x,y,ordinalex$id,iprint=0,startpar=startp)
print(ord.exch)
}
\keyword{models}
