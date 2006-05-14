\name{mprobit}
\alias{mprobit}
\alias{mprobit.formula}
\alias{mprobit.exch}
\alias{mprobit.ar}
\alias{mprobit.unstr}
\title{
Maximum Likelihood for Repeated Measures Multivariate Binary Probit Model 
}
\description{
Maximum Likelihood for Repeated Measures Multivariate Binary Probit: Exchangeable,
AR(1) and Unstructured Correlation Matrices. Quasi-Newton minimization
of negative log-likelihood is used with approximation of Joe (1995)
for rectangle probabilities for AR(1) and unstructured correlation,
and one-dimensional Romberg integration for exchangeable correlation. 
}
\usage{
mprobit(x,y,id,corstr="exch",iprint=0,startpar=0)
  or
mprobit.formula(formula,id,data,corstr="exch",iprint=0,startpar=0)
  or
mprobit.exch(x,y,id,iprint=0,startpar=0)
mprobit.ar(x,y,id,iprint=0,startpar=0)
mprobit.unstr(x,y,id,iprint=0,startpar=0)
}
\arguments{
\item{x}{
   vector or matrix of explanatory variables. Each row corresponds
   to an observation and each column to a variable.  The number of rows of x
   should equal the number of data values in y.
   Missing values are not allowed.
}
\item{y}{
   numeric vector containing the binary response (coded with values of 0,1). 
   Missing values are not allowed. 
}
\item{id}{
   group or cluster id, should be a vector of positive integers. 
   If AR(1) correlation, records are assumed to
   be ordered within each cluster id. For unstructured correlation,
   the cluster size must be constant, and
   records are assumed to be ordered the same way with each cluster id
   (i.e., jth record within each cluster refers to a common
   time/condition for the repeated measurement). 
   For AR(1) and exchangeable correlation, cluster size can vary.
   For the formula version, include the data structure.
}
\item{formula}{
  For the formula version of mprobit, a formula expression as for other
  regression models, of the form "response ~ predictors".
}
\item{data}{
  For the formula version of mprobit, the data frame which contains
  the variables in the formula, and the cluster id variable.
}
\item{corstr}{
  For mprobit as a front end to the other three functions,
  corstr="exch" means exchangeable correlation (the default); 
  corstr="ar" means AR(1) correlation;
  corstr="unstr" means unstructured correlation (in which case
  cluster size must be a constant). 
}
\item{iprint}{
  logical indicator, default is FALSE, for whether the iterations for 
  numerical maximum likelihood should be printed.
}
\item{startpar}{
  initial parameter vector in the order: regression coefficients, 
  correlation parameter(s). If not supplied, default = 0, startpar will be 
  generated automatically.
}
}
\value{
  list of MLE of parameters and their associated standard 
  errors, in the order of b0,b1,b2...b(number of covariates), rho(s);
  order of rhos is r12,r13,...,r23,...,r(d-1,d) for unstructured
\item{negloglik}{
  value of negative log-likelihood, evaluated at MLE
}
\item{beta}{
  MLE of regression parameters
}
\item{rho}{
  MLE of latent correlation parameter for AR(1) and exchangeable correlation 
}
\item{rhomat}{
  MLE of matrix of latent correlation parameter for unstructured correlation 
}
\item{mle}{
  MLE of all parameters for unstructured correlation 
}
\item{cov}{
  estimated covariance matrix of the parameters
}
}
\details{
  To get an initial version working, there are constraints:
  (a) For AR(1) and unstructured correlation, the maximum cluster size 
   is 19 (although the joint probabilities
   get to be too small well before this limit is reached);
  (b) The maximum total number of parameters (regression and latent
   correlation parameters combined is 23). So this means a smaller
   upper bound on the number of predictors for the unstructured correlation.

  Also the performance of the quasi-Newton algorithm gets worse
  as the number of parameters increase, particular with the sample
  size (number of clusters) is too small.

  The default starting point used by the code should usually be OK.
  If the returned cov=estimated covariance matrix is the identity matrix,
  then the quasi-Newton iterations did not finish cleanly with a
  gradient vector that is near zero.
  If this case, it would be useful to use iprint=1 to print the
  iterations, and try different starting points in startpar to
  check on the sensitivity to the starting point. The SE estimates
  are better if the quasi-Newton iterations finish with gradient
  vector closer to zero.
}
\references{
  Ashford, JR and Sowden, RR (1970). Multivariate probit analysis.  
     Biometrics, 26, 535-546. 

  Chaganty NR and Joe H (2004). Efficiency of the generalised estimating 
    equations for binary response.  J Royal Statist Soc B, 66, 851-860. 

  Joe, H (1995). Approximations to multivariate  normal rectangle
    probabilities based on conditional expectations, J Amer Stat Assoc, 90, 957-964.
}
\author{H. Joe, Statistics Department, UBC, with assistance of Laing Wei Chou}
\examples{
# first test data set
data(binaryex)
x=binaryex$x
y=binaryex$y
id=binaryex$id

# various ways of using mprobit are shown here
# exchangeable dependence 
out.exch = mprobit.exch(x,y,id)
print(out.exch)

# AR(1) dependence 
out.ar = mprobit(x,y,id,corstr="ar")
print(out.ar)

# unstructured correlation matrix
out.unstr = mprobit.formula(y~x,binaryex$id,data=binaryex,corstr="unstr")
print(out.unstr)

# second test data set
data(binaryar)
dat=binaryar
x=dat$x
y=dat$y
id=dat$id
x2=x*x
dat$x2=x2

# exchangeable dependence 
out2.exch = mprobit.exch(cbind(x,x2),y,id)
print(out2.exch)

# AR(1) dependence 
out2.ar = mprobit.formula(y~x+x2,dat$id,data=dat,corstr="ar")
print(out2.ar)

# varying cluster sizes
set.seed(123)
ncl=nrow(dat)/4
idel=sort(sample(ncl,ncl/4))
idel=idel*4
datsub=dat[-idel,]
out3.ar = mprobit.formula(y~x+x2,datsub$id,data=datsub,corstr="ar")
print(out3.ar)
}
\keyword{models}
