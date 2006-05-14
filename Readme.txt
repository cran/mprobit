
This mprobit package has R  functions for multivariate normal rectangle
probabilities (positive exchangeable [any dimension], 
general [up to dimension 6], approximations [up to dimension 19]); MLE of
regression and correlation parameters in the multivariate binary/ordinal
probit models: exchangeable, AR(1), and unstructured correlation matrix.
These multivariate probit models are for repeated measurements or clustered
binary/ordinal data, as the regression parameters are common to different 
components of the multivariate response.
The functions for multivariate normal rectangle probabilities can be
used completely separate from the probit models; they can also be used
to obtain joint probabilities after fitting a multivariate probit model.

References for individual functions are given in the R help files in the
'man' subdirectory. The maximum likelihood estimates and estimated
covariance matrices are obtain using the quasi-Newton method in
  Nash, J.C. (1990). Compact Numerical Methods for Computers: Linear
  Algebra and Function Minimisation,  second edition.  Hilger, New York.
This is the same method as in the function optim() in R.

The C code can also be compiled to get programs to run from the
Unix command line. See the file src/cmakefile.
The C code can be modified for multivariate probit model for a
multivariate binary response with several different binary variables.
In this case there are more regression parameters, and an estimating
equation method like "inference functions for margins" [Chapter 10 of Joe
1997, Multivariate Models and Dependence Concepts, Chapman & Hall]
may be better than maximum.likelihood when there are many parameters
(regression + latent correlation).

There are extensive tests in the subdirectory called tests.
This includes tests of the R functions and the C programs.
See the file tests/run.sh

Version of May 2006.
