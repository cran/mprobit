mprobit.exch<- function(x, y, id, iprint=0, startpar=0)
{ if(!is.vector(x))
  { if(length(y)!=length(id)|length(id)!=nrow(x)|nrow(x)!=length(y))
     stop("x, y, id not same length")
  }     
  else if(length(y)!=length(id)|length(id)!=length(x)|length(x)!=length(y))
    stop("x, y, id not same length")

  nrec=length(y)
  y=as.integer(y)
  if(sum(y>=0 & y<=1)<nrec) stop("y should be 0 or 1")

  if(is.vector(x)) npred=1
  else npred=ncol(x)
  nclbd=length(unique(id))+5
  nevals=1000+npred*100  # increasing in npred
  tolder=1.e-4
  np=npred+2   # intercept + rho

  # starting point from probit binary regression
  xx=cbind(rep(1,length(y)),x)
  names(xx)[1]<-"intcpt"
  if(length(startpar)!=np)
  { th = glm.fit(xx, y, family = binomial(link = "probit"))$coef
    th = c(th, 0.4)  # add corr parameter
  }
  else { th=startpar }
  if(iprint==1) 
  { cat("Initial theta's from glm.fit: \n")
    print(th)
  }
  out <- .C("mprobitexch",as.double(x),as.integer(y),as.integer(id),
     as.integer(nrec),as.integer(npred),as.integer(nclbd),
     as.integer(nevals),as.integer(iprint),as.double(tolder),
     th=as.double(th),
     nllkval=as.double(0),h=as.double(rep(0,np*np)) )
  hess=matrix(out$h,np,np)
  #print(out$nllkval)
  #print(out$th)
  #print(sqrt(diag(hess)))
  #print(hess)
  #list(negloglik=out$nllkval, mle=out$th, cov=hess)
  list(negloglik=out$nllkval, beta=out$th[1:(np-1)], rho=out$th[np], cov=hess)
}


