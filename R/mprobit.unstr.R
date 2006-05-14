
mprobit.unstr<- function(x,y,id,iprint=0,startpar=0)
{ if(!is.vector(x))
  { if(length(y)!=length(id)|length(id)!=nrow(x)|nrow(x)!=length(y))
     stop("x, y, id not same length")
  }     
  else if(length(y)!=length(id)|length(id)!=length(x)|length(x)!=length(y))
    stop("x, y, id not same length")

  nrec=length(y)
  y=as.integer(y)
  if(sum(y>=0 & y<=1)<nrec) stop("y should be 0 or 1")

  # check constant d
  id2=c(id[1],id[-nrec])
  idif=(id!=id2)
  idif[1]=1
  ncl=sum(idif)
  d=nrec/ncl
  tem=sum(idif[(1:ncl)*d-d+1])
  if(tem!=ncl) stop("nonconstant repeated measures")

  if(is.vector(x)) npred=1
  else npred=ncol(x)
  nclbd=length(unique(id))+5  
  nevals=1000+npred*100
  tol=1.e-3
  ncorp=d*(d-1)/2    # number parameters of the corr matrix
  np=npred+1+ncorp
  # starting point from probit binary regression
  xx=cbind(rep(1,length(y)),x)
  names(xx)[1]<-"intcpt"
  #if(missing(startpar) | length(startpar)!=np | startpar==0)
  if(length(startpar)!=np)
  { th = glm.fit(xx, y, family = binomial(link = "probit"))$coef
    th=c(th,rep(.45,ncorp))
  }
  else { th=startpar }

  if(iprint==1)
  { cat("Initial theta's from glm.fit: \n")
    print(th)
  }

  out <- .C("mprobitunstr",as.double(x),as.integer(y),as.integer(id),
          as.integer(nrec),as.integer(npred),as.integer(d),
          as.integer(nevals),as.integer(iprint),
          as.double(tol), th=as.double(th), nllkval=as.double(1),
          h=as.double(rep(0,np*np)) )
  
  hess=matrix(out$h,np,np)
  tem=out$th[(2+npred):np]
  rhomat=matrix(1,d,d)
  k=0
  for(i in 1:(d-1))
  { for(j in (i+1):d)
    { k=k+1
      rhomat[i,j]=tem[k]
      rhomat[j,i]=tem[k]
    }
  }

  #list(negloglik=out$nllkval, mle=out$th, cov=hess)
  list(negloglik=out$nllkval, mle=out$th, 
      beta=out$th[1:(npred+1)], R=rhomat, cov=hess)
}

