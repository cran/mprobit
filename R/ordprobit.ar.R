ordprobit.ar<- function(x,y,id,iprint=0,startpar=0)
{ if(!is.vector(x))
  { if(length(y)!=length(id)|length(id)!=nrow(x)|nrow(x)!=length(y))
     stop("x, y, id not same length")
  }     
  else if(length(y)!=length(id)|length(id)!=length(x)|length(x)!=length(y))
    stop("x, y, id not same length")

  nrec=length(y)
  y=as.integer(y)
  norc=max(y)
  # check that y values between 1 and norc
  if(sum(y>=1 & y<=norc)<nrec) stop("y should be 1,2,...#categories")

  if(is.vector(x)) npred=1
  else npred=ncol(x)
  nclbd=length(unique(id))+5
  nevals=1000+npred*100  
  tolder=1.e-3
  np=npred+norc  # norc-1 cutpts, rho

  # starting point from probit binary regression
  #if(missing(startpar) | length(startpar)!=np | startpar==0)
  if(length(startpar)!=np)
  { if(norc==2)
    { ybin=2-y
      xx=cbind(rep(1,nrec),x)
      names(xx)[1]<-"intcpt"
      th=glm.fit(xx,ybin,family=binomial(link="probit"))$coef
      th=c(th,.4)     # add corr parameter
    }
    else
    { cum=(1:(norc-1))
      cut=rep(0,norc-1)
      for(k in cum)
      { pr=sum(y<=k)
        if (pr==0) pr=1 
        cut[k]=qnorm(pr/nrec)
      }
      th=c(cut,rep(0,npred),.4) 
    }
  }
  else { th=startpar }

  if(iprint==1)
  { cat("Initial theta's: \n")
    print(th)
  }
  out <- .C("ordprobitar",as.double(x),as.integer(y),as.integer(id),
     as.integer(nrec),as.integer(npred),as.integer(norc),as.integer(nclbd),
     as.integer(nevals),as.integer(iprint),as.double(tolder),
     th=as.double(th), nllkval=as.double(0),h=as.double(rep(0,np*np)) )

  hess=matrix(out$h,np,np)
  #print(out$nllkval)
  #print(out$th)
  #print(sqrt(diag(hess)))
  #print(hess)
  #list(negloglik=out$nllkval, mle=out$th, cov=hess)
  list(negloglik=out$nllkval, cutpts=out$th[1:(norc-1)],
      beta=out$th[norc:(np-1)], rho=out$th[np], cov=hess)
}


