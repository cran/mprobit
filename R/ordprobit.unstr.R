ordprobit.unstr<- function(x,y,id,iprint=0,startpar=0)
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

  # check constant d 
  #tem=apply(as.array(unique(id[seq(1:nrec)[duplicated(id)]])),1,
  #          function(x){sum(ifelse(x==id,1,0))})
  #if(length(unique(tem))!=1)
  #  stop("nonconstant repeated measures")
  #d=sum(ifelse(id[1]==id,1,0))
  # simplify above
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
  tolder=1.e-3
  ncorp=d*(d-1)/2  # number parameters in the corr matrix
  np=npred+norc-1+ncorp;  # norc-1 cutpts
  #ncl=sum(ifelse(unique(id)>0,1,0))
  # starting point from probit binary regression
  #if(missing(startpar) | length(startpar)!=np | startpar==0)
  if(length(startpar)!=np)
  { if(norc==2)
    { ybin=2-y
      xx=cbind(rep(1,length(y)),x)
      names(xx)[1]<-"intcpt"
      th=glm.fit(xx,ybin,family=binomial(link="probit"))$coef
      th=c(th,rep(.4,ncorp))
    }
    else
    { cum=(1:(norc-1))
      cut=rep(0,norc-1)
      for(k in cum)
      { pr=sum(y<=k)
        if (pr==0) pr=1 
        cut[k]=qnorm(pr/nrec)
      }
      #th=c((pout$zeta)/2, (pout$coef)/2, round(runif(ncorp,0.2,0.5),2))
      th=c(cut,rep(0,npred),rep(.4,ncorp))
    }
  }
  else { th=startpar }

  if(iprint==1)
  { cat("Initial theta's: \n")
    print(th)
  }

  out <- .C("ordprobitunstr",as.double(x),as.integer(y),as.integer(id),
     as.integer(nrec),as.integer(npred),as.integer(norc),as.integer(d),
     as.integer(nevals),as.integer(iprint),
     as.double(tolder), th=as.double(th), nllkval=as.double(0),
     h=as.double(rep(0,np*np)) )
  
  hess=matrix(out$h,np,np)
  #diagE=apply(as.array(seq(1,ncol(hess))),1,function(x) ifelse(hess[x,x]==1,1,0) )
  #if(ifelse(sum(diagE)==ncol(hess),1,0))
  # warning("try different initial parameters")
  #list(negloglik=out$nllkval, mle=out$th, cov=hess)
  tem=out$th[(norc+npred):np]
  rhomat=matrix(1,d,d)
  k=0
  for(i in 1:(d-1))
  { for(j in (i+1):d)
    { k=k+1
      rhomat[i,j]=tem[k]
      rhomat[j,i]=tem[k]
    }
  }
      
  list(negloglik=out$nllkval, mle=out$th, cutpts=out$th[1:(norc-1)],
      beta=out$th[norc:(norc+npred-1)], R=rhomat, cov=hess)
}


