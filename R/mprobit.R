
mprobit<- function(x, y, id, corstr="exch", iprint=0, startpar=0)
{ choices=c("exch","ar","unstr")
  tem=match.arg(corstr,choices)
  if(tem=="ar") out=mprobit.ar(x,y,id,iprint,startpar)
  else if(tem=="unstr") out=mprobit.unstr(x,y,id,iprint,startpar)
  else out=mprobit.exch(x,y,id,iprint,startpar)
  out
}

# formula version that calls mprobit()
mprobit.formula<-function(formula,id,data, corstr="exch", iprint=0,startpar=0)
{ call <- match.call()
  m <- match.call(expand = FALSE)
  m$data=as.data.frame(data)
  m$corstr <- m$iprint <- m$startpar <- m$contrasts <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  Y <- model.extract(m, response)
  X <- model.matrix(Terms, m, contrasts)
  x<-as.matrix(X[,-1])
  y<-c(Y)
  out<-mprobit(x,y,id,corstr,iprint,startpar)
  out
}

