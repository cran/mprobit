options(error = expression(NULL))
library(mprobit)

# 4 repmeas for each cluster
# do a loop of length 2

for(iii in 1:2)
#{ if(iii==1) 
{ if(iii==2) 
  { dat=read.table("binaryex.dat",header=T)
    cat("data file binaryex.dat\n")
  }
  else 
  { dat=read.table("binaryar.dat",header=T)
    cat("data file binaryar.dat\n")
  }
  x=as.vector(dat[,2])
  y=dat[,3]
  id=dat[,1]
  ymat=matrix(y,100,4,T)
  print(cov(ymat))
  print(cor(ymat))
  x2=x*x

  # check exch, AR, unstr
  cat("mprobit.exch\n")
  out.exch <- mprobit.exch(x,y,id,iprint=1)
  print(out.exch)
  cat("============================================================\n")
  cat("mprobit.ar\n")
  out.ar <- mprobit.ar(x,y,id,iprint=1)
  print(out.ar)
  cat("============================================================\n")
  cat("mprobit.unstr\n")
  out.unstr <- mprobit.unstr(x,y,id,iprint=1)
  print(out.unstr)
  cat("============================================================\n")
  cat("mprobit.unstr: different starting point\n")
  param0=c(.4,.9,.5,.7,.7,.6,.8,.4)
  out.unstr <- mprobit.unstr(x,y,id,iprint=1,startpar=param0)
  print(out.unstr)
  cat("============================================================\n")

  # check 2 covariates
  cat("mprobit.exch: x,x^2\n")
  out2.exch <- mprobit.exch(cbind(x,x2),y,id,iprint=1)
  print(out2.exch)
  cat("============================================================\n")
  cat("mprobit.ar: x,x^2\n")
  out2.ar <- mprobit.ar(cbind(x,x2),y,id,iprint=1)
  print(out2.ar)
  cat("============================================================\n")

  # check ordinal probit with 2 categories
  yy=y+1
  cat("ordprobit.exch with 2 categories\n")
  out3.exch <- ordprobit.exch(x,yy,id,iprint=1)
  print(out3.exch)
  # signs are reversed here for beta vector
  cat("============================================================\n")
  cat("ordprobit.ar with 2 categories\n")
  out3.ar <- ordprobit.ar(x,yy,id,iprint=1)
  print(out3.ar)
  cat("============================================================\n")
  cat("ordprobit.unstr with 2 categories\n")
  out3.unstr <- ordprobit.unstr(x,yy,id,iprint=1)
  print(out3.unstr)
  cat("============================================================\n")

  # check exch, AR, with random deletion of one-fourth of 4th measurements
  set.seed(123)
  n=nrow(dat)/4
  idel=sort(sample(n,n/4))
  idel=idel*4

  datsub=dat[-idel,]
  xsub=as.vector(datsub[,3])
  ysub=as.integer(datsub[,4])
  idsub=as.integer(datsub[,1])
  cat("mprobit.exch: subset, unequal cluster sizes\n")
  sub.exch <- mprobit.exch(xsub,ysub,idsub,iprint=1)
  print(sub.exch)
  cat("============================================================\n")
  cat("mprobit.ar: subset, unequal cluster sizes\n")
  sub.ar <- mprobit.ar(xsub,ysub,idsub,iprint=1)
  print(sub.ar)
  cat("============================================================\n")
  # should be error message for unstr
  cat("Testing error messages\n")
  mprobit.unstr(xsub,ysub,idsub,iprint=1)
  cat("============================================================\n")

  # check other error messages
  id[1]=1000
  ordprobit.unstr(x,y,id,iprint=1)
  id[1]=1
  y[1]=2
  mprobit.exch(x,y,id,iprint=1)
  y=y[-1]
  mprobit.exch(x,y,id,iprint=1)
  cat("\n*** end\n\n")
}
