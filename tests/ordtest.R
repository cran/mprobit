options(error = expression(NULL))
library(mprobit)

# 4 repmeas for each cluster, 3 ordinal categories
dat=read.table("ordinalex.dat",header=T)
x=as.vector(dat[,2])
y=dat[,3]
ymat=matrix(y,100,4,T)
print(cov(ymat))
print(cor(ymat))
id=dat[,1]
x2=x*x

# check exch, AR, unstr
cat("ordprobit.exch\n")
ord.exch <- ordprobit.exch(x,y,id,iprint=1)
print(ord.exch)
cat("============================================================\n")
cat("ordprobit.ar\n")
ord.ar <- ordprobit.ar(x,y,id,iprint=1)
print(ord.ar)
cat("============================================================\n")
cat("ordprobit.unstr\n")
ord.unstr <- ordprobit.unstr(x,y,id,iprint=1)
print(ord.unstr)
cat("============================================================\n")

# check 2 covariates
cat("ordprobit.exch: x,x^2\n")
ord2.exch <- ordprobit.exch(cbind(x,x2),y,id,iprint=1)
print(ord2.exch)
cat("============================================================\n")
cat("ordprobit.ar: x,x^2\n")
ord2.ar <- ordprobit.ar(cbind(x,x2),y,id,iprint=1)
print(ord2.ar)
cat("============================================================\n")
cat("ordprobit.unstr: x,x^2\n")
ord2.unstr <- ordprobit.unstr(cbind(x,x2),y,id,iprint=1)
print(ord2.unstr)
cat("============================================================\n")


# check exch, AR, with random deletion of one-fourth of 4th measurements
set.seed(123)
n=nrow(dat)/4
idel=sort(sample(n,n/4))
idel=idel*4

datsub=dat[-idel,]
xsub=as.vector(datsub[,2])
ysub=as.integer(datsub[,3])
idsub=as.integer(datsub[,1])
cat("ordprobit.exch: subset, unequal cluster sizes\n")
osub.exch <- ordprobit.exch(xsub,ysub,idsub,iprint=1)
print(osub.exch)
cat("============================================================\n")
cat("ordprobit.ar: subset, unequal cluster sizes\n")
osub.ar <- ordprobit.ar(xsub,ysub,idsub,iprint=1)
print(osub.ar)
cat("============================================================\n")
cat("Testing error messages\n")
# should be error message for unstr
ordprobit.unstr(xsub,ysub,idsub,iprint=1)

# check other error messages
id[1]=1000
ordprobit.unstr(x,y,id,iprint=1)
id[1]=1
y[1]=0
ordprobit.exch(x,y,id,iprint=1)
y=y[-1]
ordprobit.exch(x,y,id,iprint=1)
