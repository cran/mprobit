options(error = expression(NULL))
library(mprobit)

# 4 repmeas for each cluster, 3 ordinal categories
dat=read.table("ordinalex.dat",header=T)
x=as.vector(dat[,2])
y=dat[,3]
ymat=matrix(y,100,4,T)
id=dat[,1]
x2=x*x

cat("\nordprobit.exch\n")
ord.exch <- ordprobit.exch(x,y,id,iprint=0)
print(ord.exch)
cat("============================================================\n")
cat("ordprobit.univar\n")
ord.univar<- ordprobit.univar(x,y,iprint=1)
print(ord.univar)
cat("============================================================\n")
# use of startpar
startp=c(ord.univar$cutpts,ord.univar$beta,0.5)
cat("ordprobit.exch with startpar\n")
print(startp)
ord.exch2 <- ordprobit.exch(x,y,id,iprint=0,startpar=startp)
print(ord.exch2)
cat("============================================================\n")

# compare ordinal logit or polr in library(MASS)
cat("ordlogit\n")
library(MASS)
ord.polr= polr(as.factor(y)~x)
print(ord.polr)
cat("============================================================\n")

# check 2 covariates
cat("\nordprobit.exch: x,x^2\n")
ord2.exch <- ordprobit.exch(cbind(x,x2),y,id,iprint=0)
print(ord2.exch)
cat("============================================================\n")
cat("ordprobit.univar: x,x^2\n")
ord2.univar<- ordprobit.univar(cbind(x,x2),y,iprint=1)
print(ord2.univar)
cat("============================================================\n")
cat("ordlogit: x,x^2\n")
ord2.polr= polr(as.factor(y)~x+x2)
print(ord2.polr)
cat("============================================================\n")

