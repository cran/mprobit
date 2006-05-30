#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mprobit.h"


void mvn1(int m, double w[], double x[], double r[][M], int nsim, double eps, 
     double *pr, double *sd, int *ifail)
/* first order approximation for conditional probabilities */
/* ref. Joe (1993). Approximations to Multivariate Normal Rectangle
Probabilities Based on Conditional Expectations*/
/*  INPUT
   m = dimension of multivariate normal probability
   w = vector of lower bounds
   x = vector of upper bounds
   r = correlation matrix (1's on the diagonal)
   nsim = number of random permutations used 
      [ nsim=0 for enumerating all permutations, otherwise random perms 
        nsim=0 recommended for m<=6 ]
   eps = error bound for 2-D bivariate probabilities (1.e-6 recommended)
    OUTPUT
   pr = probability of rectangle: P(w_i<X_i<z_i, i=1,...,n)
   sd = measurement of accuracy (SD from the permuted conditional prob's)
   ifail = code for whether prob was succesfully computed 
     (ifail =0 for success, ifail >0 for failure)
    ROUTINES CALLED (directly or indirectly):
   mulnor (Schervish's MVN prob)
   bivnor (Donnelly's bivariate normal upper quadrant prob)
   pnorms (univariate standard normal cdf)
   approx1, cond1 (approx to conditional probs)
   nexper (generate permutations systematically)
   isamp  (generate random permutations)
   gepp   (linear system solver, Gaussian elim with partial pivoting)
*/
/* w[0], x[0], r[0][], r[][0] not used here */
/* No checks are done on r, eg., positive definiteness not checked,
   an error may or may not occur if r is not a corr. matrix (user-beware) */
 /*   for the positive equicorrelated case, use 
        the 1-dimensional integral given on page 193 of Tong's book
        on the multivariate normal distribution */

{  double lb[8],ub[8],bound,eps1,s[22],tem;
   double p[M],pp[M][M],sc,sc2,rect;
   int i,ifault,inf[8],m1,j,ii,is,even,a[M],idet,isim,k,mm;
   int typ[M];
   void isamp(int *, int);
   double pnorms(double);
   void mulnor(double [], double [], double [], double, int,
               int [], double *, double *, int *);
   void nexper(int, int [], int *, int *);
   void approx1(double [], double [], double [][M], int, int [], 
     double [], double [][M], double *, int *);

   *ifail=0;
   for(i=1;i<=m;i++) { if(w[i]>=x[i]) { *ifail=1; *pr=0.; *sd=0.; return;}}
   if(m==1) { *pr=pnorms(x[1])-pnorms(w[1]); *sd=1.e-6; return;}
   if(m==2) 
   {  inf[0]=2; inf[1]=2;  eps1=1.e-6;
      lb[0]=w[1]; lb[1]=w[2]; ub[0]=x[1]; ub[1]=x[2]; s[0]=r[1][2];
      mulnor(ub,lb,s,eps1,m,inf,&rect,&bound,&ifault);
      if(ifault!=0) printf("error in mulnor %d\n", ifault);
      *pr=rect; *sd=bound; *ifail=ifault; return;
   }

  /* m>=3 */
/* some changes for x_i or w_i = oo, May 13, 1994 */
/* for(i=1;i<=m;i++)
   { p[i]=pnorms(x[i])-pnorms(w[i]);}
   m1=2; inf[0]=2; inf[1]=2;*/
   for(i=1;i<=m;i++)
   { if (w[i]<=-5.) { typ[i]=1; p[i]=pnorms(x[i]);}
     else if(x[i]>=5.) {typ[i]=0; p[i]=1.-pnorms(w[i]);}
     else { typ[i]=2; p[i]=pnorms(x[i])-pnorms(w[i]);}
   }
   m1=2; 
   for(i=1;i<=m;i++) pp[i][i]=p[i]*(1.-p[i]);
   for(i=1;i<m;i++)
   {  for(j=i+1;j<=m;j++)
      {  lb[0]=w[i]; lb[1]=w[j]; ub[0]=x[i]; ub[1]=x[j];
         s[0]=r[i][j];
         inf[0]=typ[i]; inf[1]=typ[j];
         mulnor(ub,lb,s,eps,m1,inf,&tem,&bound,&ifault);
         if(ifault!=0) 
         { printf("error in mulnor %d\n", ifault);
           *pr=0.; *sd=bound; *ifail=ifault; return;
         }
         pp[i][j]=tem-p[j]*p[i];
         pp[j][i]=pp[i][j];
      }
   }
   if(nsim==0) { /*enumerate all permutations */
      is=0; ii=0; sc=0.; sc2=0.;
      for(i=2,mm=1;i<=m;i++) mm*=i;  
      for(k=1;k<=mm;k++)
    /*do*/ { nexper(m,a,&is,&even);
          if(a[1]<a[2]) 
          { approx1(w,x,r,m,a,p,pp,&rect,&idet);
            if(idet==1) { ii++; sc+=rect; sc2+=rect*rect;
                /*printf("%8.4f", rect);*/
                   }
            // HJ check May 11, 2004
            //printf("%d %8.4f\n", k,rect);
          }
      } /*while(is==1);*/
    /*printf("\n");*/
   }
   else /* random permutations */
   { ii=0; sc=0.; sc2=0.;
     for(isim=1;isim<=nsim;isim++)
     { isamp(a,m);
       approx1(w,x,r,m,a,p,pp,&rect,&idet);
       if(idet==1) { ii++; sc+=rect; sc2+=rect*rect;}
     } 
   }
   sc/=ii; sc2=(sc2-ii*sc*sc)/(ii-1.); 
   *pr=sc; if(sc2>0.) *sd=sqrt(sc2); else *sd=0.;
/* printf("  mean approx=%9.5f %9.5f %3d\n\n",sc,*sd,ii);*/
}

void approx1(double w[], double x[], double r[][M], int m, int a[], 
     double p[], double pp[][M], double *rec, int *idet)
{  double pra,rec2;
   double lb[8],ub[8],eps,bound,s[22];
   int ii,a1,a2,idt,ifault,m1,inf[8];
   void mulnor(double [], double [], double [], double, int,
               int [], double *, double *, int *);
   double cond1(int, int [], double [], double [][M], int *);

   /* approx , biv prob and product of cond prob*/
   a1=a[1]; a2=a[2];
   m1=2; inf[0]=2; inf[1]=2;  eps=1.e-6;
   lb[0]=w[a1]; lb[1]=w[a2]; ub[0]=x[a1]; ub[1]=x[a2];
   s[0]=r[a1][a2];
   mulnor(ub,lb,s,eps,m1,inf,&rec2,&bound,&ifault);
   if(ifault!=0) printf("error in mulnor %d\n", ifault);
   *idet=1;
   for(ii=3;ii<=m;ii++)
   { pra=cond1(ii,a,p,pp,&idt);
     rec2*=pra;
     if(idt==0) *idet=0;
   }
   *rec=rec2;
}

double cond1(int m, int a[], double p[], double pp[][M], int *idt)
{  double A[N][N];
   double det,tem,tol;
   int i,j;
   void gepp(double [][N], int, int, double *, double);

   //tol=1.e-4;
   // HJ change May 11, 2004
   tol=1.e-8; 
   for(i=1;i<m;i++)
   {  A[i][m]=pp[a[i]][a[m]];
      A[i][i]=pp[a[i]][a[i]];
   }
   for(i=1;i<m-1;i++)
   {  for(j=i+1;j<m;j++)
      {  A[i][j]=pp[a[i]][a[j]];
         A[j][i]=A[i][j];
      }
   }
   *idt=1;
   gepp(A,m-1,m,&det,tol);
   if(fabs(det)==0.) { *idt=0;}
   for(i=1,tem=0.;i<m;i++) tem+=A[i][m]*(1.-p[a[i]]);
   tem+=p[a[m]]; if(tem>1.) tem=1;  if(tem<0.) tem=0.;
   return(tem);
}

void gepp(double A[][N], int n, int nk, double *det, double tol)
/* routine for solving a linear system of equations */
{   int parity,i,j,k,l,mrow;
    double mx,t,sum;
    parity=1;
    for(k=1;k<n;k++)
    {   mx=fabs(A[k][k]); mrow=k;
	for(l=k+1;l<=n;l++)
	{  if(fabs(A[l][k])>mx) {mx=fabs(A[l][k]); mrow=l;} }
	if(mx<=tol) {*det=0.; return;};
	if(mrow>k) {parity*= -1;
		    for(i=1;i<=nk;i++)
		      {t=A[mrow][i];A[mrow][i]=A[k][i];A[k][i]=t;}
                   }
        for(i=k+1;i<=n;i++)
	{      A[i][k]/=A[k][k];
	       for(j=k+1;j<=nk;j++) A[i][j]-=A[i][k]*A[k][j];
        }
    }
    if(fabs(A[n][n])<=tol) {*det=0.; return;}
    for(i=1,*det=parity;i<=n;i++) *det*=A[i][i];
    if(n==nk) return;
    for(l=n+1;l<=nk;l++)
    {   A[n][l]/=A[n][n];
	for(i=n-1;i>=1;i--)
	{   for(j=i+1,sum=0.;j<=n;j++) sum+=A[i][j]*A[j][l];
	    A[i][l]=(A[i][l]-sum)/A[i][i];
        }
    }
}

void mvn2(int m, double w[], double x[], double r[][M], int nsim, double eps, 
     double *pr, double *sd, int *ifail)
/*     int igenz, double *pr, double *sd, int *ifail)
   igenz = 1 if genz's method used for 4-dim MVN prob,
             otherwise schervish's routine is used.
   (code allowing igenz=1 has been commented out, 
            genz's method is slower for 4D)
*/
/* second order approximation for conditional probabilities */
/* ref. Joe (1993). Approximations to Multivariate Normal Rectangle
Probabilities Based on Conditional Expectations*/
/*  INPUT
   m = dimension of multivariate normal probability
   w = vector of lower bounds
   x = vector of upper bounds
   r = correlation matrix (1's on the diagonal)
   nsim = number of random permutations used 
      [ nsim=0 for enumerating all permutations, otherwise random perms
        nsim=0 recommended for m<=6 ]
   eps = error bound for 4-D bivariate probabilities (1.e-4 or 1.e-5
       recommended depending on m)
    OUTPUT
   pr = probability of rectangle: P(w_i<X_i<z_i, i=1,...,n)
   sd = measurement of accuracy (SD from the permuted conditional prob's)
   ifail = code for whether prob was succesfully computed
     (ifail =0 for success, ifail >0 for failure)
    ROUTINES CALLED (directly or indirectly):
   mulnor (Schervish's MVN prob)
   bivnor (Donnelly's bivariate normal upper quadrant prob)
   pnorms (univariate standard normal cdf)
   approx2, cond2 (approx to conditional probs)
   iorder (auxiliary routine)
   nexper (generate permutations systematically)
   isamp  (generate random permutations)
   gepp   (linear system solver, Gaussian elim with partial pivoting)
*/
/* w[0], x[0], r[0][], r[][0] not used here */
/* No checks are done on r, eg., positive definiteness not checked,
   an error may or may not occur if r is not a corr. matrix (user-beware) */
/*     for the positive equicorrelated case, use 
        the 1-dimensional integral given on page 193 of Tong's book
        on the multivariate normal distribution  */

{  double rect;
   double lb[8],ub[8],eps1,bound,s[22];
   double p[M],pp[M][M],sc,sc2;
   double p3[M][M][M],p4[M][M][M][M],q4[M][M][M][M],tem;
   double rr[M][M];
   int i,ifault,inf[8],m1,k,j,ii,is,even,a[M],idet,isim;
   int m2,mm,l,i1,j1,k1,l1;
   int typ[M];
   void isamp(int *, int);
   double pnorms(double);
   void mulnor(double [], double [], double [], double, int,
               int [], double *, double *, int *);
   void iorder(int, int, int, int, int *, int *, int *, int *);
   void nexper(int, int [], int *, int *);
   void approx2(int, int [], 
     double [], double [][M], double [][M][M], double [][M][M][M], 
     double *, int *);

   *ifail=0;
   for(i=1;i<=m;i++) { if(w[i]>=x[i]) { *ifail=1; *pr=0.; *sd=0.; return;}}
   if(m==1) { *pr=pnorms(x[1])-pnorms(w[1]); *sd=1.e-6; return;}
   if(m==2) 
   {  inf[0]=2; inf[1]=2;  eps1=1.e-6;
      lb[0]=w[1]; lb[1]=w[2]; ub[0]=x[1]; ub[1]=x[2]; s[0]=r[1][2];
      mulnor(ub,lb,s,eps1,m,inf,&rect,&bound,&ifault);
      if(ifault!=0) printf("error in mulnor %d\n", ifault);
      *pr=rect; *sd=bound; *ifail=ifault; return;
   }
   if(m==3 || m==4)
   {  eps1=1.e-6;  /* changed on Apr 21, 1994*/ eps1=1.e-5;
      for(i=0;i<m;i++) {inf[i]=2; lb[i]=w[i+1]; ub[i]=x[i+1];}
      for(i=2,k=0;i<=m;i++) 
      {  for(j=1;j<i;j++) { s[k]=r[i][j]; k++;}}
   /* s[0]=r[1][2]; s[1]=r[1][3]; s[2]=r[2][3];*/
   /* addition on Apr 24, 1994 */
      for(i=0;i<m;i++) 
      { if(ub[i]>=5.) inf[i]=0;
        if(lb[i]<=-5.) inf[i]=1;
      }
      mulnor(ub,lb,s,eps1,m,inf,&rect,&bound,&ifault);
      if(ifault!=0) printf("error in mulnor %d\n", ifault);
      *pr=rect; *sd=bound; *ifail=ifault; return;  
   }

   /* m>=5 */
   eps1=1.e-6; m1=2;
   for(i=1;i<=4;i++) rr[i][i]=1.;
   for(i=1;i<=m;i++)
   { p[i]=pnorms(x[i])-pnorms(w[i]); pp[i][i]=p[i];}
   for(i=0;i<4;i++) inf[i]=2;
   for(i=1;i<m;i++)
   {  for(j=i+1;j<=m;j++)
      {  lb[0]=w[i]; lb[1]=w[j]; ub[0]=x[i]; ub[1]=x[j];
         s[0]=r[i][j];
         mulnor(ub,lb,s,eps1,m1,inf,&tem,&bound,&ifault);
         if(ifault!=0) printf("error in mulnor %d\n", ifault);
         pp[i][j]=tem;
         pp[j][i]=pp[i][j];
      }
   }
   /* 3 and 4 dimensional arrays (of moments) */
   eps1=1.e-5; m1=3; m2=4;
/* addition on Apr 25, 1994 */
   for(i=1;i<=m;i++)
   { if(x[i]>=5.) typ[i]=0;
     else if(w[i]<=-5.) typ[i]=1;
     else typ[i]=2;
   }
   for(i=1;i<=m;i++)
   {  for(j=1;j<m;j++)
      {  for(k=j+1;k<=m;k++)
         {  if(i==j || i==k) { p3[i][j][k]=pp[j][k]; p3[i][k][j]=pp[j][k];}
            else
            {  ub[0]=x[i]; ub[1]=x[j]; ub[2]=x[k];
               lb[0]=w[i]; lb[1]=w[j]; lb[2]=w[k];
/* addition on Apr 25, 1994 */
               inf[0]=typ[i]; inf[1]=typ[j]; inf[2]=typ[k];
               s[0]=r[i][j]; s[1]=r[i][k]; s[2]=r[j][k];
               mulnor(ub,lb,s,eps1,m1,inf,&tem,&bound,&ifault);
               if(ifault!=0) printf("error in mulnor %d\n", ifault);
               p3[i][j][k]=tem; p3[i][k][j]=tem;
            }
         }
      }
   }

   for(i=1;i<m;i++)
   {  for(j=i+1;j<=m;j++)
      {  for(k=1;k<m;k++)
         {  for(l=k+1;l<=m;l++)
            {  if(i==k && j==l) p4[i][j][k][l]=pp[i][j]*(1.-pp[i][j]);
               else if (i==k) p4[i][j][k][l]=p3[l][i][j]-pp[i][j]*pp[k][l];
               else if (j==l) p4[i][j][k][l]=p3[k][i][j]-pp[i][j]*pp[k][l];
               else if (j==k) p4[i][j][k][l]=p3[l][i][j]-pp[i][j]*pp[k][l];
               else if (i==l) p4[i][j][k][l]=p3[k][i][j]-pp[i][j]*pp[k][l];
            /* else if (m==4) p4[i][j][k][l]=0.;*/
               else if(i<j && j<k && k<l)
               {  /* schervish's routine */
                  ub[0]=x[i]; ub[1]=x[j]; ub[2]=x[k]; ub[3]=x[l];
                  lb[0]=w[i]; lb[1]=w[j]; lb[2]=w[k]; lb[3]=w[l];
                  s[0]=r[i][j]; s[1]=r[i][k]; s[2]=r[j][k];
                  s[3]=r[i][l]; s[4]=r[j][l]; s[5]=r[k][l];
/* addition on Apr 25, 1994 */
                  inf[0]=typ[i]; inf[1]=typ[j]; inf[2]=typ[k]; inf[3]=typ[l];
                  mulnor(ub,lb,s,eps,m2,inf,&tem,&bound,&ifault);
                  if(ifault!=0) printf("error in mulnor %d\n", ifault);
                  q4[i][j][k][l]=tem;
                  p4[i][j][k][l]=tem-pp[i][j]*pp[k][l];
               }
               else
               {  iorder(i,j,k,l,&i1,&j1,&k1,&l1);
                  p4[i][j][k][l]=q4[i1][j1][k1][l1]-pp[i][j]*pp[k][l];
               }
               p4[j][i][k][l]=p4[i][j][k][l];
               p4[i][j][l][k]=p4[i][j][k][l];
               p4[j][i][l][k]=p4[i][j][k][l];
            }
         }
      }
   }
   for(i=1;i<=m;i++)
   {  for(j=1;j<m;j++)
      {  for(k=j+1;k<=m;k++)
         { p3[i][j][k]-=p[i]*pp[j][k]; p3[i][k][j]=p3[i][j][k]; 
         }
      }
   }
   if(nsim==0) { /*enumerate all permutations */
      for(i=2,mm=1;i<=m;i++) mm*=i;  
      is=0; ii=0; sc=0.; sc2=0.;
      for(k=1;k<=mm;k++)
      {  nexper(m,a,&is,&even);
         if(a[1]<a[2] && a[2]<a[3] && a[3]<a[4])
         { /*approx2(w,x,r,m,a,p,pp,p3,p4,&rect,&idet);*/
           approx2(m,a,p,pp,p3,p4,&rect,&idet);
           if(idet==1) { ii++; sc+=rect; sc2+=rect*rect;
          /* printf("%8.4f", rect);*/
              }
         }
     }/* while(is==1);*/
   /*printf("\n");*/
   }
   else /* random permutations */
   { ii=0; sc=0.; sc2=0.;
     for(isim=1;isim<=nsim;isim++)
     { isamp(a,m);
       approx2(m,a,p,pp,p3,p4,&rect,&idet);
       if(idet==1) { ii++; sc+=rect; sc2+=rect*rect;}
     }
   }
   sc/=ii; sc2=(sc2-ii*sc*sc)/(ii-1.); 
   *pr=sc; if(sc2>0.) *sd=sqrt(sc2); else *sd=0.;
/* printf("  mean approx=%9.5f %9.5f %3d\n\n",sc,*sd,ii);*/
}

/*void approx2(double w[], double x[], double r[][M], int m, int a[], */
void approx2(int m, int a[], 
     double p[], double pp[][M], double p3[][M][M], double p4[][M][M][M], 
     double *rec, int *idet)
{  double pra,rec2;
   double cond2(int, int [], double [], double [][M], double [][M][M], 
             double [][M][M][M], int *);
   int ii,a1,a2,a3,idt,m1,a4;

   /* approx , triv prob and product of cond prob*/
   a1=a[1]; a2=a[2]; a3=a[3]; a4=a[4];
   m1=3; m1=4;
   rec2=p4[a1][a2][a3][a4]+pp[a1][a2]*pp[a3][a4]; 
   *idet=1;
   for(ii=5;ii<=m;ii++)
   { pra=cond2(ii,a,p,pp,p3,p4,&idt);
     rec2*=pra;
     if(idt==0) *idet=0;
   }
   *rec=rec2;
}

double cond2(int m, int a[], double p[], double pp[][M], double p3[][M][M], 
             double p4[][M][M][M], int *idt)
{  double A[N][N];
   double det,tem;
   int i,j,mm,ii,jj,k,l;

   mm=(m*(m-1))/2;
   /* setting up linear system */
   for(i=1;i<m;i++)
   {  for(j=1;j<m;j++) A[i][j]=pp[a[i]][a[j]]-p[a[i]]*p[a[j]];
      A[i][mm+1]=pp[a[i]][a[m]]-p[a[i]]*p[a[m]];;
   }
   for(i=1;i<m;i++)
   {  for(j=1,jj=m-1;j<m-1;j++)
      {  for(k=j+1;k<m;k++)
         {  jj++; A[i][jj]=p3[a[i]][a[j]][a[k]]; A[jj][i]=A[i][jj]; }
      }
   }

   for(i=1,ii=m-1;i<m-1;i++)
   {  for(j=i+1;j<m;j++)
      {  ii++;
         for(k=1,jj=m-1;k<m-1;k++)
         {  for(l=k+1;l<m;l++)
            {  jj++; A[ii][jj]=p4[a[i]][a[j]][a[k]][a[l]];}
         }
         A[ii][mm+1]=p3[a[m]][a[i]][a[j]];
      }
   }

   *idt=1;
   gepp(A,mm,mm+1,&det,1.e-8);
   if(fabs(det)==0.) { *idt=0;}
   for(i=1,tem=0.;i<m;i++) tem+=A[i][mm+1]*(1.-p[a[i]]);
   for(i=1,ii=m-1;i<m-1;i++)
   {  for(j=i+1;j<m;j++)
      {  ii++; tem+=A[ii][mm+1]*(1.-pp[a[i]][a[j]]); }
   }
   tem+=p[a[m]]; if(tem>1.) tem=1;  if(tem<0.) tem=0.;
   return(tem);
}

void iorder(int i, int j, int k, int l, int *i1, int *j1, int *k1, int *l1)
{  int y[5],tem,ii,jj,kk;
   y[1]=i; y[2]=j; y[3]=k; y[4]=l;
   for(jj=2;jj<=4;jj++)
   { for(ii=1;ii<=jj-1;ii++)
     { if (y[ii]>y[jj]) {
        tem=y[jj];
        for(kk=jj;kk>=ii+1;kk--) y[kk]=y[kk-1];
        y[ii]=tem; }
     }
   }
   *i1=y[1]; *j1=y[2]; *k1=y[3]; *l1=y[4];
}

/* pnorms replaced by code from R source, which is more accurate: March 1998 */


void nexper(int n, int a[], int *mtc, int *even)
/*  routine for generating the permutations of 1,...,n in systematic
          order, a has length n */
{ int s,d,i,j,nm3,j1,j2,ia,l,i1,m;
  if(*mtc==0)
  {  nm3=n-3; for(i=1;i<=n;i++) a[i]=i;
     *mtc=1; *even=1; return;
  }
  if(*even == 1) 
  {  ia=a[1]; a[1]=a[2]; a[2]=ia; *even=0; 
     if(a[n] != 1 || a[1] != 2+(n%2)) return;
     if(n<=3) {*mtc=0; return;}
     for(i=1;i<=n-3;i++)
     { if(a[i+1] != a[i]+1) return; }
     *mtc=0; return;
  }
  else
  {  s=0;
     for(i1=2;i1<=n;i1++)
     { ia=a[i1]; i=i1-1;
       for(j=1,d=0;j<=i;j++) {  if(a[j] > ia) d++;}
       s+=d;
       if(d != i*(s%2)) goto A;
     }
     *mtc=0; return;
  A: m=((s+1)%2)*(n+1);
     for(j=1;j<=i;j++)
     { j1=a[j]-ia; j2=a[j]-m;
       if((j1 >= 0 && j2 < 0) || (j1<0 && j2>=0) )
       { m=a[j]; l=j;}
     }
     a[l]=ia; a[i1]=m; *even=1; return;
  }
}

void isamp(int x[], int n)
/* random permutation of size n without replacement */
{ int j,k,temp,rndgen(int, int);
  for(j=1;j<=n;j++) x[j]=j;
  for(j=1; j<=n; j++)
  {  k=rndgen(j,n);	/*uniform on [j,n]*/
     temp=x[j]; x[j]=x[k]; x[k]=temp;
  }
}

int rndgen(int a, int b)
/* random integer between a and b inclusive */
{ double r;
/* random() is built-in random number generator with numbers from 0 to
    2147483647 */
  // r=random()/2147483648.0;
  // mingw for Windows doesn't have random()
  r=rand()/2147483648.0;
  return((int)(ceil((b-a+1)*r))+a-1);
}

