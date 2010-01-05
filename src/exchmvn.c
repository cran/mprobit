#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mprobit.h"

/* gcc -DMAIN -o exchmvn exchmvn.c pnorms.c romberg.c phi.c -lm */
/* mvn rectangle probability for positive exch case, */

#ifdef MAIN
main()
{ int m,i;
  double rh,*x,*w,pr,eps;
  double exchmvn(int, double *, double *, double, double);
  
  eps=1.e-6;
  scanf("%d", &m);
  while(m>0)
  { scanf("%lf", &rh);
    x=(double *) malloc((m+1) * sizeof(double));
    w=(double *) malloc((m+1) * sizeof(double));
    for(i=1;i<=m;i++)  scanf("%lf", &w[i]);
    for(i=1;i<=m;i++)  scanf("%lf", &x[i]);
    printf("m=%3d, rh=%f\n", m,rh);
    for(i=1;i<=m;i++) printf("%8.4f", w[i]);  printf("\n");
    for(i=1;i<=m;i++) printf("%8.4f", x[i]);  printf("\n");
    pr=exchmvn(m,w,x,rh,eps);
    printf("exch.  : %.10f\n", pr);
    free(x); free(w);
    scanf("%d", &m);
  }
}
#endif

//double exchmvn(int m, double w[], double x[], double rh, double eps)
double exchmvn(int m, double *w, double *x, double rh, double eps)
{ double g(double),pr;
  double romberg(double (*)(double), double, double, double);
  int i;
  extern int mm;
  //extern double ww[M],xx[M],rs,r1;
  extern double *ww,*xx,rs,r1;
  mm=m; rs=sqrt(rh); r1=sqrt(1.-rh);
  xx=(double *) malloc((mm+1) * sizeof(double));
  ww=(double *) malloc((mm+1) * sizeof(double));
  for(i=1;i<=mm;i++) { ww[i]=w[i]; xx[i]=x[i];}
  pr=romberg(g,-UB,UB,eps);
  free(xx); free(ww);
  return(pr);
}

double g(double z)
{ double pnorms(double),phi(double),a,b;
  extern int mm;
  //extern double ww[M],xx[M],rs,r1;
  extern double *ww,*xx,rs,r1;
  int i;
  double tem;
  for(i=1,tem=1.;i<=mm;i++)
  { a=(ww[i]-rs*z)/r1; b=(xx[i]-rs*z)/r1;
    tem*=pnorms(b)-pnorms(a);
  }
  tem*=phi(z);
  return(tem);
}

#ifdef ALONE
double phi(double z)
{ return( 0.3989422804014327*exp(-.5*z*z));}

//#define S 10
// not getting convergence for some integrals (of derivatives) when rho>=.6?
#define S 12

/* Romberg integration (one-dimensional integral) with relative convergence
   criterion, this integration method should be mentioned in any
   text on numerical methods; see also Numerical Recipes.
   Maybe can use the version which return a code for convergence.
*/
double romberg(double (*g)(double), double a, double b, double eps)
{ double t[S+1][S+1],h,fourj,sum,integ;
/* the 0th row and column of t[][] are not used */
  int m,k,i,j;
  h=b-a; m=1;
  /* t[1][1]=h*(g(a)+g(b))/2.;*/ /* either of the two forms work, why? */
  t[1][1]=h*((*g)(a)+(*g)(b))/2.;
  for(k=2;k<=S;k++)
  { h/=2.; m*=2;
    /* for(i=1,sum=0.;i<=m;i+=2) sum+= g(a+i*h);*/
    for(i=1,sum=0.;i<=m;i+=2) sum+= (*g)(a+i*h);
    t[k][1]=t[k-1][1]*.5+sum*h;
    for(j=2,fourj=1.;j<=k;j++)
    { fourj*=4.;
      t[k][j]=t[k][j-1] + (t[k][j-1]-t[k-1][j-1])/(fourj-1.);
    }
    if(fabs((t[k][k]-t[k-1][k-1])/t[k][k]) <= eps)
    { integ=t[k][k]; return(integ);}
  }
  integ=t[S][S];
#ifdef DIAG 
  printf("*** convergence not reached ***\n");
#endif
  return(integ);
}
#endif

#undef UB
//#undef M
