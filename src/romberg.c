#include <math.h>
#include <stdio.h>
#include "mprobit.h"

// not getting convergence for some integrals (of derivatives) when rho>=.6?

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
