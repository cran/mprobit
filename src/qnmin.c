#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include "mprobit.h"

/* quasi-newton minimizer from Nash, "Compact Numerical
   Methods for Computers" (1979), alg 21.
*/
/* malloc version */
#ifdef NOALLOC
void qnmin(int np, double b[], double h[][MAXP], double *p0, int nevals, 
  int *ifail, int mon, void (*funct1)(int, double *, double *), double tolder)
#else
void qnmin(int np, double b[], double **h, double *p0, int nevals, 
  int *ifail, int mon, void (*funct1)(int, double *, double *), double tolder)
#endif
/*
   funct1(np,b,p0) : returns value p0 at parameter vector b of dimension np,
   grad(np,b,g,p0,tolder) : gradient function, follows funct1 evaluation
     puts gradient vector at b in g; tolder is stepsize for
     numerical derivative

   Inputs:

   funct1 = function to be minimized
   np = dimension of function to be minimized
   b = initial guess for point of minimum
   nevals = maximum number of function evals. permitted
   mon = 1 for monitoring (printing of iterations)
       = 0 for no printing
   tolder = stepsize for numerical derivative
 
   Outputs:
   b = point of minimum
   h = estimated hessian matrix of second derivatives
     (actually inverse of Hessian)
     it is the estimate of covariance matrix of the MLE,,
     when funct1 is a negative log-likelihood
   p0 = value of the function minimum
   ifail = returned failure indicator 
       0 means ok
       1 if np>MAXP (dimension too large)
       2 for too many function evaluations
       3 for initial infeasible point
 
   h is allocated as (np+1)*(np+1) matrix if not -DNOALLOC on compiling 
 
   If the initial function evaluation is larger than 
     1.0e308 the program exits at once - this should be
     used as an indicator of an infeasible point.        
 
   Note that 0 index of arrays not used
*/

{ double k,w,tol,eps,p,d1,sn,s,d2;
  int count,i,j,ig,ifn;
  void grad(int, double *, double *, double,
      void (*funct1)(int, double *, double *), double);
#ifdef NOALLOC
  double x[MAXP], c[MAXP], g[MAXP], t[MAXP];
  if (np<=0 || np>MAXP-1) 
  { *ifail = 1; printf(" np out of range\n"); goto R;}
  
#else
  double *x,*c,*g,*t;
  /* allocate space */
  x=(double *) malloc((np+1) * sizeof(double));
  c=(double *) malloc((np+1) * sizeof(double));
  g=(double *) malloc((np+1) * sizeof(double));
  t=(double *) malloc((np+1) * sizeof(double));
  
  if (np<=0) { *ifail = 1; printf(" np out of range\n"); goto R;}
#endif
  
  w=.2; tol=1.e-4; /*eps=1.e-6;*/ eps=1.e-5;
  ifn = np+1; ig = 1;
  (*funct1)(np,b,p0);
  if(mon==1)
  { printf("curr val %8.4f at ", *p0);
    for(i=1;i<=np;i++) printf(" %8.4f", b[i]);
    printf("\n");
  }
  if(*p0>1.0e308)
  { *ifail = 3; printf("**Initial point infeasible\n"); goto R;}
  grad(np,b,g,*p0,funct1,tolder);
  if(mon==1)
  { printf("grad vec:");
    for(i=1;i<=np;i++) printf(" %11.4f", g[i]);
    printf("\n");
  }
  
  /* reset hessian */
  
  m10: for(i=1;i<=np;i++)
  { for(j=1;j<=np;j++) h[i][j]=0.0;
    h[i][i] = 1.0;
  }
  /* ilast not used */
  //ilast = ig;
  
  /* top of iteration  */
  
  m40: for(i=1;i<=np;i++) { x[i] = b[i]; c[i] = g[i]; }
  
  /* find search direction t */
  
  for(i=1,d1=0.,sn=0.;i<=np;i++)
  { for(j=1,s=0.;j<=np;j++) s-=h[i][j]*g[j];
    t[i]=s; sn+=s*s; d1-=s*g[i];
  }
  
  /* check if downhill */
  //if(d1<=0.0) goto m10;
  /* change to avoid infinite loop if g[]=0 */
  if (d1<0.0) goto m10;
  
  /* search along t */
  
  sn = 0.5/sqrt(sn);
  k=1.; if(sn<k) k=sn; k/=w;
  p= *p0+5.;
  while (p>= *p0-d1*k*tol)
  { k*=w;
    for(i=1,count=0;i<=np;i++)
    { b[i]=x[i]+k*t[i];
      if(fabs(b[i]-x[i])<eps) count+=1;
    }
    
  /* check if converged, ifail=0 means successful conclusion */
    
    if(count==np) { *ifail = 0; goto R;}
    (*funct1)(np,b,&p);
    if(mon==1)
    { printf("curr val %8.4f at ", p);
      for(i=1;i<=np;i++) printf(" %8.4f", b[i]);
      printf("\n");
    }
    ifn++;
    if (ifn>=nevals)
    { *ifail = 2;
      printf("qnmin: too many function evaluations\n"); goto R;
    }
  }
  
  /* new lowest value */
  
  *p0 = p; ig++;
  grad(np,b,g,p,funct1,tolder);
  if(mon==1)
  { printf("grad vec:");
    for(i=1;i<=np;i++) printf(" %11.4f", g[i]);
    printf("\n");
  }
  ifn+=np;
  
  /* update hessian */
  
  for(i=1,d1=0.0;i<=np;i++)
  { t[i]*=k; c[i]=g[i]-c[i]; d1+=t[i]*c[i];}
  
  /* check if +ve def addition */
  
  if (d1<=0.0) goto m10;
  for(i=1,d2=0.0;i<=np;i++)
  { for(j=1,s=0.0;j<=np;j++) s+=h[i][j]*c[j];
    x[i] = s;
    d2+=s*c[i];
  }
  d2 = 1.+d2/d1;
  for(i=1;i<=np;i++)
  { for(j=1;j<=np;j++)
        h[i][j]-= (t[i]*x[j]+t[j]*x[i]-d2*t[i]*t[j])/d1;
  }
  goto m40;

#ifdef NOALLOC
 R: return;
#else
 /* free space before returning */
 R: free(x); free(c); free(g); free(t); return;
#endif
}

#ifndef NONUMGR
void grad(int np, double *b, double *g, double p0,
    void (*funct1)(int, double *, double *), double tolder)
/* optional routine to approximate gradient of funct1
  (use analytic derivatives if they are easily computed) 
*/
{ double tol,p1;
  int i;
  tol=tolder; 
  for(i=1;i<=np;i++)
  { b[i]+=tol;
    (*funct1)(np,b,&p1);
    g[i]=(p1-p0)/tol;
    b[i]-=tol;
  }
}
#endif
