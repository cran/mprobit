#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mprobit.h"

/* gcc -DMAIN -o r_exchmvn r_exchmvn.c pnorms.c phi.c romberg.c -lm */
/* mvn rectangle probability for positive exch case, */
/* version with pointers for link to R, zero indexes are used */
#ifdef MAIN
main()
{ int m,i;
  double rh,*x,*w,pr,eps;
  void r_exchmvn(int *, double *, double *, double *, double *, double *);
  
  eps=1.e-6;
  scanf("%d", &m);
  while(m>0)
  { scanf("%lf", &rh);
    x=(double *) malloc(m * sizeof(double));
    w=(double *) malloc(m * sizeof(double));
    for(i=0;i<m;i++)  scanf("%lf", &w[i]);
    for(i=0;i<m;i++)  scanf("%lf", &x[i]);
    printf("m=%3d, rh=%f\n", m,rh);
    for(i=0;i<m;i++) printf("%8.4f", w[i]);  printf("\n");
    for(i=0;i<m;i++) printf("%8.4f", x[i]);  printf("\n");
    r_exchmvn(&m,w,x,&rh,&eps, &pr);
    printf("exch.  : %.10f\n", pr);
    free(x); free(w);
    scanf("%d", &m);
  }
}
#endif

/* version with all pointers and zero indices for interface to R directly */
void r_exchmvn(int *m, double *w, double *x, double *rh, double *eps, double *pr)
{ double r_g(double);
  double romberg(double (*)(double), double, double, double);
  int i;
  extern int mm;
  extern double *ww,*xx,rs,r1;
  mm=*m; rs=sqrt(*rh); r1=sqrt(1.-(*rh));
  xx=(double *) malloc(mm * sizeof(double));
  ww=(double *) malloc(mm * sizeof(double));
  for(i=0;i<mm;i++) { ww[i]=w[i]; xx[i]=x[i];}
  *pr=romberg(r_g,-UB,UB,*eps);
  free(xx); free(ww);
}

double r_g(double z)
{ double pnorms(double),phi(double),a,b;
  extern int mm;
  extern double *ww,*xx,rs,r1;
  int i;
  double tem;
  for(i=0,tem=1.;i<mm;i++)
  { a=(ww[i]-rs*z)/r1; b=(xx[i]-rs*z)/r1;
    tem*=pnorms(b)-pnorms(a);
  }
  tem*=phi(z);
  return(tem);
}

/*
double phi(double z)
{ return( 0.3989422804014327*exp(-.5*z*z));}
*/
#undef UB

