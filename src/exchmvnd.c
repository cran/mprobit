#include <stdio.h>
#include <math.h>
//#define M 20
//#define UB 7.
#define UB 6.
#define EPS 1.e-7
/* gcc -DMAIN2 -o exchmvnd exchmvnd.c exchmvn.c pnorms.c phi.c romberg.c -lm */
/* mvn rectangle probability and derivatives for positive exch case, 
   with Romberg integration */
int mm,kk,ksign;
//double ww[M],xx[M],rs,r1,r32;
double *ww,*xx,rs,r1,r32;
#ifdef MAIN2
main()
{ int m,i;
  //double rh,a[M],b[M],pr,eps;
  double rh,*a,*b,pr,eps;
  //double exchmvn(int, double [], double [], double, double);
  //double emvnd(int, double [], double [], double, int, int, double);
  //double emvndrh(int, double [], double [], double, double);
  double exchmvn(int, double *, double *, double, double);
  double emvnd(int, double *, double *, double, int, int, double);
  double emvndrh(int, double *, double *, double, double);
  double heps,pr2,dera,derb,drh;
  
  eps=1.e-6;
  heps=1.e-4;
  scanf("%d", &m);
  while(m>0)
  { scanf("%lf", &rh);
    a=(double *) malloc((m+1) * sizeof(double));
    b=(double *) malloc((m+1) * sizeof(double));
    for(i=1;i<=m;i++)  scanf("%lf", &a[i]);
    for(i=1;i<=m;i++)  scanf("%lf", &b[i]);
    printf("m=%3d, rh=%6.2f\n", m,rh);
    for(i=1;i<=m;i++) printf("%8.4f", a[i]);  printf("\n");
    for(i=1;i<=m;i++) printf("%8.4f", b[i]);  printf("\n");
    pr=exchmvn(m,a,b,rh,eps);
    printf("exch.  : %9.5f\n", pr);
    a[1]+=heps;
    pr2=exchmvn(m,a,b,rh,eps);
    dera=(pr2-pr)/heps;
    a[1]-=heps;
    b[1]+=heps;
    pr2=exchmvn(m,a,b,rh,eps);
    derb=(pr2-pr)/heps;
    b[1]-=heps;
    printf("num deriv: dera1=%f, derb1=%f\n", dera,derb);
    dera=emvnd(m,a,b,rh,1,-1,eps);
    derb=emvnd(m,a,b,rh,1,1,eps);
    printf("integ:     dera1=%f, derb1=%f\n", dera,derb);
    
    a[m]+=heps;
    pr2=exchmvn(m,a,b,rh,eps);
    dera=(pr2-pr)/heps;
    a[m]-=heps;
    b[m]+=heps;
    pr2=exchmvn(m,a,b,rh,eps);
    derb=(pr2-pr)/heps;
    b[m]-=heps;
    printf("num deriv: dera1=%f, derb1=%f\n", dera,derb);
    dera=emvnd(m,a,b,rh,m,-1,eps);
    derb=emvnd(m,a,b,rh,m,1,eps);
    printf("integ:     dera1=%f, derb1=%f\n", dera,derb);

    rh+=heps;
    pr2=exchmvn(m,a,b,rh,eps);
    drh=(pr2-pr)/heps;
    rh-=heps;
    printf("num deriv: derrh=%f\n", drh);
    drh=emvndrh(m,a,b,rh,eps);
    printf("integ.   : derrh=%f\n", drh);
    free(a); free(b);
    scanf("%d", &m);
  }
}
#endif

/* P(Z_j\in (a_j,b_j)): deriv wrt a_k or b_k,
   ks=-1 for a_k, ks=1 for b_k
*/
//double emvnd(int m, double w[], double x[], double rh, int k, int ks, 
double emvnd(int m, double *w, double *x, double rh, int k, int ks, 
  double eps)
{ double gd(double),der;
  double romberg(double (*)(double), double, double, double);
  int i;
  extern int mm,kk,ksign;
  //extern double ww[M],xx[M],rs,r1;
  extern double *ww,*xx,rs,r1;
  mm=m; kk=k; rs=sqrt(rh); r1=sqrt(1.-rh); ksign=ks;
  xx=(double *) malloc((mm+1) * sizeof(double));
  ww=(double *) malloc((mm+1) * sizeof(double));
  for(i=1;i<=mm;i++) { ww[i]=w[i]; xx[i]=x[i]; }
  der=romberg(gd,-UB,UB,eps);
  free(xx); free(ww);
  return(ks*der);
}

/* integrand for emvnd */
double gd(double z)
{ double pnorms(double),phi(double),a,b;
  extern int mm,kk;
  //extern double ww[M],xx[M],rs,r1;
  extern double *ww,*xx,rs,r1;
  int i;
  double tem;
  for(i=1,tem=1.;i<=mm;i++)
  { if(i!=kk)
    { a=(ww[i]-rs*z)/r1; b=(xx[i]-rs*z)/r1;
      tem*=pnorms(b)-pnorms(a);
    }
    else if(ksign==-1)
    { a=(ww[i]-rs*z)/r1; tem*=phi(a)/r1; }
    else
    { b=(xx[i]-rs*z)/r1; tem*=phi(b)/r1; }
  }
  tem*=phi(z);
  return(tem);
}


/* P(Z_j\in (a_j,b_j)): deriv wrt rho */
//double emvndrh(int m, double w[], double x[], double rh, double eps)
double emvndrh(int m, double *w, double *x, double rh, double eps)
{ double grh(double),der,*t,tem,sum;  // t[M]
  double romberg(double (*)(double), double, double, double);
  double pnorms(double),phi(double);
  int i,k;
  extern int mm;
  //extern double ww[M],xx[M],rs,r1,r32;
  extern double *ww,*xx,rs,r1,r32;
  mm=m; rs=sqrt(rh); r1=sqrt(1.-rh); r32=r1*(1.-rh);
  xx=(double *) malloc((mm+1) * sizeof(double));
  ww=(double *) malloc((mm+1) * sizeof(double));
  t=(double *) malloc((mm+1) * sizeof(double));
  for(i=1;i<=mm;i++) { ww[i]=w[i]; xx[i]=x[i]; }
  if(rh>=0.) der=romberg(grh,-UB,UB,eps);
  else /* rho=0 */
  { for(i=1;i<=mm;i++) t[i]=pnorms(x[i])-pnorms(w[i]);
    for(k=1,sum=0.;k<=mm;k++)
    { for(i=1,tem=1.;i<=mm;i++)
      { if(i!=k) { tem*=t[i]; }
        else tem*=(x[i]*phi(x[i])-w[i]*phi(w[i]));
        // maybe check if x>10 or w<-10?
      }
      sum+=tem;
    }
    der=.5*sum;
  }
  free(xx); free(ww); free(t);
  return(der);
}

/* integrand for emvndrh */
double grh(double z)
{ double pnorms(double),phi(double),a,b;
  extern int mm;
  //extern double ww[M],xx[M],rs,r1,r32;
  extern double *ww,*xx,rs,r1,r32;
  int i,k;
  //double tem,sum,tem2,t[M];
  double tem,sum,tem2,*t;
  // how to handle the limit as rho->0?
  t=(double *) malloc((mm+1) * sizeof(double));
  for(i=1;i<=mm;i++) 
  { a=(ww[i]-rs*z)/r1; b=(xx[i]-rs*z)/r1;
    t[i]=pnorms(b)-pnorms(a);
  }
  for(k=1,sum=0.;k<=mm;k++)
  { for(i=1,tem=1.;i<=mm;i++)
    { if(i!=k) { tem*=t[i]; }
      else
      { a=(ww[i]-rs*z)/r1; b=(xx[i]-rs*z)/r1;
        tem2=phi(b)*(xx[i]-z/rs)-phi(a)*(ww[i]-z/rs);
        tem2*=.5/r32;
        tem*=tem2;
      }
    }
    sum+=tem;
  }
  tem=sum*phi(z);
  free(t);
  return(tem);
}

