#include <stdio.h>
#include <math.h>
#define M 20
//#define UB 7.
#define UB 6.
#define EPS 1.e-7
/* gcc -DMAIN2 -o r_exchmvnd r_exchmvnd.c r_exchmvn.c pnorms.c phi.c romberg.c -lm */
/* mvn rectangle probability and derivatives for positive exch case, 
   with Romberg integration */
/* version with pointers for link to R, zero indexes are used */
int mm,kk,ksign;
double *ww2,*xx2,rs,r1,r32;
#ifdef MAIN2
main()
{ int m,i,k,ks;
  double rh,*a,*b,pr,eps;
  void r_exchmvn(int *, double *, double *, double *, double *, double *);
  void r_emvnd(int *, double *,double *,double *, int *,int *, double *,double *);
  void r_emvndrh(int *, double *, double *, double *, double *, double *);
  double heps,pr2,dera,derb,drh;
  
  eps=1.e-6;
  heps=1.e-4;
  scanf("%d", &m);
  while(m>0)
  { scanf("%lf", &rh);
    a=(double *) malloc(m * sizeof(double));
    b=(double *) malloc(m * sizeof(double));
    for(i=0;i<m;i++)  scanf("%lf", &a[i]);
    for(i=0;i<m;i++)  scanf("%lf", &b[i]);
    printf("m=%3d, rh=%6.2f\n", m,rh);
    for(i=0;i<m;i++) printf("%8.4f", a[i]);  printf("\n");
    for(i=0;i<m;i++) printf("%8.4f", b[i]);  printf("\n");
    r_exchmvn(&m,a,b,&rh,&eps,&pr);
    printf("exch.  : %9.5f\n", pr);
    a[0]+=heps;
    r_exchmvn(&m,a,b,&rh,&eps,&pr2);
    dera=(pr2-pr)/heps;
    a[0]-=heps;
    b[0]+=heps;
    r_exchmvn(&m,a,b,&rh,&eps,&pr2);
    derb=(pr2-pr)/heps;
    b[0]-=heps;
    printf("num deriv: dera1=%f, derb1=%f\n", dera,derb);
    k=1; ks=-1;
    r_emvnd(&m,a,b,&rh,&k,&ks,&eps,&dera);
    k=1; ks=1;
    r_emvnd(&m,a,b,&rh,&k,&ks,&eps,&derb);
    printf("integ:     dera1=%f, derb1=%f\n", dera,derb);
    
    a[m-1]+=heps;
    r_exchmvn(&m,a,b,&rh,&eps,&pr2);
    dera=(pr2-pr)/heps;
    a[m-1]-=heps;
    b[m-1]+=heps;
    r_exchmvn(&m,a,b,&rh,&eps,&pr2);
    derb=(pr2-pr)/heps;
    b[m-1]-=heps;
    printf("num deriv: dera1=%f, derb1=%f\n", dera,derb);
    k=m; ks=-1;
    r_emvnd(&m,a,b,&rh,&k,&ks,&eps,&dera);
    k=m; ks=1;
    r_emvnd(&m,a,b,&rh,&k,&ks,&eps,&derb);
    printf("integ:     dera1=%f, derb1=%f\n", dera,derb);

    rh+=heps;
    r_exchmvn(&m,a,b,&rh,&eps,&pr2);
    drh=(pr2-pr)/heps;
    rh-=heps;
    printf("num deriv: derrh=%f\n", drh);
    r_emvndrh(&m,a,b,&rh,&eps,&drh);
    printf("integ.   : derrh=%f\n", drh);
    free(a); free(b);
    scanf("%d", &m);
  }
}
#endif

/* P(Z_j\in (a_j,b_j)): deriv wrt a_k or b_k,
   ks=-1 for a_k, ks=1 for b_k
*/
void r_emvnd(int *m, double *w, double *x, double *rh, int *k, int *ks, 
  double *eps, double *deriv)
{ double r_gd(double),der;
  double romberg(double (*)(double), double, double, double);
  int i;
  extern int mm,kk,ksign;
  extern double *ww2,*xx2,rs,r1;
  mm=*m; kk=(*k-1); rs=sqrt(*rh); r1=sqrt(1.-(*rh)); ksign=*ks;
  xx2=(double *) malloc(mm * sizeof(double));
  ww2=(double *) malloc(mm * sizeof(double));
  for(i=0;i<mm;i++) { ww2[i]=w[i]; xx2[i]=x[i]; }
  der=romberg(r_gd,-UB,UB,*eps);
  free(xx2); free(ww2);
  *deriv= ksign*der;
}

/* integrand for emvnd */
double r_gd(double z)
{ double pnorms(double),phi(double),a,b;
  extern int mm,kk;
  extern double *ww2,*xx2,rs,r1;
  int i;
  double tem;
  for(i=0,tem=1.;i<mm;i++)
  { if(i!=kk)
    { a=(ww2[i]-rs*z)/r1; b=(xx2[i]-rs*z)/r1;
      tem*=pnorms(b)-pnorms(a);
    }
    else if(ksign==-1)
    { a=(ww2[i]-rs*z)/r1; tem*=phi(a)/r1; }
    else
    { b=(xx2[i]-rs*z)/r1; tem*=phi(b)/r1; }
  }
  tem*=phi(z);
  return(tem);
}


/* P(Z_j\in (a_j,b_j)): deriv wrt rho */
void r_emvndrh(int *m, double *w, double *x, double *rh, double *eps, 
  double *deriv)
{ double r_grh(double),der,tem,sum;
  double romberg(double (*)(double), double, double, double);
  double pnorms(double),phi(double);
  int i,k;
  double *t;
  extern int mm;
  extern double *ww2,*xx2,rs,r1,r32;
  mm=*m; rs=sqrt(*rh); r1=sqrt(1.-(*rh)); r32=r1*(1.-(*rh));
  xx2=(double *) malloc(mm * sizeof(double));
  ww2=(double *) malloc(mm * sizeof(double));
  t=(double *) malloc(mm * sizeof(double));
  for(i=0;i<mm;i++) { ww2[i]=w[i]; xx2[i]=x[i]; }
  if((*rh)>=0.) der=romberg(r_grh,-UB,UB,*eps);
  else /* rho=0 */
  { for(i=0;i<mm;i++) t[i]=pnorms(x[i])-pnorms(w[i]);
    for(k=0,sum=0.;k<mm;k++)
    { for(i=0,tem=1.;i<mm;i++)
      { if(i!=k) { tem*=t[i]; }
        else tem*=(x[i]*phi(x[i])-w[i]*phi(w[i]));
        // maybe check if x>10 or w<-10?
      }
      sum+=tem;
    }
    der=.5*sum;
  }
  free(xx2); free(ww2); free(t);
  *deriv=der;
}

/* integrand for emvndrh */
double r_grh(double z)
{ double pnorms(double),phi(double),a,b;
  extern int mm;
  extern double *ww2,*xx2,rs,r1,r32;
  int i,k;
  double tem,sum,tem2,*t;

  t=(double *) malloc(mm * sizeof(double));
  for(i=0;i<mm;i++) 
  { a=(ww2[i]-rs*z)/r1; b=(xx2[i]-rs*z)/r1;
    t[i]=pnorms(b)-pnorms(a);
  }
  for(k=0,sum=0.;k<mm;k++)
  { for(i=0,tem=1.;i<mm;i++)
    { if(i!=k) { tem*=t[i]; }
      else
      { a=(ww2[i]-rs*z)/r1; b=(xx2[i]-rs*z)/r1;
        tem2=phi(b)*(xx2[i]-z/rs)-phi(a)*(ww2[i]-z/rs);
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

