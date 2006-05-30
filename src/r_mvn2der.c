#include <stdio.h>
#include <math.h>
#include "mprobit.h"
/* derivatives of mvn rectangle probability */
/* version with pointers for link to R, zero indexes are used */
/* gcc -DMAIN2 -o r_mvn2der r_mvn2der.c pnorms.c mvnsub.c pmnorm.c pbnorm.c phi.c -lm */
/* second order approx */
#ifdef MAIN2
main()
{ int m,i,j,nsim,ifail,m2,k,ksign;
  //double r[M][M],a[M],b[M],pr,sd,eps;
  double *rhmat,*a,*b,pr,sd,eps;
  /*
  void mvn2(int, double [], double [], double [][M], int, double,
       double *, double *, int *);
  */
  void r_mvndu2(int *, double *, double *, double *, int *,
       int *, double *, int *, double *);
  void r_mvndrh2(int *, double *, double *, double *, int *, int *,
       int *, double *, int *, double *);
  double heps,pr2,dera,derb,drh;
  
  eps=1.e-6; // heps=1.e-4;
  nsim=0;
  scanf("%d", &m);
  while(m>0)
  { a=(double *) malloc((m+1) * sizeof(double));
    b=(double *) malloc((m+1) * sizeof(double));
    m2=m*m;
    rhmat=(double *) malloc((m2+1) * sizeof(double));

    //for(i=1;i<=m;i++)
    //{ for(j=1;j<=m;j++) scanf("%lf", &r[i][j]);}
    printf("\nm=%3d\n", m);
    for(k=0;k<m2;k++) scanf("%lf", &rhmat[k]);
    for(k=0;k<m2;k++) 
    { printf("%f ", rhmat[k]); if(k%m==(m-1)) printf("\n"); }
    for(i=0;i<m;i++)  scanf("%lf", &a[i]);
    for(i=0;i<m;i++)  scanf("%lf", &b[i]);
    for(i=0;i<m;i++) printf("%8.4f", a[i]);  printf("\n");
    for(i=0;i<m;i++) printf("%8.4f", b[i]);  printf("\n");

    printf("derivative wrt a_k,b_k, k=1,...,%d\n", m);
    for(k=1;k<=m;k++)
    { printf("  k=%d\n", k);
      ksign=-k;
      r_mvndu2(&m,a,b,rhmat,&ksign,&nsim,&eps,&ifail,&dera);
      ksign=k;
      r_mvndu2(&m,a,b,rhmat,&ksign,&nsim,&eps,&ifail,&derb);
      printf("   integ: dera=%f, derb=%f\n", dera,derb);
    }

    printf("derivative wrt rho(j,k)\n"); 
    for(j=1;j<m;j++)
    { for(k=j+1;k<=m;k++)
      { r_mvndrh2(&m,a,b,rhmat,&j,&k,&nsim,&eps,&ifail,&drh);
        printf("  (j,k)=%d %d:  derrh=%f\n", j,k,drh);
      }
    }

    free(a); free(b); free(rhmat);
    scanf("%d", &m);
  }
}
#endif

// P(Z_j\in (a_j,b_j);R): deriv wrt a_k or b_k, R is corr matrix
// mean vector is zero
// ksign=k for deriv wrt b_k, =-k for deriv wrt a_k
//double mvndu(int m, double w[], double x[], double r[][M], int ksign,
void r_mvndu2(int *m0, double *w, double *x, double *rhmat, int *ksign,
 int *nsim, double *eps, int *ifail, double *der0)
{ double der,tem,bb,pr,sd;
  double rr[M][M],mu[M],a[M],b[M],rc[M][M],dd[M];
  int i,j,m1,k,m,ii;
  void mvn2(int, double [], double [], double [][M], int, double,
       double *, double *, int *);
  double phi(double);

  m= *m0;
  k=abs(*ksign);
  for(i=1,ii=0;i<=m;i++)
  { for(j=1;j<=m;j++) 
    { rr[i][j]=rhmat[ii]; ii++; }
    a[i]=w[i-1]; b[i]=x[i-1];
  }
  /*for(i=1;i<=m;i++)
  { for(j=1;j<=m;j++) printf("%f ", rr[i][j]);
    printf("\n");
  }*/
  if(k<m) // switch indices m<->k, first flip rows, then columns
  { for(i=1;i<=m;i++)
    { tem=rr[i][m]; rr[i][m]=rr[i][k]; rr[i][k]=tem; }
    for(i=1;i<=m;i++)
    { tem=rr[m][i]; rr[m][i]=rr[k][i]; rr[k][i]=tem; }
    tem=a[m]; a[m]=a[k]; a[k]=tem;
    tem=b[m]; b[m]=b[k]; b[k]=tem;
  }
  /*for(i=1;i<=m;i++)
  { for(j=1;j<=m;j++) printf("%f ", rr[i][j]);  printf("\n"); }
  for(i=1;i<=m;i++) printf("%f ", b[i]); printf("\n");
  */

  // deriv wrt b[m]
  // S_{11.2}= S_{11}-S_{12}S_{21}/s_{22}
  // mu_{11.2}= mu1+S_{12}(b-mu2)/s22
  m1=m-1; 
  if((*ksign)>0) bb=b[m]; else bb=a[m];
  for(i=1;i<=m1;i++)
  { for(j=1;j<=m1;j++) rr[i][j]-=rr[i][m]*rr[m][j]; 
    mu[i]=rr[i][m]*bb;
  }

  /* for(i=1;i<=m1;i++)
  { for(j=1;j<=m1;j++) printf("%f ", rr[i][j]);  printf("\n"); }
  for(i=1;i<=m1;i++) printf("%f ", mu[i]); printf("\n");
  */

  // need to rescale to new correlation matrix
  for(i=1;i<=m1;i++)
  { dd[i]=sqrt(rr[i][i]);
    a[i]=(a[i]-mu[i])/dd[i];
    b[i]=(b[i]-mu[i])/dd[i];
  }
  for(i=1;i<=m1;i++)
  { for(j=1;j<=m1;j++) rc[i][j]=rr[i][j]/(dd[i]*dd[j]); }

  /*for(i=1;i<=m1;i++)
  { for(j=1;j<=m1;j++) printf("%f ", rc[i][j]);  printf("\n"); }
  for(i=1;i<=m1;i++) printf("%f ", b[i]); printf("\n");
  */

  mvn2(m1,a,b,rc,*nsim,*eps,&pr,&sd,ifail);

  der=phi(bb)*pr;
  if(*ksign<0) der=-der;
  *der0=der;
}


// P(Z_j\in (a_j,b_j);R): deriv wrt rh_{jk}, R is corr matrix
// mean vector is zero
//double mvndrh(int m, double w[], double x[], double r[][M], int j1, int k1,
void r_mvndrh2(int *m0, double *w, double *x, double *rhmat, int *j0, int *k0,
 int *nsim, double *eps, int *ifail, double *drh)
{ double der,tem,pr,sd,rh,r1,rjm,rjm1;
  double rr[M][M],a[M],b[M],rc[M][M],dd[M];
  double mu11[M],mu12[M],mu21[M],mu22[M],ac[M],bc[M];
  int i,j,m1,k,temi,m2,ii,m,j1,k1;
  void mvn2(int, double [], double [], double [][M], int, double,
       double *, double *, int *);
  double phi2(double,double,double);

  m= *m0; j1= *j0; k1= *k0;
  if(j1>k1) { temi=j1, j1=k1; k1=temi; }
  for(i=1,ii=0;i<=m;i++)
  { for(j=1;j<=m;j++) 
    { rr[i][j]=rhmat[ii]; ii++; }
    a[i]=w[i-1]; b[i]=x[i-1];
  }
  m1=m-1; m2=m-2;

  // handling special case of m=2
  if(m==2)
  { rh=rr[1][2];
    der=phi2(a[1],a[2],rh);
    der-=phi2(a[1],b[2],rh);
    der-=phi2(b[1],a[2],rh);
    der+=phi2(b[1],b[2],rh);
    *drh=der;
  }

  // if k1=m, j1=m-1, no switch
  // if k1=m, j1<m-1, switch j1 and m-1
  // if k1=m-1, j1<m-1, switch j1 and m
  // if k1,j1<m-1, switch j1 and m-1, k1 and m
  if(k1==m)
  { if(j1<m1)
    { for(i=1;i<=m;i++)
      { tem=rr[i][m1]; rr[i][m1]=rr[i][j1]; rr[i][j1]=tem; }
      for(i=1;i<=m;i++)
      { tem=rr[m1][i]; rr[m1][i]=rr[j1][i]; rr[j1][i]=tem; }
      tem=a[m1]; a[m1]=a[j1]; a[j1]=tem;
      tem=b[m1]; b[m1]=b[j1]; b[j1]=tem;
    }
  }
  else if(k1==m1)
  { for(i=1;i<=m;i++)
    { tem=rr[i][m]; rr[i][m]=rr[i][j1]; rr[i][j1]=tem; }
    for(i=1;i<=m;i++)
    { tem=rr[m][i]; rr[m][i]=rr[j1][i]; rr[j1][i]=tem; }
    tem=a[m]; a[m]=a[j1]; a[j1]=tem;
    tem=b[m]; b[m]=b[j1]; b[j1]=tem;
  }
  else /* k1<m-1 */
  { for(i=1;i<=m;i++)
    { tem=rr[i][m1]; rr[i][m1]=rr[i][j1]; rr[i][j1]=tem; }
    for(i=1;i<=m;i++)
    { tem=rr[m1][i]; rr[m1][i]=rr[j1][i]; rr[j1][i]=tem; }
    tem=a[m1]; a[m1]=a[j1]; a[j1]=tem;
    tem=b[m1]; b[m1]=b[j1]; b[j1]=tem;
    for(i=1;i<=m;i++)
    { tem=rr[i][m]; rr[i][m]=rr[i][k1]; rr[i][k1]=tem; }
    for(i=1;i<=m;i++)
    { tem=rr[m][i]; rr[m][i]=rr[k1][i]; rr[k1][i]=tem; }
    tem=a[m]; a[m]=a[k1]; a[k1]=tem;
    tem=b[m]; b[m]=b[k1]; b[k1]=tem;
  }
  /*for(j=1;j<=m;j++)
  { for(k=1;k<=m;k++) printf("%f ", rr[j][k]);  printf("\n"); }
  for(j=1;j<=m;j++) printf("%f ", b[j]); printf("\n");
  */

  // rescale corr matrix for first m-2 
  rh=rr[m1][m]; r1=1./(1.-rh*rh);
  for(j=1;j<=m2;j++)
  { rjm=rr[j][m]; rjm1=rr[j][m1];
    for(k=1;k<=m2;k++) 
    { rr[j][k]-=r1*(rjm1*rr[m1][k]+rjm*rr[m][k]
                   -rh*rjm1*rr[m][k]-rh*rjm*rr[m1][k]); 
    }
    // four sets of mus
    mu11[j]=r1*(rjm1*a[m1]+rjm*a[m]-rh*rjm1*a[m]-rh*rjm*a[m1]);
    mu12[j]=r1*(rjm1*a[m1]+rjm*b[m]-rh*rjm1*b[m]-rh*rjm*a[m1]);
    mu21[j]=r1*(rjm1*b[m1]+rjm*a[m]-rh*rjm1*a[m]-rh*rjm*b[m1]);
    mu22[j]=r1*(rjm1*b[m1]+rjm*b[m]-rh*rjm1*b[m]-rh*rjm*b[m1]);
  }
  /*for(j=1;j<=m2;j++)
  { printf("%f %f %f %f\n", mu11[j],mu12[j],mu21[j],mu22[j]); }
  */

  /*for(j=1;j<=m2;j++)
  { for(k=1;k<=m2;k++) printf("%f ", rr[j][k]);  printf("\n"); }
  for(j=1;j<=m2;j++) printf("%f ", mu22[j]); printf("\n");
  */

  // rescaling mus and rr, 4 mvn calcs
  for(j=1;j<=m2;j++)
  { dd[j]=sqrt(rr[j][j]);
    ac[j]=(a[j]-mu11[j])/dd[j];
    bc[j]=(b[j]-mu11[j])/dd[j];
  }
  for(j=1;j<=m2;j++)
  { for(k=1;k<=m2;k++) rc[j][k]=rr[j][k]/(dd[j]*dd[k]); }

  /*for(j=1;j<=m2;j++)
  { for(k=1;k<=m2;k++) printf("%f ", rc[j][k]);  printf("\n"); }
  for(j=1;j<=m2;j++) printf("%f ", bc[j]); printf("\n");
  */

  der=0.;
  // simplify this by first checking if a<-6 or b>6
  if(a[m1]>-6. && a[m]>-6.)
  { mvn2(m2,ac,bc,rc,*nsim,*eps,&pr,&sd,ifail);
    der+=phi2(a[m1],a[m],rh)*pr;
  }
  //printf("%f ", der);
  // other terms
  for(j=1;j<=m2;j++)
  { ac[j]=(a[j]-mu12[j])/dd[j]; bc[j]=(b[j]-mu12[j])/dd[j]; }
  if(a[m1]>-6. && b[m]<6.)
  { mvn2(m2,ac,bc,rc,*nsim,*eps,&pr,&sd,ifail);
    der-=phi2(a[m1],b[m],rh)*pr;
  }
  //printf("%f ", der);
  for(j=1;j<=m2;j++)
  { ac[j]=(a[j]-mu21[j])/dd[j]; bc[j]=(b[j]-mu21[j])/dd[j]; }
  if(a[m]>-6. && b[m1]<6.)
  { mvn2(m2,ac,bc,rc,*nsim,*eps,&pr,&sd,ifail);
    der-=phi2(b[m1],a[m],rh)*pr;
  }
  //printf("%f ", der);
  for(j=1;j<=m2;j++)
  { ac[j]=(a[j]-mu22[j])/dd[j]; bc[j]=(b[j]-mu22[j])/dd[j]; }
  if(b[m]<6. && b[m1]<6.)
  { mvn2(m2,ac,bc,rc,*nsim,*eps,&pr,&sd,ifail);
    der+=phi2(b[m1],b[m],rh)*pr;
  }
  //printf("%f\n", der);
  *drh=der;
}

#ifdef ALONE
double phi(double z)
{ return( 0.398942280401*exp(-.5*z*z));}

double phi2(double z1, double z2, double rh)
{ double r1,tem;
  r1=1.-rh*rh;
  tem=(z1*z1+z2*z2-2.*rh*z1*z2)/r1;
  return( 0.1591549430918953*exp(-.5*tem)/sqrt(r1));
}
#endif
