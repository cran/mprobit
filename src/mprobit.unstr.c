/* 
 *  MLE and standard error estimation for multivariate probit model. 
 *  Constant cluster size, unstructured correlation matrix.
 *  Assumptions: 
 *    cluster id's and binary responses are at the 1st, and last column
 *  respectively. Other columns are covariates.  The cluster id's are positive
 *  integers and in consecutive order.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mprobit.h"

#ifdef MAIN3
main(int argc, char *argv[])
{ extern double **x;
  extern int *y, nn, nc, ncl, d;
  int i,j,ieof, nevals, ifail, iprint, np; 
  int nrec, npred, *ydat, *id0;
  double *xdat, *h, th[MAXP], nllkval, tolder;
  char line[200];
  void mprobitunstr(double *,int *, int *, int *, int *, int *, int *, 
              int *, double *, double *, double *, double *);
  FILE *in;

  tolder=1.e-3;
  if(argc==1) 
  { printf("Usage: %s datafile nrec npred d th[1] ... th[np]\n", argv[0]);
    printf("   where np=#parameters= npred+1+d(d-1)/2, npred=#predictors\n");
    printf("   and   d=#repeated measures per cluster\n");
    printf("   Datafile has one header line\n");
    exit(0);
  }
  else
  { in=fopen(argv[1],"r");
    nrec=atoi(argv[2]);    // excludes the header
    npred=atoi(argv[3]); 
    d=atoi(argv[4]); 
    np=npred+1+d*(d-1)/2; 
    for(i=0;i<np;i++)  th[i]=atof(argv[5+i]); 
  }  
  setbuf(stdout,NULL); 
  // allocate memories
  xdat=(double *)malloc((nrec*npred+1) * sizeof(double));
  ydat=(int *)malloc((nrec+1) * sizeof(int));
  id0=(int *)malloc((nrec+1) * sizeof(int));
  h=(double *)malloc((np+1)*(np+1) * sizeof(double)); 

  // read in data set with 1 header line
  fgets(line,200,in);     
  for(i=1;i<=nrec;i++)
  { ieof=fscanf(in,"%d", &id0[i-1]); 
    if(ieof==EOF) break;
    for(j=1;j<=npred;j++)  fscanf(in,"%lf", &xdat[nrec*(j-1)+i-1]); 
    fscanf(in,"%d", &ydat[i-1]); 
  }

  /*--- quasi-newton minimization of -loglik ---*/
  nevals=1000;  iprint=1; 
  mprobitunstr(xdat,ydat,id0,&nrec,&npred,&d,&nevals,&iprint,
               &tolder,th,&nllkval,h);
  free(xdat); free(ydat); free(id0); free(h);
  fclose(in);
  exit(0);
}
#endif


void mprobitunstr(double *xdat,int *ydat, int *id0, int *nrec, int *npred,
                  int *d0, int *nevals, int *iprint, double *tolder,
                  double *th0, double *nllkval, double *h0) 
/* xdat -- covariate data
   ydat -- response variable
   id   -- cluster id
   nrec -- total number of records
   npred -- number of covariates
   iprint   -- index of printing out intermediate results
   nevals -- number of function evaluations
   tolder -- tolder = stepsize for numerical derivative
   th0 -- initial values of parameters; at the end it returns the MLE.
   nllkval -- returned value of neg loglik at MLE
   h0 -- returned estimated inverse hessian or cov matrix at MLE
*/
{ extern double **x;
  extern int *y, nn, nc, ncl, d;
  int i,j, np, ifail, *id; 
  double **h, *th;
  void nllkunstr(int, double *, double *);
  void qnmin(int, double *, double **, double *, int,
  int *, int, void (*funct1)(int, double *, double *), double);
  double **dmatrix(int, int);
  
  // rename pointers to local use
  nn=*nrec; 
  nc=*npred; // number of covariates
  d=*d0;
  ncl=nn/d;
  np=nc+1+d*(d-1)/2;    // total number of parameters to be estimated

  x=dmatrix(nn+1,nc+1);
  y=(int *)malloc((nn+1) * sizeof(int*));
  id=(int *)malloc((nn+1) * sizeof(int*));
  th=(double *)malloc((np+1) * sizeof(double));
  h=dmatrix(np+1,np+1);

  for(i=1;i<=nn;i++) 
  { y[i]=ydat[i-1]; id[i]=id0[i-1];  
    for(j=1;j<=nc;j++)  
    { x[i][j]=xdat[nn*(j-1)+i-1];   
    }
  }
  for(i=1;i<=np;i++)  
  { th[i]=th0[i-1];
    for(j=1;j<=np;j++) h[i][j]=h0[np*(j-1)+i-1];  
  }

  if(*iprint==1) printf("#clusters=%d\n",ncl);

  /* --- quasi-newton minimization of -loglik --- */
  qnmin(np,th,h,nllkval,*nevals,&ifail,*iprint,nllkunstr,*tolder);
  if(*iprint==1)  
  { printf("\nNeg log-likelihood = %8.4f, end code = %3d\n",*nllkval,ifail);
    printf("\nParameter estimates and standard errors\n");
    for(i=1;i<=np;i++)
     printf("  %8.4f %8.4f\n", th[i],sqrt(h[i][i]));  
    // inverse Hessian matrix
    printf("\nEstimated asymptotic covariance matrix\n");
    for(i=1;i<=np;i++)
    { for(j=1;j<=np;j++) printf(" %10.6f", h[i][j]); 
      printf("\n");
    }
  }

  for(i=1;i<=np;i++)
  { th0[i-1]=th[i];
    for(j=1;j<=np;j++) h0[np*(j-1)+i-1]=h[i][j]; 
  }
  free(x[0]); free(x); free(y); free(id);
  free(h[0]); free(h); free(th); 
  return;  
} 


/* nllk for mult-variate binary probit with unstructured correlation matrix;
   common beta regression coefficients
*/
void nllkunstr(int np, double *th, double *fnv)
{ extern int *y,nn,nc,ncl,d;
  extern double **x;
  double rh[M][M],L[M][M],*pa,*pb;
  double pr,nlk,eps,sd;
  double b0,*b,*rr,tem;
  int i,j,k,jc;
  int nsim,ifail;
  double pnorms(double);
  int pdef(double[][M] , double[][M] , int, double);

#ifdef TWO
   void mvn2(int, double [], double [], double [][M], int, double,
       double *, double *, int *);
#else
   void mvn1(int, double [], double [], double [][M], int, double,
       double *, double *, int *);
#endif
  pa=(double *) malloc((d+1) * sizeof(double));
  pb=(double *) malloc((d+1) * sizeof(double));
  b=(double *) malloc((nc+1) * sizeof(double));
  rr=(double *) malloc((d+1) * sizeof(double));
  
  nlk=0.; eps=1.e-6; 
  b0=th[1]; for(j=1;j<=nc;j++) b[j]=th[j+1];
  for(j=nc+2;j<=np;j++)
  { if(th[j]>=1. || th[j]<=-1.) 
    { *fnv=1.e9; free(b); free(pa); free(pb); free(rr); return; }
  }

  for(i=1,k=1;i<=d;i++)
  { for(j=i;j<=d;j++)
    { if(i==j)  rh[i][j]=1;
      else
      { rh[i][j]=th[nc+1+k]; k++; 
        if(rh[i][j]>1 || rh[i][j]<-1) 
        {  *fnv=1.e9; free(b); free(pa); free(pb); free(rr); return; }
      }
    }
   for(j=1;j<=i-1;j++)  rh[i][j]=rh[j][i]; 
  }

  if(!pdef(rh,L,d,1.e-10)) 
  { *fnv=1.e9; free(b); free(pa); free(pb); free(rr); return;
  } 
  for(i=1;i<=ncl;i++)
  { if(d<=7) nsim=0; else nsim=2000; 
    for(j=1;j<=d;j++)
    { for(jc=1,tem=b0;jc<=nc;jc++) tem+=b[jc]*x[d*(i-1)+j][jc];
      if(y[d*(i-1)+j]==1) { pa[j]=-6.; pb[j]=tem; }
      else { pa[j]=tem; pb[j]=6.; }
    }
#ifdef TWO
    if(d>=2)
    { eps=1.e-5; 
      mvn2(d,pa,pb,rh,nsim,eps,&pr,&sd,&ifail);
      if(ifail>0) 
      { *fnv=1.e9; free(b); free(pa); free(pb); free(rr); return; }
    }
#else
    if(d>=2)
    { mvn1(d,pa,pb,rh,nsim,eps,&pr,&sd,&ifail);
      if(ifail>0)
      { *fnv=1.e9; free(b); free(pa); free(pb); free(rr); return; }
    }
#endif
    else pr=pnorms(pb[1])-pnorms(pa[1]);
    if(pr<=-.1)
    { *fnv=1.e9;
      for(j=1;j<=d;j++)
          printf("pa[%d]=%6.3f, pb[%d]=%6.3f\n", j,pa[j],j,pb[j]);
      return;
    }
    if(pr<=0.) pr=1.e-15;
    nlk-=log(pr);
  }
  *fnv=nlk;
  free(b); free(pa); free(pb); free(rr);
  return;
}

