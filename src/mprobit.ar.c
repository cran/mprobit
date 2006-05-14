/*
 *  Numerical MLE and SE estimation for Multivariate Binary AR(1) Probit model. 
 *  Varying cluster size; 
 *  Assumption: cluster id and binary response are located at the 1st and last column 
 *  respectively. Cluster id's are positive integers.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "mvn.h"

#define MAXP 24

double **x;     
int nn,nc,dmax,ncl,*y,*dvec,*dstart;        

/* x    -- data matrix
   y    -- response variable
   nn   -- total number of individuals 
   nc   -- number of covariate
   dmax  --  max cluster size
   ncl  -- number of unique clusters 
   dvec -- a cluster size vector
   dstart -- the starting point of each unique cluster
*/

#ifdef MAIN3
main(int argc, char *argv[])
{ extern double **x;
  extern int nn,nc,damx,ncl,*y,*dvec,*dstart;
  int i,j,ieof;
  int nevals, ifail, iprint, np; 
  int nrec, npred, nclbd, *ydat,*id0;
  double *xdat, *h, th[MAXP], nllkval, tolder,tlap;
  char line[200];
  void mprobitar(double *,int *, int *, int *, int *,int *, int *,
                   int *, double *, double *, double *, double *);

  FILE *in;
  // nrec = #records in file  or total number of individuals
  // tlap = time elapsed 
  // first column is clusterid; assume binary response is at the last column
  // other columns are covariates (floating points)
  // output buffer is zero, important if computations not instant

  tolder=1.e-3;
  if(argc==1) 
  { printf("Usage: %s datafile nrec nclbd npred th[1] ... th[np]\n", argv[0]); 
    printf("   where np=#parameters= npred+2, npred=#predictors\n");
    printf("   Datafile has one header line\n");
    exit(0);
  }
  else
  { in=fopen(argv[1],"r");
    nrec=atoi(argv[2]);    // excludes the header line
    nclbd=atoi(argv[3]);   // excludes the header line
    npred=atoi(argv[4]);  np=npred+2; 
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
  mprobitar(xdat,ydat,id0,&nrec,&npred,&nclbd,&nevals,&iprint,&tolder,
              th,&nllkval,h);
  free(xdat); free(ydat); free(id0); free(h);
  fclose(in);
  exit(0);
}
#endif


/* multivariate binary probit AR(1) model  */
void mprobitar(double *xdat,int *ydat, int *id0, int *nrec, int *npred,
       int *nclbd, int *nevals, int *iprint, double *tolder,
       double *th0, double *nllkval, double *h0) 
/* xdat -- covariate data
   ydat -- response variable
   id   -- cluster id
   nrec -- total number of records
   npred -- number of covariates
   nclbd -- max cluster size
   iprint   -- index of printing out intermediate results
   nevals -- number of function evaluations
   tolder -- tolder = stepsize for numerical derivative
   th -- initial values of parameters; at the end it returns the MLE.
   nllkval -- returned value of neg loglik at MLE
   h0 -- returned estimated inverse hessian or cov matrix at MLE
*/
{ extern int nn,nc,dmax;
  extern double **x;
  extern int *y, *dvec, *dstart, ncl;
  int i,j, icl, cnt, target, np, ifail, *id; 
  double **h, *th;
  void nllkar1(int, double *, double *);
  void qnmin(int, double *, double **, double *, int,
  int *, int, void (*funct1)(int, double *, double *), double);
  double **dmatrix(int, int);
  
  // rename pointers to local use
  nn=*nrec; 
  nc=*npred; // number of covariates
  np=nc+2;   // total number of parameters to be estimated

  x=dmatrix(nn+1,nc+1);
  y=(int *)malloc((nn+1) * sizeof(int*));
  id=(int *)malloc((nn+1) * sizeof(int*));
  dvec  =(int *)malloc(((*nclbd)+1) * sizeof(double));
  dstart=(int *)malloc(((*nclbd)+1) * sizeof(double));
  th=(double *)malloc((np+1) * sizeof(double));
  h=dmatrix(np+1,np+1);

  for(i=1;i<=nn;i++) 
  { y[i]=ydat[i-1]; id[i]=id0[i-1];  
    for(j=1;j<=nc;j++)  x[i][j]=xdat[nn*(j-1)+i-1];   
  }
  for(i=1;i<=np;i++)  
  { th[i]=th0[i-1];
    for(j=1;j<=np;j++) h[i][j]=h0[np*(j-1)+i-1];  
  }

  dmax=0;
  icl=0;              
  dstart[icl]=1;      
  target=-99; cnt=0;           // assumes id always > 0
  for(i=1;i<=nn;i++)
  { if(id[i] == target) cnt++; //checks each cluster size 
    else 
    { if(icl>0) 
      { dvec[icl]=cnt; if(cnt>dmax) dmax=cnt; }
      cnt=1; target=id[i]; icl++; 
      if(icl>*nclbd) 
      { printf("Error: insufficient bound, make nclbd larger\n"); return; } 
      dstart[icl]=i;
    } 
  }
  dvec[icl]=cnt;              // last cluster
  if(cnt>dmax) dmax=cnt;      
  ncl=icl;                   

  if(*iprint==1) printf("#clusters=%d before quasi-Newton\n", ncl);

  /* --- quasi-newton minimization of -loglik --- */
  qnmin(np,th,h,nllkval,*nevals,&ifail,*iprint,nllkar1,*tolder);

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
  free(x[0]); free(x); 
  free(y); free(id); free(dvec); free(dstart);
  free(h[0]); free(h); free(th); 
  return;  
} 


/* nllk for mult AR(1) probit, common beta regression coefficients */
void nllkar1(int np, double *th, double *fnv)
{ extern int nn,nc,ncl,dmax,*y,*dvec,*dstart;
  extern double **x;
  double rh[M][M],*pa,*pb;
  double pr,nlk,eps,sd;
  double b0,*b,rho,*rr,tem;
  int i,j,k,jk,jc;
  int nsim,ifail,d,start;
  double pnorms(double);
#ifdef TWO
   void mvn2(int, double [], double [], double [][M], int, double,
       double *, double *, int *);
#else
   void mvn1(int, double [], double [], double [][M], int, double,
       double *, double *, int *);
#endif
  pa=(double *) malloc((dmax+1) * sizeof(double));
  pb=(double *) malloc((dmax+1) * sizeof(double));
  b=(double *) malloc((nc+1) * sizeof(double));
  rr=(double *) malloc((dmax+1) * sizeof(double));
  
  nlk=0.; eps=1.e-6; 
  b0=th[1]; for(j=1;j<=nc;j++) b[j]=th[j+1];
  rho=th[np];
 
  for(jk=1,rr[0]=1.;jk<dmax;jk++) rr[jk]=rr[jk-1]*rho;
  for(j=1;j<=dmax;j++)
  { for(k=1;k<=dmax;k++) rh[j][k]=rr[abs(k-j)]; }

  for(i=1;i<=ncl;i++)
  { d=dvec[i];
    start=dstart[i];
    if(d<=7) nsim=0; else nsim=2000; // bound of the number of permutations 
    for(j=1;j<=d;j++)
    { for(jc=1,tem=b0;jc<=nc;jc++) tem+=b[jc]*x[start+(j-1)][jc]; 
      if(y[start+j-1]==1) { pa[j]=-6.; pb[j]=tem; }
      else { pa[j]=tem; pb[j]=6.; }
    }
#ifdef TWO
    if(d>=2)
    { eps=1.e-3; 
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

