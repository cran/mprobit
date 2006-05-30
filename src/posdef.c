#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "mprobit.h"

/* check for positive definiteness of a symmetric matrix */
#ifdef MAIN2
main(int argc, char *argv[])
{ double rh[M][M], a[M][M];
  int i, j, d;
  int pdef(double[][M] , double[][M] , int, double);
  int ipdef;

  // matrix index starts at 1 
  scanf("%d", &d);
  while(d>0)
  { for(i=1;i<=d;i++) 
    { for(j=1;j<=d;j++) scanf("%lf", &rh[i][j]); }
    for(i=1;i<=d;i++) 
    { for(j=1;j<=d;j++) printf(" %8.4f", rh[i][j]); 
      printf("\n");
    }
    ipdef=pdef(rh,a,d,1.e-10);
    if(ipdef==1) printf("matrix is positive-definite \n");
    else printf("matrix is not positive-definite \n");
    scanf("%d", &d);
  }
  exit(0);
}
#endif


// Using chol() to check positive-definite of input correlation matrix //
// return 0, if not positive-definite
// a[][]=input n*n matrix, assumed symmetric 
// l[][]=output lower-triangular Cholesky matrix
// toler=tolerance value in one comparison
int pdef(double a[][M], double l[][M], int n, double toler)
{ int i,j,k;
  double sum;

  if(a[1][1]<=0) { return 0; }
  else l[1][1]=sqrt(a[1][1]);
  for(k=2;k<=n;k++)
  { for(i=1;i<k;i++)
    { sum=0.;
      for(j=1;j<i;j++) sum+=l[i][j]*l[k][j];
      if(fabs(a[k][i]-sum)>toler) 
      { l[k][i]=(a[k][i]-sum)/l[i][i]; }
      else l[k][i]=0.;
      l[i][k]=0.;
    }
    for(j=1,sum=0.;j<k;j++)  sum+=l[k][j]*l[k][j];
    if(a[k][k]-sum<=0.) { l[k][k]=0.; return 0; }
    else l[k][k]=sqrt(a[k][k]-sum);
  }
#ifdef MAIN2
  printf("Cholesky matrix\n");
  for(i=1;i<=n;i++)
  { for(j=1;j<=n;j++) printf(" %8.4f ",l[i][j]);
    printf("\n");
  }
#endif
  return 1;
}

