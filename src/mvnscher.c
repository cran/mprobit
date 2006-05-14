#include <stdio.h>
#include <math.h>
/*  program for link from R/Splus to mulnor */
void mvnscher(double ub[], double lb[], double sig0[], double *eps0, int *n0, 
   int inf0[], double *prob, double *bound, int *ifault)
{  double eps;
   int n;
 /*int i,j; */
   void mulnor(double [], double [], double [], double, int,
               int [], double *, double *, int *);
   eps= *eps0; n= *n0;

 /* for(i=0;i<n;i++) printf("%8.4f", lb[i]);  printf("\n");
    for(i=0;i<n;i++) printf("%8.4f", ub[i]);  printf("\n");
    for(j=0;j<(n*(n-1))/2;j++) printf("%8.4f", sig0[j]); printf("\n"); 
 */
   mulnor(ub,lb,sig0,eps,n,inf0,prob,bound,ifault);
 /* printf("%f %f %d\n", *prob,*bound,*ifault); */
}
