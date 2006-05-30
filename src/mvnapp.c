#include <stdio.h>
#include <math.h>
#include "mprobit.h"
/* first order approx., cond. exp. */
/* second order approx., cond. exp. */
void mvnapp(int *type, int *m0, double w0[], double x0[], double corr[],
    int *nsim, double *eps,
    double *pr, double *sd, int *ifail)
{  double r[M][M];
   double w[M],x[M];
   double p,perr;
   int i,m,j,k,ierr;
   void mvn1(int, double [], double [], double [][M], int, double,
       double *, double *, int *);
   void mvn2(int, double [], double [], double [][M], int, double, 
       double *, double *, int *);

   m=*m0;
   for(i=1;i<=m;i++) { w[i]=w0[i-1]; x[i]=x0[i-1];}
   for(i=1,k=0;i<=m;i++)
   {  for(j=1;j<=m;j++) { r[i][j]=corr[k]; k++;}}

/* for(i=1;i<=m;i++)
   {  for(j=1;j<=m;j++) printf("%8.4f", r[i][j]);
      printf("\n");
   }
   for(i=1;i<=m;i++) printf("%8.4f", w[i]);  printf("\n");
   for(i=1;i<=m;i++) printf("%8.4f", x[i]);  printf("\n");
*/

   /* first order binary approx */
   if(*type==1) mvn1(m,w,x,r,*nsim,*eps,&p,&perr,&ierr);
   /* second order binary approx */ 
   else  mvn2(m,w,x,r,*nsim,*eps,&p,&perr,&ierr);
   /* printf("appr. %d: %9.5f %9.5f %d\n", *type,p,perr,ierr);*/
   *pr=p; *sd=perr; *ifail=ierr;
}
