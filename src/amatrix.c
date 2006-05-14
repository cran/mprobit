#include <malloc.h>
#include <stdlib.h>

/* allocate matrices such that rows occupy contiguous space */
double **dmatrix(int nrow, int ncol)
{ int i; double **amat,*avec;
  avec=(double *)malloc((unsigned) (nrow*ncol)*sizeof(double));
  amat=(double **)malloc((unsigned) nrow*sizeof(double*));
  for(i=0;i<nrow;i++) amat[i]=avec+i*ncol;
  //free(avec);
  return amat;
}

int **imatrix(int nrow, int ncol)
{ int i; int **amat,*avec;
  avec=(int *)malloc((unsigned) (nrow*ncol)*sizeof(int));
  amat=(int **)malloc((unsigned) nrow*sizeof(int*));
  for(i=0;i<nrow;i++) amat[i]=avec+i*ncol;
  //free(avec);
  return amat;
}

/* 3-dim integer array */
/*int **iarray3(int n1, int n2, int n3)
{ int i,j; int ***arr,*avec;
  avec=(int *)malloc((unsigned) (n1*n2*n3)*sizeof(int));
  arr=(int ***)malloc((unsigned) n1*sizeof(int**));
  for(i=0;i<n1;i++) arr[i]=avec+i*n2*n3;
  for(i=0;i<n1;i++) 
  { for(j=0;j<n2;j++) arr[i][j]=avec+i*n2*n3+j*n3; }
  return arr;
}*/

/* 3-dim double array */
/*double **darray3(int n1, int n2, int n3)
{ int i,j; double ***arr,*avec;
  avec=(double *)malloc((unsigned) (n1*n2*n3)*sizeof(double));
  arr=(double ***)malloc((unsigned) n1*sizeof(double**));
  for(i=0;i<n1;i++) arr[i]=avec+i*n2*n3;
  for(i=0;i<n1;i++) 
  { for(j=0;j<n2;j++) arr[i][j]=avec+i*n2*n3+j*n3; }
  return arr;
}*/

