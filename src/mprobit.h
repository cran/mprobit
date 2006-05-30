
#define UB 6.
#define EPS 1.e-7
#define MAXP 24
#define S 12
/* M-1 is an upper bound on the dimension 
   N is an upper bound on dimension for linear system solver
   for mvn1, N>=M;
   for mvn2, N>=M*(M-1)/2
*/
#define M 20
#define N 300

int mm,kk,ksign;
double *ww,*xx,rs,r1,r32;

double **x;     
int nn,nc,dmax,ncl,*y,*dvec,*dstart, d, norc;        
/* x    -- data matrix
   y    -- response variable
   nn   -- total number of individuals 
   nc   -- number of covariate
   dmax  --  max cluster size
   ncl  -- number of unique clusters 
   dvec -- a cluster size vector
   dstart -- the starting point of each unique cluster
   d    -- number of repeated measures of each culster
   norc   -- number of ordinal categories
*/
