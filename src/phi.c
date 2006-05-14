#include <math.h>

double phi(double z)
{ return( 0.3989422804014327*exp(-.5*z*z)); }

double phi2(double z1, double z2, double rh)
{ double r1,tem;
  r1=1.-rh*rh;
  tem=(z1*z1+z2*z2-2.*rh*z1*z2)/r1;
  return( 0.1591549430918953*exp(-.5*tem)/sqrt(r1));
}

