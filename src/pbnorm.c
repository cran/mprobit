#include <math.h>
/* 
   Bivariate normal routine from: 
   Donnelly, T.G. [1973],  Algorithm 462: bivariate normal distribution, 
   Communications of the association for computing machinery, 16, 638. 
*/
double bivnor(double ah, double ak, double r)
{  double twopi,b,gh,gk,rr,con,sqr,wh,wk,h2,a2,ex,h4,w2,sn,sp,ap,cn,t;
   double gw,g2,s2,s1,conex,sgn;
   double alnorm(double, int);
   int idig,i,is;
   twopi=6.283185307179587;
   b=0.;
   idig=9;
   gh=alnorm(ah,1)/2.; gk=alnorm(ak,1)/2.;
   if(r!=0.0)  rr=1.-r*r;
   else 
   { b=4.*gh*gk; if(b<0.) b=0.; if(b>1.) b=1.; return(b); }
   if(rr<0.) { return (-1.);}
   if(rr==0.0)  
   { if(r<0.) /* r=-1 */
     { if(ah+ak<0.) b=2.*(gh+gk)-1.;
       if(b<0.) b=0.; if(b>1.) b=1.; return(b); 
     }
     else /* r=1 */
     {  if(ah-ak<0.) b=2.*gk;
        else b=2.*gh;
        if(b<0.) b=0.; if(b>1.) b=1.; return(b);
     } 
   }
   /* rr != 0 , r!=0 */
   sqr=sqrt(rr);
   con=twopi*.5;
   for(i=1;i<=idig;i++) con/=10.;
   if(ah==0.0) 
   { if(ak==0.) 
     { b=atan(r/sqr)/twopi+.25; if(b<0.) b=0.; if(b>1.) b=1.; return(b); }
     else
     {  b+=gk;
        if(ah!=0.) { wh=-ah; wk=(ak/ah-r)/sqr; gw=2.*gh; is=-1;}
        else { wh=-ak; wk=(ah/ak-r)/sqr; gw=2.*gk; is=1;}
        goto L210;
     }
   }
   
      b=gh;
      if(ah*ak<0) { b-=.5;} 
      if(ah*ak!=0.) 
      { b+=gk;
        if(ah!=0.) { wh=-ah; wk=(ak/ah-r)/sqr; gw=2.*gh; is=-1;}
        else { wh=-ak; wk=(ah/ak-r)/sqr; gw=2.*gk; is=1;}
      }
      else { wh=-ah; wk=(ak/ah-r)/sqr; gw=2.*gh; is=-1;}
 L210: sgn=-1.; t=0.;
      if(wk!=0.) 
      { if(fabs(wk)==1.) { t=wk*gw*(1.-gw)*.5; goto L310;}
        else 
        { if(fabs(wk)>1.)
	  { sgn=-sgn; wh=wh*wk; g2=alnorm(wh,0); wk=1./wk;
	    if (wk<0) { b=b+.5;}
	    b=b-(gw+g2)*.5+ gw*g2;
	  }
	  h2=wh*wh; a2=wk*wk; h4=h2*.5; ex=0.0; if(h4<87.0) ex=exp(-h4);
	  w2=h4*ex; ap=1.; s2=ap-ex; sp=ap; s1=0.; sn=s1; conex=fabs(con/wk);
	  goto L290;
	}
      }
      else 
      {  goto L320; }

 L290: cn=ap*s2/(sn+sp); s1+=cn;
      while(fabs(cn)-conex>0) 
      {   sn=sp; sp+=1.; s2-=w2;
          if(fabs(w2)<=1.0e-15 || fabs(h4)<=1.0e-15) w2=0.0;
	  else { w2*=h4/sp;}
	  /*    underflow prevention   */
	  if(fabs(ap)<=1.0e-15 || fabs(a2)<=1.0e-15) { ap=0.0;}
	  else { ap=-ap*a2; }
          //goto L290;
	  cn=ap*s2/(sn+sp); s1+=cn;
      }
      //else
      t=(atan(wk)-wk*s1)/twopi; //goto L310;}
 L310: b+=sgn*t;
 L320: if(is<0) 
       {  if(ak!=0.) 
          {  wh=-ak; wk=(ah/ak-r)/sqr; gw=2.*gk; is=1; goto L210; }
          else { if(b<0.) b=0.; if(b>1.) b=1.; return(b); }
       }
       else { if(b<0.) b=0.; if(b>1.) b=1.; return(b); }
    
}

double  alnorm(double x, int upper)
/* algorithm as 66 by i.d. hill */
{  int up;
   double ltone,utzero,con,prob,z,y;
   con=1.28; ltone=5.0; utzero=12.5;
   up=upper; z=x;
   if(z<0.0) { up=1-up; z=-z;}
   if(z<=ltone || (up==1 && z<=utzero)) 
   { y=0.5*z*z;
     if(z<=con) 
     { prob=0.5-z*(0.398942280444-0.399903438504*y/ 
         (y+5.75885480458-29.8213557808/ 
         (y+2.62433121679+48.6959930692/ 
         (y+5.92885724438))));
     }
     else
     { prob=0.398942280385*exp(-y)/ 
         (z-3.8052e-8+1.00000615302/ 
         (z+3.98064794e-4+1.98615381364/ 
         (z-0.151679116635+5.29330324926/ 
         (z+4.8385912808-15.1508972451/ 
         (z+0.742380924027+30.789933034/(z+3.99019417011))))));
     }
   }
   else { prob = 0.0;}
   if(up==0) prob=1.0-prob;
   return(prob);
}
