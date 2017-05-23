#include "stdio.h"
#include "nrutil.h"
#include "math.h"

#define TINY 1.0e-20
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void spline(float x[],float y[],int n,float yp1,float ypn,float y2[])
{
	int i,k;
	float p,qn,sig,un,*u;
	
	u=vector(1,n-1);
	if(yp1>0.99e30) y2[1]=u[1]=0.0;
	else
	{
		y2[1]= -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for(i=2;i<=n-1;i++)
	{
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if(ypn>0.99e30) qn=un=0.0;
	else
	{
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for(k=n-1;k>=1;k--) y2[k]=y2[k]*y2[k+1]+u[k];
	free_vector(u,1,n-1);
}

void splint(float xa[],float ya[],float y2a[],int n,float x,float *y)
{
	int klo,khi,k;
	float h,b,a;
	
	klo=1;
	khi=n;
	while(khi-klo>1)
	{
		k=(khi+klo)>>1;
		if(xa[k]>x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if(h==0.0) nrerror("Bad xa input to routine splint");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void ludcmp(double **a,int n,int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;
	
	vv=dvector(1,n);
	*d=1.0;
	for(i=1;i<=n;i++)
	{
		big=0.0;
		for(j=1;j<=n;j++) if((temp=fabs(a[i][j]))>big) big=temp;
		if(big==0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for(j=1;j<=n;j++)
	{
		for(i=1;i<j;i++)
		{
			sum=a[i][j];
			for(k=1;k<i;k++) sum-=a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for(i=j;i<=n;i++)
		{
			sum=a[i][j];
			for(k=1;k<j;k++) sum-=a[i][k]*a[k][j];
			a[i][j]=sum;
			if((dum=vv[i]*fabs(sum))>=big)
			{
				big=dum;
				imax=i;
			}
		}
		if(j!=imax)
		{
			for(k=1;k<=n;k++)
			{
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d= -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if(a[j][j]==0.0) a[j][j]=TINY;
		if(j!=n)
		{
			dum=1.0/(a[j][j]);
			for(i=j+1;i<=n;i++) a[i][j]*=dum;
		}
	}
	free_dvector(vv,1,n);
}

void lubksb(double **a,int n,int *indx,double b[])
{
	int i,ii=0,ip,j;
	double sum;
	
	for(i=1;i<=n;i++)
	{
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if(ii) for(j=ii;j<=i-1;j++) sum-=a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for(i=n;i>=1;i--)
	{
		sum=b[i];
		for(j=i+1;j<=n;j++) sum-=a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}


void invert_matrix(double **a, double **y, int n)
{
   int *indx;
   double *col, d;
   int i, j;
   col = dvector(1,n);
   indx = ivector(1,n);
   ludcmp(a,n,indx,&d);
   for (j=1;j<=n;j++)
   {
      for (i=1; i<=n; i++) col[i] = 0.0;
      col[j] = 1.0;
      lubksb(a,n,indx,col);
      for (i=1; i<=n; i++) y[i][j]=col[i];
   }
   free_dvector(col, 1,n);
   free_ivector(indx, 1,n);
}



#define NMAX 50000

void amoeba(double **p,double y[],int ndim, double ftol, double(*funk)(double []),int *nfunk, int limit)
{
	double amotry(double **p,double y[],double psum[],int ndim,double (*funk)(double []),int ihi, double fac);
	int i,ihi,ilo,inhi,j,mpts=ndim+1;
	double rtol,sum,swap,ysave,ytry,*psum;
	
	psum=dvector(1,ndim);
	*nfunk=0;
	for(j=1;j<=ndim;j++)
	{
		for(sum=0.0,i=1;i<=mpts;i++) sum+=p[i][j];
		psum[j]=sum;
	}
	for(;;)
	{
                if (*nfunk > limit)
                {
                        printf("Amoeba: limit to evaluations reached %d\n", limit);
   		        break;
                }
		ilo=1;
		ihi=y[i]>y[2]?(inhi=2,1):(inhi=1,2);
		for(i=1;i<=mpts;i++)
		{
			if(y[i]<=y[ilo]) ilo=i;
			if(y[i]>y[ihi])
			{
				inhi=ihi;
				ihi=i;
			}
			else if (y[i]>y[inhi]&&i!=ihi) inhi=i;
		}
/*		printf("%.7lf ",y[ilo]);
		for(i=1;i<=ndim;i++) printf("%.4lf ",p[ilo][i]);
		printf("\n"); */
		rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
		if(rtol<ftol)
		{
			SWAP(y[1],y[ilo])
			for(i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i])
			break;
		}
		if(*nfunk>=NMAX) nrerror("NMAX exceeded");
		*nfunk+=2;
		ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0);
		if(ytry<=y[ilo]) ytry=amotry(p,y,psum,ndim,funk,ihi,2.0);
		else if(ytry>=y[inhi])
		{
			ysave=y[ihi];
			ytry=amotry(p,y,psum,ndim,funk,ihi,0.5);
			if(ytry>=ysave)
			{
				for(i=1;i<=mpts;i++)
				{
					if(i!=ilo)
					{
						for(j=1;j<=ndim;j++) p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
						y[i]=(*funk)(psum);
					}
				}
				*nfunk+=ndim;
				for(j=1;j<=ndim;j++)
				{
					for(sum=0.0,i=1;i<=mpts;i++) sum+=p[i][j];
					psum[j]=sum;
				}
			}
		}
		else --(*nfunk);
	}
	free_dvector(psum,1,ndim);
}

double amotry(double **p,double y[],double psum[],int ndim,double (*funk)(double []),int ihi,double fac)
{
	int j;
	double fac1,fac2,ytry,*ptry;

	ptry=dvector(1,ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for(j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	ytry=(*funk)(ptry);
	if(ytry<y[ihi])
	{
		y[ihi]=ytry;
		for(j=1;j<=ndim;j++)
		{
			psum[j]+=ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	free_dvector(ptry,1,ndim);
	return(ytry);
}



/* maximizebysimplex-------------------------------------------------------
Parameters. n = number of variable parameters
            param_vector = vector containing the n parameter values
            maxevals = maximum evaluations of function allowed
            funk = function taking a vector as parameter. This contains
                   the current values of the n parameters which must be unloaded
                   into their corresponding variables
            fac  = initial factor for step size (e.g., 1.1)
            convcrit = convergence criterion, ie. 1e10-9
--------------------------------------------------------------------------*/
double maximizebysimplex(int n, double *param_vector,
   int maxevals,  double(*funk)(double []), double fac,
   double convcrit)
{
   double **p,*y,min;
   int i,j,nfunc;

/*   printf("No. variable parameters %d Max. no. evaluations %d\n", n, maxevals);*/
   p=dmatrix(1,n+1,1,n);
   y=dvector(1,n+1);

   for (i=1; i<=n; i++)
   {
     p[1][i] = param_vector[i];
/*     printf("p[1][%d] = %lf\n", i, p[1][i]);*/
   }

   y[1] = funk(p[1]);

   for(i=2;i<=(n+1);i++)
   {
      for(j=1;j<=n;j++) p[i][j]=p[1][j];
      p[i][i-1]*=fac;
      y[i]=funk(p[i]);
   }
   amoeba(p,y,n,convcrit,funk,&nfunc, maxevals);
   min=y[1];
   free_dmatrix(p,1,n+1,1,n);
   free_dvector(y,1,n+1);
   return(min);

}


void powell(double p[],double **xi,int n,double ftol,int *iter,double *fret,double(*func)(double []))
{
	void linmin(double p[],double xi[],int n,double *fret,double (*func)(double []));
	int i,ibig,j;
	double del,fp,fptt,t,*pt,*ptt,*xit;
	
	pt=dvector(1,n);
	ptt=dvector(1,n);
	xit=dvector(1,n);
	*fret=(*func)(p);
	for(j=1;j<=n;j++) pt[j]=p[j];
	for(*iter=1;;++(*iter))
	{
		fp=(*fret);
		ibig=0;
		del=0.0;
		for(i=1;i<=n;i++)
		{
			for(j=1;j<=n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			linmin(p,xit,n,fret,func);
			if(fabs(fptt-(*fret))>del)
			{
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
		if(2.0*fabs(fp-(*fret))<=ftol*(fabs(fp)+fabs(*fret)))
		{
			free_dvector(xit,1,n);
			free_dvector(ptt,1,n);
			free_dvector(pt,1,n);
			return;
		}
		if(*iter==NMAX) nrerror("powell exceeding maximum iterations.");
		for(j=1;j<=n;j++)
		{
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt);
		if(fptt<fp)
		{
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
			if(t<0.0)
			{
				linmin(p,xit,n,fret,func);
				for(j=1;j<=n;j++)
				{
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
/*		printf("%.7lf ",(*fret));
		for(i=1;i<=n;i++) printf("%.4lf ",p[i]);
		printf("\n"); */
	}
}

#define TOL 1.0e-8

int ncom;
double *pcom,*xicom,(*nrfunc)(double []);

void linmin(double p[],double xi[],int n,double *fret,double (*func)(double []))
{
	double brent(double ax,double bx,double cx,double (*f)(double),double tol,double *xmin);
	double f1dim(double x);
	void mnbrak(double *ax,double *bx,double *cx,double *fa,double *fb,double *fc,double (*func)(double));
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	pcom=dvector(1,n);
	xicom=dvector(1,n);
	nrfunc=func;
	for(j=1;j<=n;j++)
	{
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	for (j=1;j<=n;j++)
	{
		xi[j]*=xmin;
		p[j]+=xi[j];
	}
	free_dvector(xicom,1,n);
	free_dvector(pcom,1,n);
}

double f1dim(double x)
{
	int j;
	double f,*xt;
	
	xt=dvector(1,ncom);
	for(j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	free_dvector(xt,1,ncom);
	return(f);
}

#define GOLD 1.61803399
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(double *ax,double *bx,double *cx,double *fa,double *fb,double *fc,double (*func)(double))
{
	double ulim,u,r,q,fu,dum;
	
	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if(*fb>*fa)
	{
		SHFT(dum,*ax,*bx,dum);
		SHFT(dum,*fb,*fa,dum);
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while(*fb>*fc)
	{
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if((*bx-u)*(u-*cx)>0.0)
		{
			fu=(*func)(u);
			if(fu<*fc)
			{
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			}
			else if (fu>*fb)
			{
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		else if((*cx-u)*(u-ulim)>0.0)
		{
			fu=(*func)(u);
			if(fu<*fc)
			{
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		}
		else if((u-ulim)*(ulim-*cx)>=0.0)
		{
			u=ulim;
			fu=(*func)(u);
		}
		else
		{
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u);
		SHFT(*fa,*fb,*fc,fu);
	}
}

#define R 0.61803399
#define C (1.0 - R)
#define SHFT2(a,b,c) (a) = (b);(b)=(c);
#define SHFT3(a,b,c,d) (a) = (b);(b)=(c);(c)=(d);

double golden(double ax, double bx, double cx, double (*f)(double),
   double tol, double *xmin, int *neval)
{
   double f1, f2, x0, x1, x2, x3;
   x0 = ax;
   x3 = cx;
   if (fabs(cx-bx) > fabs(bx-ax))
   {
      x1 = bx;
      x2 = bx + C*(cx-bx);
   }
   else
   {
      x2 = bx;
      x1 = bx - C*(bx-ax);
   }
   f1 = (*f)(x1);
   f2 = (*f)(x2);
   *neval = 2;
   while (fabs(x3 - x0) > tol*(fabs(x1)+fabs(x2)))
   {
//      printf("x0 %lf x1 %lf x2 %lf x3 %lf\n", x0, x1, x2, x3);
//      printf("f1 %lf f2 %lf\n", f1, f2);
//      monitorinput();
      if (f2 < f1)
      {
         SHFT3(x0, x1, x2, R*x1+C*x3);
         SHFT2(f1, f2, (*f)(x2));
      }
      else
      {
         SHFT3(x3, x2, x1, R*x2+C*x0);
         SHFT2(f2, f1, (*f)(x1));
      }
      (*neval)++;
   }
   if (f1 < f2)
   {
      *xmin = x1;
      return f1;
   }
   else
   {
      *xmin = x2;
      return f2;
   }
}


#define ITMAX 500
#define CGOLD 0.3819660
#define ZEPS 1.0e-20

double brent(double ax,double bx,double cx,double (*f)(double),double tol,double *xmin)
{
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;
		
	a=(ax<cx?ax:cx);
	b=(ax>cx?ax:cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
	for(iter=1;iter<=ITMAX;iter++)
	{
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if(fabs(x-xm)<=(tol2-0.5*(b-a)))
		{
			*xmin=x;
			return fx;
		}
		if(fabs(e)>tol1)
		{
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if(q>0.0) p=-p;
			q=fabs(q);
			etemp=e;
			e=d;
			if(fabs(p)>=fabs(0.5*q*etemp)||p<=q*(a-x)||p>=q*(b-x)) d=CGOLD*(e=(x>=xm?a-x:b-x));
			else
			{
				d=p/q;
				u=x+d;
				if(u-a<tol2||b-u<tol2) d=SIGN(tol1,xm-x);
			}
		}
		else d=CGOLD*(e=(x>=xm?a-x:b-x));
		u=(fabs(d)>=tol1?x+d:x+SIGN(tol1,d));
		fu=(*f)(u);
		if(fu<=fx)
		{
			if(u>=x) a=x;
			else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		}
		else
		{
			if(u<x) a=u;
			else b=u;
			if(fu<=fw||w==x)
			{
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			}
			else if(fu<=fv||v==x||v==w)
			{
				v=w;
				fv=fu;
			}
		}
	}
	nrerror("Too many iterations in brent");
	*xmin=x;
	return fx;
}


float gammln(float xx)
{
   double x, y, tmp, ser;
   float res;
   static double cof[6]={76.18009172947146,-86.50532032941677,
      24.01409824083091, -1.231739572450155,
      0.1208650973866179e-2, -0.5395239384953e-5};
   int j;
   y=x=xx;
   tmp = x+5.5;
   tmp -= (x+0.5)*log(tmp);
   ser = 1.000000000190015;
   for (j=0; j<=5; j++) ser += cof[j]/++y;
   res = -tmp+log(2.5066282746310005*ser/x);
/*   printf("gammln(%f) %f\n", xx, res);*/
   return(res);
}


#include <math.h>


double qgaus(double (*func)(double), double a, double b)
{
   int j;
   double xr, xm, dx, s, temp;
   static double x[]={0.0,0.09501,0.28160,0.45802,0.61788,0.75540,0.86563,0.94458,0.98940};
   static double w[]={0.0,0.18945,0.18260,0.16916,0.14960,0.12463,0.09516,0.06225,0.02715};
   xm = 0.5*(b+a);
   xr = 0.5*(b-a);
   s = 0;
   for (j=1; j<=8; j++)
   {
      dx = xr*x[j];
      temp =  w[j]*((*func)(xm+dx)+(*func)(xm-dx));
/*      printf("qgaus: j %d, temp %f, xm+dx %f, xm-dx %f\n", j, temp, xm+dx, xm-dx);*/
      s += temp;
   }
   return s *= xr;
}



#define EPS 1.0e-5
#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5
#define FUNC(x) ((*func)(x))



double midpnt(double (*func)(double), double a, double b, int n)
{
   double x, tnm, sum, del, ddel;
   static double s;
   int it, j;
   if (n==1)
   {
      return (s=(b-a)*FUNC(0.5*(a+b)));
   }
   else
   {
      for (it=1, j=1; j<n-1; j++) it *= 3;
      tnm = it;
      del = (b-a)/(3.0*tnm);
      ddel = del + del;
      x = a+0.5*del;
      sum = 0.0;
      for (j=1; j<=it; j++)
      {
         sum += FUNC(x);
         x+= ddel;
         sum += FUNC(x);
         x += del;
      }
      s = (s+(b-a)*sum/tnm)/3.0;
      return s;
   }
}



void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
   int i, m, ns=1;
   double den, dif, dift, ho, hp, w;
   double *c, *d;
   dif = fabs(x-xa[1]);
   c = dvector(1, n);
   d = dvector(1, n);
   for (i=1; i<=n; i++)
   {
      if ((dift=fabs(x-xa[i])) < dif)
      {
         ns = i;
         dif = dift;
      }
      c[i] = ya[i];
      d[i] = ya[i];
   }
   *y = ya[ns--];
   for (m=1; m<n; m++)
   {
      for (i=1; i<=n-m; i++)
      {
         ho = xa[i] - x;
         hp = xa[i+m] - x;
         w = c[i+1] - d[i];
         if ((den=ho-hp) == 0.0) nrerror("Error in routine polint");
         den = w/den;
         d[i] = hp*den;
         c[i] = ho*den;
      }
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
   }
   free_dvector(d, 1, n);
   free_dvector(c, 1, n);
}





double qromo(double (*func)(double), double a, double b, double
   (*choose)(double(*)(double), double, double, int))
{
   int j;
   double ss, dss, h[JMAXP+1], s[JMAXP+1];
   h[1] = 1.0;
   for (j=1; j<=JMAX; j++)
   {
      s[j] = (*choose)(func, a, b, j);
      if (j>=K)
      {
         polint(&h[j-K], &s[j-K], K, 0.0, &ss, &dss);
         if (fabs(dss) < EPS*fabs(ss)) return ss;
      }
      s[j+1] = s[j];
      h[j+1] = h[j]/9.0;
   }
   gabort("Too many steps in qromo", 0);
   return 0.0;
}


double trapzd (double (*func)(double), double a, double b, int n)
{
   double x, tnm, sum ,del;
   static double s;
   int it, j;
   if (n==1)
   {
      return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
   }
   else
   {
      for (it =1, j=1; j<n-1; j++) it <<=1;
      tnm = it;
      del = (b-a)/tnm;
      x = a + 0.5*del;
      for (sum=0.0, j=1; j<=it; j++, x += del) sum += FUNC(x);
      s = 0.5*(s+(b-a)*sum/tnm);
      return s;
   }
}


double qromb(double (*func)(double), double a, double b)
{
   double ss, dss;
   double s[JMAXP+1], h[JMAXP+1];
   int j;
   
   h[1] = 1.0;
   for (j=1; j<=JMAX; j++)
   {
      s[j] = trapzd(func, a, b, j);
      if (j >= K)
      {
         polint(&h[j-K], &s[j-K], K, 0.0, &ss, &dss);
         if (fabs(dss) < EPS*fabs(ss)) return ss;
      }
      s[j+1] = s[j];
      h[j+1] = h[j]/9.0;
   }
   gabort("Too many steps in qromo", 0);
   return 0.0;
}


// Returns ln(n!)

double factln(int n)
{
//	double gammln(double xx);
	void nrerror(char error_text[]);
	static double a[101];

	if (n < 0) nrerror("Negative factorial in routine factln");
	if (n <= 1) return 0.0;
	if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
	else return gammln(n+1.0);
}

//Returns the binomial coefficient (n k)

double bico(int n, int k)
{
	double factln(int n);

	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

//--------------------------------------------------------------------------

#define GET_PSUM \
					for (n=1;n<=ndim;n++) {\
					for (sum=0.0,m=1;m<=mpts;m++) sum += p[m][n];\
					psum[n]=sum;}
double tt;

// NOTE - CONTAINS A PROBLEM WHEN USED WITHIN A LIBRARY!
void amebsa(double **p, double y[], int ndim, double pb[], double *yb,
   double ftol, double (*funk)(), int *iter, double temptr)
{
	double amotsa();
	int i,ihi,ilo,j,m,n,mpts=ndim+1;
	double rtol,sum,swap,yhi,ylo,ynhi,ysave,yt,ytry,*psum;

	psum=dvector(1,ndim);
	tt = -temptr;
	GET_PSUM
	for (;;) {
		ilo=1;
		ihi=2;
		ynhi=ylo=y[1]+tt*log(uniform());
		yhi=y[2]+tt*log(uniform());
		if (ylo > yhi) {
			ihi=1;
			ilo=2;
			ynhi=yhi;
			yhi=ylo;
			ylo=ynhi;
		}
		for (i=3;i<=mpts;i++) {
			yt=y[i]+tt*log(uniform());
			if (yt <= ylo) {
				ilo=i;
				ylo=yt;
			}
			if (yt > yhi) {
				ynhi=yhi;
				ihi=i;
				yhi=yt;
			} else if (yt > ynhi) {
				ynhi=yt;
			}
		}
		rtol=2.0*fabs(yhi-ylo)/(fabs(yhi)+fabs(ylo));
		if (rtol < ftol || *iter < 0) {
			swap=y[1];
			y[1]=y[ilo];
			y[ilo]=swap;
			for (n=1;n<=ndim;n++) {
				swap=p[1][n];
				p[1][n]=p[ilo][n];
				p[ilo][n]=swap;
			}
			break;
		}
		*iter -= 2;
		ytry=amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,-1.0);
		if (ytry <= ylo) {
			ytry=amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,2.0);
		} else if (ytry >= ynhi) {
			ysave=yhi;
			ytry=amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,0.5);
			if (ytry >= ysave) {
				for (i=1;i<=mpts;i++) {
					if (i != ilo) {
						for (j=1;j<=ndim;j++) {
							psum[j]=0.5*(p[i][j]+p[ilo][j]);
							p[i][j]=psum[j];
						}
						y[i]=(*funk)(psum);
					}
				}
				*iter -= ndim;
				GET_PSUM
			}
		} else ++(*iter);
	}
	free_dvector(psum,1,ndim);
}
#undef GET_PSUM


double amotsa(double **p, double y[], double psum[], int ndim, double pb[],
   double *yb, double (*funk)(), int ihi, double *yhi, double fac)
{
	int j;
	double fac1,fac2,yflu,ytry,*ptry;

	ptry=dvector(1,ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=1;j<=ndim;j++)
		ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	ytry=(*funk)(ptry);
	if (ytry <= *yb) {
		for (j=1;j<=ndim;j++) pb[j]=ptry[j];
		*yb=ytry;
	}
	yflu=ytry-tt*log(uniform());
	if (yflu < *yhi) {
		y[ihi]=ytry;
		*yhi=yflu;
		for (j=1;j<=ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	free_dvector(ptry,1,ndim);
	return yflu;
}
