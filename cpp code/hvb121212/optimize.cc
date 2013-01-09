// Optimization routines
// 100317
#include"optimize.h"

// This is a "metafunction"
// FuncGrad returns the value of a function and its gradient
// The first arguments are the outputs. Caution: G must be allocated!
double FuncGrad(Vector &G, const Vector &X, 
		Opt_VectorFType Function, void *params)
{
     double F=Function(X,params);
     double dx=1e-8;
     Vector X2=X;
     for (long i=1;i<=X.N;i++)
     {
	  double x=X2(i);
          X2(i)=x+dx;
          double F2=Function(X2,params);
	  G(i)=(F2-F)/dx;
          X2(i)=x;
     }
     return F;
}

typedef struct 
{
     Opt_VectorFType Function;
     Vector X0;
     Vector G;
     void *params;
}Func1Dparams_type;
// Needed to convert a vector function to 1d function for line search


// Another "metafunction"
// 1D-mensionalization of a vector function along a line
double Func1D(double x, void *Q)
{
     Func1Dparams_type Qex=*((Func1Dparams_type*)Q);
     return Qex.Function(Qex.X0+x*Qex.G,Qex.params);
}

double Line_Optimize(Vector &X, const Vector &G, 
		     Opt_VectorFType Function, void *params, double tol)
{
     // First, prepare parameters for Func1D
     Func1Dparams_type Qex;
     Qex.Function=Function;
     Qex.X0=X;
     Qex.G=G;
     Qex.params=params;
     
     double a=0.0,c=G.Norm(),b=c/2.0;
     Bracket_Minimum(a,b,c,Func1D,&Qex);
//     printf("minimum bracketed: %g %g %g\n",a,b,c);
     double x;
     double f=Brent_Optimize(x,a,b,c,Func1D,&Qex,tol);
//     printf("minimum found: %g\n",x);
     X=X+x*G;
     return f;
}

#define shift(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
void Bracket_Minimum(double &a, double &b, double &c, Opt_ScalarFType F,
		     void *p)
{
     static double golden=(1.0+sqrt(5.0))/2.0; 
     static double glimit=100.0;
     static double epsilon=1e-20;
     
     double ulim,u,r,q,fu;

     double fa=F(a,p);
     double fb=F(b,p);
     if (fb > fa) 
     {
	  SWAP(a,b);
	  SWAP(fa,fb);
     }
     c=b+golden*(b-a);
     double fc=F(c,p);
     while (fb > fc) 
     {
	  r=(b-a)*(fb-fc);
	  q=(b-c)*(fb-fa);
	  u=b-((b-c)*q-(b-a)*r)/
	       (2.0*SIGN(MAX(fabs(q-r),epsilon),q-r));
	  ulim=b+glimit*(c-b);
	  if ((b-u)*(u-c) > 0.0) 
	  {
	       fu=F(u,p);
	       if (fu < fc) 
	       {
		    a=b;
		    b=u;
		    fa=fb;
		    fb=fu;
		    return;
	       } 
	       else if (fu > fb) 
	       {
		    c=u;
		    fc=fu;
		    return;
	       }
	       u=c+golden*(c-b);
	       fu=F(u,p);
	  } 
	  else if ((c-u)*(u-ulim) > 0.0) 
	  {
	       fu=F(u,p);
	       if (fu < fc) 
	       {
		    shift(b,c,u,c+golden*(c-b));
		    shift(fb,fc,fu,F(u,p));
	       }
	  } else if ((u-ulim)*(ulim-c) >= 0.0) 
	  {
	       u=ulim;
	       fu=F(u,p);
	  } else 
	  {
	       u=c+golden*(c-b);
	       fu=F(u,p);
	  }
	  shift(a,b,c,u);
	  shift(fa,fb,fc,fu);
     }

}


double Brent_Optimize(double &x, 
		      double ax, double bx, double cx, 
		      Opt_ScalarFType F, void *p, double tol)
{
     static double cgold=(3.0-sqrt(5.0))/2.0;
     static long itmax= 200;
     static double epsilon=1e-10;     

     double a,b,etemp,fu,f3,f2,fx,tol1,tol2,u,x3,x2,xm;
     double d=0.0,e=0.0;
     
     a=(ax < cx ? ax : cx);
     b=(ax > cx ? ax : cx);
     x=x2=x3=bx;
     f2=f3=fx=F(x,p);
     long iter=0;
     do
     {
	  xm=0.5*(a+b);
	  tol2=2.0*(tol1=tol*fabs(x)+epsilon);
	  if (fabs(x-xm) <= (tol2-0.5*(b-a))) 			
	       return fx;
	  
	  if (fabs(e) > tol1) 
	  {
	       double r=(x-x2)*(fx-f3);
	       double q=(x-x3)*(fx-f2);
	       double p=(x-x3)*q-(x-x2)*r;
	       q=2.0*(q-r);
	       if (q > 0.0) p = -p;
	       q=fabs(q);
	       etemp=e;
	       e=d;
	       if (fabs(p) >= fabs(0.5*q*etemp) || 
		   p <= q*(a-x) || p >= q*(b-x))
		    d=cgold*(e=(x >= xm ? a-x : b-x));
	       else 
	       {
		    d=p/q;
		    u=x+d;
		    if (u-a < tol2 || b-u < tol2)
			 d=SIGN(tol1,xm-x);
	       }
	  } 
	  else 
	  {
	       d=cgold*(e=(x >= xm ? a-x : b-x));
	  }
	  u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
	  fu=F(u,p);
	  if (fu <= fx) 
	  {
	       if (u >= x) a=x; else b=x;
	       shift(x3,x2,x,u)
	       shift(f3,f2,fx,fu)		       
	  } 
	  else 
	  {
	       if (u < x) a=u; else b=u;
	       if (fu <= f2 || x2 == x) 
	       {
		    x3=x2;
		    x2=u;
		    f3=f2;
		    f2=fu;
	       } 
	       else if (fu <= f3 || x3 == x || x3 == x2) 
	       {
		    x3=u;
		    f3=fu;
	       }
	  }
	  iter++;
     }while(iter<itmax);
     printf("Too many iterations in brent");
     return fx;
}	

#undef shift

double CG_Optimize(Vector &X, Opt_VectorFType F, void *p, 
		   double tol)
{
     static long itmax=300;
     static double epsilon=1e-10;
     long N=X.N;
     
     double gg,fp,dgg,fret;
     
     Vector G(N);
     Vector H(N);
     Vector Xi(N);
     
     fp=FuncGrad(Xi,X,F,p);
     G=(-1.0)*Xi;
     Xi=H=G;
     for (long its=1;its<=itmax;its++) 
     {
	  fret=Line_Optimize(X,Xi,F,p,tol); 
	  if (2.0*fabs(fret-fp) <= tol*(fabs(fret)+fabs(fp)+epsilon)) 
	       return fret;
	  fp=FuncGrad(Xi,X,F,p);
	  gg=Dot(G,G);
	  dgg=Dot(Xi+G,Xi);
	  if (gg==0.0) return fp;
	  double gam=dgg/gg;
	  G=(-1.0)*Xi;
	  Xi=H=G+gam*H;
     }
     merror("Too many iterations in CG_Optimize\n");
     return fp;

}


double Annealing_Optimize(Vector &X, Opt_VectorFType F, void *p,
			  void *annealing_params)
{
     double betamin=1.0;
     double betamax=1000.0;
     double rbeta=1.02;
     long ntimes=1000;
     double A0=0.001;

     long N=X.N;
     
     if (annealing_params!=NULL)
     {
	  Annealing_Params *Q=(Annealing_Params*)annealing_params;
	  betamin=Q->betamin;
	  betamax=Q->betamax;
	  rbeta=Q->rbeta;
	  ntimes=Q->ntimes;
	  A0=Q->A0;
     }

     double f=F(X,p);
     double fbest=f;
     Vector Xbest=X;

     for (double beta=betamin;beta<=betamax;beta*=rbeta)
     {
	  long naccepted=0;
	  for (long nt=1;nt<=ntimes;nt++)
	  {
	       Vector X2(X);
	       long idx=rand_int(1,N);
	       X2(idx)+=A0*rand_double(-1.0,1.0);
	       double f2=F(X2,p);
	       if (f2<fbest)
	       {
		    fbest=f2;
		    Xbest=X2;
	       }
	       if (f2<f || rand_double()<exp(-beta*(f2-f)))
	       {
		    X=X2;
		    f=f2;
		    naccepted++;
	       }
	       if (nt == ntimes/10) // make it a parameter?
	       {
		    double acceptance_ratio=(double)naccepted/(double)(nt);
		    // printf("Acc: %g\n",acceptance_ratio);
		    if (acceptance_ratio<0.4) A0/=1.2;
		    if (acceptance_ratio>0.6) A0*=1.2;
	       }
	  }
     }
     X=Xbest;
     return fbest;
     
    // printf("Final epsilon: %g\n",eps);
     //error=Conjugate_Fit(D,X,Y,f);
     //return error;
}
