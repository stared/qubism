// Fitting to a general curve
// 100318


#include"fitting.h"

typedef struct
{
     Fit_FType F;
     const Vector* X;
     const Vector* Y;
} FitError_Params;

double Fit_Error(const Vector &D, void *p)
{
     FitError_Params *P=(FitError_Params *)p;
     Fit_FType F=P->F;
     const Vector *X=P->X;
     const Vector *Y=P->Y;
     long N=X->N;
     double error=0.0;
     for (long i=1;i<=N;i++)
	  error+=sqr( (*Y)(i) - F((*X)(i),D) );
     return error;
}

// A continuous broken straight line
double Fit_BrokenLine(double x, const Vector &D)
{
     double a1=D(1), b1=D(2), c=D(3), a2=D(4), b2=D(5);
     // double b2=(a1-a2)*c + b1;
     return (x<c ? a1*x+b1 : a2*x+b2);
}

double Fit_Line(double x, const Vector &D)
{
     double a=D(1), b=D(2);
     return a*x+b;
}

double Fit(Vector &D, Fit_FType F, const Vector &X, const Vector &Y)
{
     FitError_Params P;
     P.F=F;
     P.X=&X;
     P.Y=&Y;
     
//     Annealing_Optimize(D,Fit_Error,&P,NULL);
     return CG_Optimize(D,Fit_Error,&P,1e-5); // tolerance out?
}

