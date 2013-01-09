// Optimization routines
// JaviRL, 100317
#ifndef OPTIMIZE_HEADER
#define OPTIMIZE_HEADER

#include"matrix.h"

typedef double(*Opt_VectorFType)(const Vector &, void *);
typedef double(*Opt_ScalarFType)(double x, void *);


// This is a "metafunction"
// FuncGrad returns the value of a function and its gradient
// The first arguments are the outputs. Caution: G must be allocated!
double FuncGrad(Vector &G, const Vector &X, Opt_VectorFType Function,
     void *params);

double Line_Optimize(Vector &X, const Vector &G, Opt_VectorFType Function,
		     void *params, double tol);

void Bracket_Minimum(double &a, double &b, double &c, Opt_ScalarFType F,
		     void *params);

double Brent_Optimize(double &x, double a, double b, double c, 
		      Opt_ScalarFType F, void *params, double tol);

double CG_Optimize(Vector &X, Opt_VectorFType Function, void *params, 
		   double tol);

typedef struct 
{
     double betamin;
     double betamax;
     double rbeta;
     long ntimes;
     double A0;
} Annealing_Params;

double Annealing_Optimize(Vector &X, Opt_VectorFType Function, void *params,
			  void *annealing_params);
#endif
