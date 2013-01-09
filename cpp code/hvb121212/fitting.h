// Fitting to a general curve
// 100318

#ifndef FITTING_HEADER
#define FITTING_HEADER

#include"optimize.h"

typedef double(*Fit_FType)(double x, const Vector &D);

double Fit_Error(Vector &D, void *p);

double Fit(Vector &D, Fit_FType F, const Vector &X, const Vector &Y);

double Fit_Line(double x, const Vector &D);
double Fit_BrokenLine(double x, const Vector &D);

#endif
