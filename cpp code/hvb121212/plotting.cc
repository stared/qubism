// Plotting routines
// Part of the hvb++ project
// 100323

#include"plotting.h"


Plotting_Frame *Current_Plotting_Frame=NULL;

void Start_Current_Frame()
{
     if (Current_Plotting_Frame!=NULL) return;
     Current_Plotting_Frame=(Plotting_Frame*)malloc(sizeof(Plotting_Frame));
     Current_Plotting_Frame->open=false;
     Current_Plotting_Frame->frozen=false;
}

Plotting_Frame *Get_Frame(const Vector &X, const Vector &Y)
{
     Plotting_Frame *F=(Plotting_Frame*)malloc(sizeof(Plotting_Frame));
     F->open=true;
     F->frozen=false;
     F->xmin=X.Min();
     F->xmax=X.Max();
     F->ymin=Y.Min();
     F->ymax=Y.Max();
     F->color=EXAllocNamedColor("white");
     F->PS=P_Points;
     return F;
}

void Adapt_Frame(Plotting_Frame *F, const Vector &X, const Vector &Y)
{
     if (F->frozen) return;
     if (F->open)
     {
	  F->xmin=MIN(F->xmin,X.Min());
	  F->xmax=MAX(F->xmax,X.Max());
	  F->ymin=MIN(F->ymin,Y.Min());
	  F->ymax=MAX(F->ymax,Y.Max());
     }
     else 
     {
	  F->xmin=X.Min();
	  F->xmax=X.Max();
	  F->ymin=Y.Min();
	  F->ymax=Y.Max();
	  F->open=true;
     }
}

void Set_Plotting_Color(const palette color)
{
     if (!Current_Plotting_Frame) Start_Current_Frame();
     Current_Plotting_Frame->color=color;
}

void Set_Plotting_Color(const char *colorname)
{
     if (!Current_Plotting_Frame) Start_Current_Frame();
     Current_Plotting_Frame->color=EXAllocNamedColor(colorname);
}

void Set_Plotting_Line(const Plotting_Styles PS)
{
     if (!Current_Plotting_Frame) Start_Current_Frame();
     Current_Plotting_Frame->PS=PS;
}

void Plot(const Vector &X, const Vector &Y)
{
     if (!Current_Plotting_Frame)
	  Current_Plotting_Frame=Get_Frame(X,Y);
     else Adapt_Frame(Current_Plotting_Frame,X,Y);
     Plot_Using_Frame(X,Y,Current_Plotting_Frame);
     Current_Plotting_Frame->frozen=true;
}

void Plot_Using_Frame(const Vector &X, const Vector &Y, 
		      const Plotting_Frame *F)
{

     // Currently, only plot using lines... we'll do more in future
     EXSetColor(F->color);
     long width, height;
     EXGetWindowSize(width,height);
     long N=X.N;
     double Dx=F->xmax-F->xmin;
     double Dy=F->ymax-F->ymin;

     double xnold=(X(1)-F->xmin)/Dx;
     double ynold=(Y(1)-F->ymin)/Dy;
     int pxold=(int)floor(xnold*(double)width);
     int pyold=height-(int)floor(ynold*(double)height);

     for (long i=2;i<=N;i++)
     {
	  double xn=(X(i)-F->xmin)/Dx;
	  double yn=(Y(i)-F->ymin)/Dy;
	  int px=(int)floor(xn*(double)width);
	  int py=height-(int)floor(yn*(double)height);
//	  printf("%d %d\n",px,py);
	  EXLine(pxold,pyold,px,py);
	  pxold=px; pyold=py;
     }
     EXFlush();
}

void Plot_Curves(const Vector &X, const Matrix &Y, const palette* P)
{
     // First, get frame for all curves
     if (Current_Plotting_Frame) free(Current_Plotting_Frame);
     Current_Plotting_Frame=Get_Frame(X,Y.Col(1));
     for (long i=2;i<=Y.N2;i++)
	  Adapt_Frame(Current_Plotting_Frame,X,Y.Col(i));
     Current_Plotting_Frame->frozen=true;
     
     for (long i=1;i<=Y.N2;i++)
     {
	  Set_Plotting_Color(P[i]);
	  Plot(X,Y.Col(i));
     }
}

