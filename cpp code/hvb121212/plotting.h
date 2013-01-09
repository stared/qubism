// Plotting routines
// Part of the hvb++ project
// 100323

#ifndef PLOTTING
#define PLOTTING
#include"easyx.h"
#include"matrix.h"

enum Plotting_Styles {P_Points,P_Lines,P_Linepoints}; 

typedef struct
{
     bool open; // if true, no frame at all! build!
     bool frozen; // if true, do not modify!
     double xmin, xmax;
     double ymin, ymax;
     palette color;
     Plotting_Styles PS;
} Plotting_Frame;

extern Plotting_Frame *Current_Plotting_Frame;

void Start_Current_Frame();
Plotting_Frame *Get_Frame(const Vector &X, const Vector &Y);
void Adapt_Frame(Plotting_Frame *F, const Vector &X, const Vector &Y);
void Set_Plotting_Color(const palette color);
void Set_Plotting_Color(const char *colorname);
void Set_Plotting_Line(const Plotting_Styles PS);

void Plot_Using_Frame(const Vector &X, const Vector &Y, 
		      const Plotting_Frame *F);
void Plot(const Vector &X, const Vector &Y);

void Plot_Curves(const Vector &X, const Matrix &Y, const palette* P);

#endif
