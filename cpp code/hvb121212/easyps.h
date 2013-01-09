////////////////////////////////////////////////////////
// hvb++ 1.0
// Copyleft: Javier Rodríguez Laguna
// 100119

// Simple driver to create PostScript graphics from C
// Headers
// JaviRL, 0004, 0208, 0505, 0607, 0807
#ifndef EASYPS_HEADER
#define EASYPS_HEADER
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

FILE *PS_Open(const char *name, int x0, int y0, int x1, int y1);

void PS_Close(FILE *psfile);

void PS_Line(FILE *psfile, double x0, double y0, double x1, double y1);

void PS_Circle(FILE *psfile, double cx0, double cy0, double R);
     
void PS_FillCircle(FILE *psfile, double cx0, double cy0, double R, 
		     double filling);

void PS_ColorCircle(FILE *psfile, double cx0, double cy0, double R,
		    double red, double green, double blue);

void PS_Rectangle(FILE *psfile, double x0, double y0, double wx, double wy);

void PS_FillRectangle(FILE *psfile, double x0, double y0, 
			double wx, double wy, double filling);

void PS_ColorRectangle(FILE *psfile, double x0, double y0, 
		       double wx, double wy,
		       double red, double green, double blue);

void PS_NormalLine(FILE *psfile);

void PS_DashedLine(FILE *psfile);

void PS_SetLineWidth(FILE *psfile, double f);

void PS_Color(FILE *psfile, double red, double green, double blue);

extern const char *PS_DEFAULT_FONT;

// Remember: angle in degrees!!!!
void PS_Arc(FILE *psfile, double cx0, double cy0, double R,
	    double a0, double a1);

void PS_Curve(FILE *psfile, double x1, double y1, double x2, double y2,
	      double x3, double y3, double x4, double y4);

void PS_PrepareFont(FILE *psfile, const char *fontname, int size);

void PS_Text(FILE *psfile, double x0, double y0, const char *S);

#endif


