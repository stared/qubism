////////////////////////////////////////////////////////
// hvb++ 1.0
// Copyleft: Javier Rodríguez Laguna
// 080725

// Simple driver to create PostScript graphics from C
// Headers
// JaviRL, 0004, 0208, 0505, 0607
#ifndef EASYPS
#define EASYPS
#include"easyps.h"


FILE *PS_Open(const char *name,int x0, int y0, int x1, int y1)
{
     FILE *psfile;
     psfile=fopen(name,"wt");
     if (psfile==NULL) 
     {
	  printf("Error opening file %s\n",name);
	  exit(1);
     }
     fprintf(psfile,"%%!PS-Adobe-2.0\n");
     fprintf(psfile,"%%%%Title: %s\n",name);
     fprintf(psfile,"%%%%Creator: EasyPS version 100119\n");
     fprintf(psfile,"%%%%BoundingBox: %d %d %d %d\n",x0,y0,x1,y1);
     fprintf(psfile,"%%%%EndComments\n");
     fprintf(psfile,"/n {newpath} def /m {moveto} def /l {lineto} def\n");
     fprintf(psfile,"/s {stroke} def /sl {setlinewidth} def \n");
     fprintf(psfile,"/sr {setrgbcolor} def /sg {setgray} def /f {fill} def\n");
     fprintf(psfile,"/a {0 360 arc} def /c {closepath} def ");
     fprintf(psfile,"/sd {setdash} def\n");
     fprintf(psfile,"0.2 sl\n");
     return psfile;
}

void PS_Close(FILE *psfile)
{
     fprintf(psfile,"showpage\n");
     fclose(psfile);
}

void PS_Line(FILE *psfile, double x0, double y0, double x1, double y1)
{
     fprintf(psfile,"n\n");
     fprintf(psfile,"%.4g %.4g m\n",x0,y0);
     fprintf(psfile,"%.4g %.4g l\n",x1,y1);
     fprintf(psfile,"s\n");
}

void PS_Circle(FILE *psfile, double cx0, double cy0, double R)
{
     fprintf(psfile,"n\n");
     fprintf(psfile,"%.4g %.4g %.4g a\n",cx0,cy0,R);
     fprintf(psfile,"s\n");
}

void PS_FillCircle(FILE *psfile, double cx0, double cy0, double R, 
		     double filling)
{
     fprintf(psfile,"n\n");
     fprintf(psfile,"%.4g %.4g %.4g a\n",cx0,cy0,R);
     fprintf(psfile,"%.4g sg\n",filling);
     fprintf(psfile,"f\n");
     fprintf(psfile,"0 sg\n");
}

void PS_ColorCircle(FILE *psfile, double cx0, double cy0, double R,
		    double red, double green, double blue)
{
     fprintf(psfile,"n\n");
     fprintf(psfile,"%.4g %.4g %.4g a\n",cx0,cy0,R);
     fprintf(psfile,"%.4g %.4g %.4g sr\n",red, green, blue);
     fprintf(psfile,"f\n");
     fprintf(psfile,"0 0 0 sr\n"); 
}

void PS_Rectangle(FILE *psfile, double x0, double y0, double wx, double wy)
{
     fprintf(psfile,"n\n");
     fprintf(psfile,"%.4g %.4g m\n",x0,y0);
     fprintf(psfile,"%.4g %.4g l\n",x0+wx,y0);
     fprintf(psfile,"%.4g %.4g l\n",x0+wx,y0+wy);
     fprintf(psfile,"%.4g %.4g l\n",x0,y0+wy);
     fprintf(psfile,"c\n");
     fprintf(psfile,"s\n");
}

void PS_FillRectangle(FILE *psfile, double x0, double y0, 
			double wx, double wy, double filling)
{
     fprintf(psfile,"n\n");
     fprintf(psfile,"%.4g %.4g m\n",x0,y0);
     fprintf(psfile,"%.4g %.4g l\n",x0+wx,y0);
     fprintf(psfile,"%.4g %.4g l\n",x0+wx,y0+wy);
     fprintf(psfile,"%.4g %.4g l\n",x0,y0+wy);
     fprintf(psfile,"c\n");
     fprintf(psfile,"%.4g sg\n",filling);
     fprintf(psfile,"f\n");
     fprintf(psfile,"0 sg\n");
}

void PS_ColorRectangle(FILE *psfile, double x0, double y0, 
		       double wx, double wy,
		       double red, double green, double blue)
{
     fprintf(psfile,"n\n");
     fprintf(psfile,"%.4g %.4g m\n", x0, y0);
     fprintf(psfile,"%.4g %.4g l\n", x0+wx, y0);
     fprintf(psfile,"%.4g %.4g l\n", x0+wx, y0+wy);
     fprintf(psfile,"%.4g %.4g l\n", x0, y0+wy);
     fprintf(psfile,"c\n");
     fprintf(psfile,"%.4g %.4g %.4g sr\n",red, green, blue);
     fprintf(psfile,"f\n");
     fprintf(psfile,"0 0 0 sr\n");
}

void PS_NormalLine(FILE *psfile)
{
     fprintf(psfile,"[] 0 sd\n");
}

void PS_DashedLine(FILE *psfile)
{
     fprintf(psfile,"[3 2] 0 sd\n");
}

void PS_SetLineWidth(FILE *psfile, double f)
{
     fprintf(psfile,"%.4g sl\n",f);
}

void PS_Color(FILE *psfile, double red, double green, double blue)
{
     fprintf(psfile,"%.4g %.4g %.4g sr\n", red, green, blue);

}

void PS_Arc(FILE *psfile, double cx0, double cy0, double R,
	    double a0, double a1)
{
     fprintf(psfile,"n\n");
     fprintf(psfile,"%.4g %.4g %.4g %.4g %.4g arc\n",cx0,cy0,R,a0,a1);
     fprintf(psfile,"s\n");
}

void PS_Curve(FILE *psfile, double x1, double y1, double x2, double y2,
	      double x3, double y3, double x4, double y4)
{
     fprintf(psfile,"n\n");
     fprintf(psfile,"%.4g %.4g m\n",x1,y1);
     fprintf(psfile,"%.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g curveto\n",
	     x1,y1,x2,y2,x3,y3,x4,y4);
     fprintf(psfile,"s\n");
}

const char *PS_DEFAULT_FONT="Times-Roman";

void PS_PrepareFont(FILE *psfile, const char *fontname, int size)
{
     fprintf(psfile,"/%s findfont\n",fontname);
     fprintf(psfile,"%d scalefont\n",size);
     fprintf(psfile,"setfont\n");
}

void PS_Text(FILE *psfile, double x0, double y0, const char *S)
{
     fprintf(psfile,"%.4g %.4g m\n",x0,y0);
     fprintf(psfile,"(%s) show\n",S);
}



#endif


