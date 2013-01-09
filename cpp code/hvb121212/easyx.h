////////////////////////////////////////////////////////
// hvb++ 1.0
// Copyleft: Javier Rodríguez Laguna
// 080725

// easyx.h version 1.0
// JaviRL, Apr 2006; Aug 2006

#ifndef EASYX_HEADER
#define EASYX_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include "common.h"

typedef struct
{
     Display *display;
     int screen;
     Window rootwindow;
     Visual *visual;
     int depth;
     GC gc;
     XFontStruct *font;
     Colormap colormap;
} EXInfoType;

typedef struct
{
     bool use_buffer;
     Window window;
     Pixmap buffer;
     long width, height;
}EXWindow;

extern EXInfoType EXInfo;
extern EXWindow* EXCW; // EX Current Window

EXWindow* EXStart(int x, int y, int width, int height);
EXWindow* EXCreateWindow(int x, int y, int width, int height);
void EXSetName(const char *name);
void EXClose();

void EXEnableBuffer();
void EXDisableBuffer();

void EXGetWindowSize(long &width, long &height);

void EXPixel(int x, int y);
void EXFlush();
void EXLine(int x0, int y0, int x, int y);
void EXCircle(int x, int y, int R);
void EXFillCircle(int x, int y, int R);
void EXFillRectangle(int x0, int y0, int w, int h);
void EXDrawRectangle(int x0, int y0, int w, int h);
void EXClear();

void EXSetLineWidth(int lw);
void EXFillPolygon(const List &X, const List &Y);

void EXGetEvent(XEvent *);
KeySym EXKey2Keysym (XEvent *event);
int EXKeyPressed();
char EXReadKey();

// Struct used to return the position where pointer button was pressed
typedef struct
{
     int x;
     int y;
     int button;
} EXPointerState;

int EXPointer();
int EXPointerPressed();
EXPointerState EXReadPointer();

extern const char Default_Font[];
void EXPrepareFont(const char *font_id);
void EXDrawString(int x, int y, const char *cadena);
long EXTextWidth(const char *cadena, long n);
long EXTextAscent();
long EXTextDescent();

XImage* EXGetImage(int x0, int y0, int xsize, int ysize);
void EXPutImage(XImage* imagen, int x0, int y0, int xsize, int ysize);

typedef unsigned long int palette;

void EXSetColor(palette colornum);
palette EXAllocNamedColor(const char *colorname);
palette EXAllocRGBColor(double red, double green, double blue);
palette* EXPalette(double R1, double G1, double B1,
		   double R2, double G2, double B2, int Ncolors);
	 
#endif
