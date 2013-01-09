////////////////////////////////////////////////////////
// hvb++ 1.0
// Copyleft: Javier Rodríguez Laguna
// 080725-111201

// Routines for easy connection and usage of X-Windows server. 
// JaviRL Jan-Dec99. 
// All you have to do is: include "easyx.h" in your header, and use the 
// appropriate compilation order adding -L/usr/X11R6/lib -lX11.  
// modifications: April 16, 2000; Feb 26, 2004.
// last modification: version 1.0, April 16, 2006

#ifndef EASYX
#define EASYX

#include "easyx.h"

EXInfoType EXInfo;
EXWindow* EXCW; // EX Current Window
const char Default_Font[]="-misc-fixed-medium-*-*-*-*-*-*-*-c-100-iso8859-*";

EXWindow* EXCreateWindow(int x, int y, int width, int height)
{
     unsigned long event_mask; 
     event_mask=StructureNotifyMask;

     unsigned long attr_mask=CWEventMask | CWBackPixel | 
	  CWBorderPixel | CWOverrideRedirect;
     
     XSetWindowAttributes attributes;
     attributes.event_mask=event_mask;
     attributes.border_pixel=WhitePixel(EXInfo.display,EXInfo.screen);
     attributes.background_pixel=BlackPixel(EXInfo.display,EXInfo.screen);
     attributes.override_redirect=false;
          
     Window window=XCreateWindow(EXInfo.display, EXInfo.rootwindow, 
				 x, y, width, height,
				 2, CopyFromParent, InputOutput,
				 EXInfo.visual, attr_mask, &attributes);
     //printf("created window: %d\n",window);
     // Size hints
     XSizeHints *size_hints = XAllocSizeHints();
     size_hints->x = x;            size_hints->y = y;
     size_hints->height = height;  size_hints->width = width;
     size_hints->min_height = height;
     size_hints->min_width = width;
     size_hints->flags = USPosition | USSize | PMinSize;
     size_hints->base_width = width;
     size_hints->base_height = height;
     size_hints->flags |= PBaseSize;
     XSetWMNormalHints(EXInfo.display, window, size_hints);
     XFree((char*)size_hints);

     // Class hints (name)
     XClassHint class_hints;
     class_hints.res_class=(char*)NULL;
     class_hints.res_name=(char*)NULL;
     XSetClassHint (EXInfo.display, window, &class_hints);
 
     // Window-manager hints
     XWMHints wm_hints;
     wm_hints.flags = InputHint | StateHint;
     wm_hints.initial_state = NormalState;
     wm_hints.input = true;
     XSetWMHints (EXInfo.display, window, &wm_hints);

     XMapRaised(EXInfo.display,window);

     XEvent evt;
     do
     {
	  XNextEvent( EXInfo.display , &evt );   // calls XFlush
	  //printf("event for window: %d\n",evt.xany.window);
     }while( evt.type != MapNotify || evt.xany.window!=window);

     usleep(1);

     event_mask=ExposureMask | PointerMotionMask | KeyPressMask | 
	  StructureNotifyMask | ButtonPressMask | ButtonMotionMask | 
	  ButtonReleaseMask;
     XSelectInput(EXInfo.display,window,event_mask);
     
//     XFlush(EXInfo.display);
     XSync(EXInfo.display, false);

     EXWindow* W=new EXWindow;
     W->window=window;
     W->use_buffer=false;
     W->width=width;
     W->height=height;	  
     EXCW=W;
     return W;
}

void EXSetName(const char *name)
{
     XStoreName(EXInfo.display, EXCW->window, name);
     XMapRaised(EXInfo.display, EXCW->window);
     XFlush(EXInfo.display);
}

// returns the index of the window
EXWindow* EXStart(int x, int y, int width, int height, const char *name)
{
     EXInfo.display=XOpenDisplay((char*)NULL);
     if (EXInfo.display==(Display*)NULL)
	  merror("I couldn't connect to Xserver");
     EXInfo.screen=DefaultScreen(EXInfo.display);
     EXInfo.rootwindow=RootWindow(EXInfo.display,EXInfo.screen);

     EXInfo.visual=DefaultVisual(EXInfo.display,EXInfo.screen);
     EXInfo.depth=DefaultDepth(EXInfo.display,EXInfo.screen);
     EXInfo.colormap=DefaultColormap(EXInfo.display,EXInfo.screen);

     EXWindow* W = EXCreateWindow(x, y, width, height);  
     XGCValues xgcvalues;
     xgcvalues.foreground=WhitePixel(EXInfo.display,EXInfo.screen);
     xgcvalues.background=BlackPixel(EXInfo.display,EXInfo.screen);
     EXInfo.gc=XCreateGC(EXInfo.display, W->window, 
			 (GCForeground | GCBackground), &xgcvalues);
     EXCW=W;
     return W;
}

// For backwards-compatibility reasons
EXWindow* EXStart(int x, int y, int width, int height)
{
     return EXStart(x,y,width,height,"");
}

void EXEnableBuffer()
{
     if (EXCW->use_buffer) return;
     XWindowAttributes wa;
     XGetWindowAttributes(EXInfo.display, EXCW->window, &wa);
     EXCW->buffer = XCreatePixmap(EXInfo.display, EXInfo.rootwindow,
				  wa.width, wa.height, wa.depth);
//     printf("Created buffer, with size: %ld\n",sizeof(EXCW->buffer));
     XSetForeground(EXInfo.display, EXInfo.gc, 
		    BlackPixelOfScreen(DefaultScreenOfDisplay(EXInfo.display)));
     XFillRectangle(EXInfo.display, EXCW->buffer, EXInfo.gc, 
		    0, 0, wa.width, wa.height);
     EXCW->use_buffer=true;

}

void EXDisableBuffer()
{
     XFreePixmap(EXInfo.display, EXCW->buffer);
     EXCW->use_buffer=false;
}

// returns the width of the WINDOW we're using
void EXGetWindowSize(long &width, long &height)
{
     Window rootw;
     int x, y; unsigned int w, h, border, depth;
     XGetGeometry(EXInfo.display,EXCW->window,&rootw,&x,&y,
		  &w,&h,&border,&depth);
     width=w;
     height=h;
}

void EXPixel(int x, int y)
{
     if (EXCW->use_buffer)
	  XDrawPoint(EXInfo.display,EXCW->buffer,EXInfo.gc,x,y);
     else
	  XDrawPoint(EXInfo.display,EXCW->window,EXInfo.gc,x,y);
}

void EXFlush()
{
     if (EXCW->use_buffer)
     {
	  XWindowAttributes wa;
	  XGetWindowAttributes(EXInfo.display, EXCW->window, &wa);
	  XCopyArea(EXInfo.display, EXCW->buffer, 
		    EXCW->window, EXInfo.gc,
		    0, 0, wa.width, wa.height, 0, 0);
     }
     // XMapRaised(EXInfo.display,EXCW->window);
     XFlush(EXInfo.display);
}

void EXClose()
{
     XCloseDisplay(EXInfo.display);
}

void EXLine(int x0, int y0, int x, int y)
{
     if (EXCW->use_buffer)
	  XDrawLine(EXInfo.display,EXCW->buffer,EXInfo.gc,x0,y0,x,y);
     else
	  XDrawLine(EXInfo.display,EXCW->window,EXInfo.gc,x0,y0,x,y);
}

void EXCircle(int x, int y, int R)
{
     if (EXCW->use_buffer)
	  XDrawArc(EXInfo.display,EXCW->buffer,EXInfo.gc,
		   x-R,y-R,2*R,2*R,0,64*360);
     else
	  XDrawArc(EXInfo.display,EXCW->window,EXInfo.gc,
		   x-R,y-R,2*R,2*R,0,64*360);
}

void EXFillCircle(int x, int y, int R)
{
     if (EXCW->use_buffer)
	  XFillArc(EXInfo.display,EXCW->buffer,EXInfo.gc,
		   x-R,y-R,2*R,2*R,0,64*360);
     else
	  XFillArc(EXInfo.display,EXCW->window,EXInfo.gc,
		   x-R,y-R,2*R,2*R,0,64*360);
}

void EXFillRectangle(int x0, int y0, int x1, int y1)
{
     if (EXCW->use_buffer)
	  XFillRectangle(EXInfo.display,EXCW->buffer,EXInfo.gc,
			 x0,y0,x1,y1);
     else
	  XFillRectangle(EXInfo.display,EXCW->window,EXInfo.gc,
			 x0,y0,x1,y1);
}

void EXDrawRectangle(int x0, int y0, int x1, int y1)
{
     if (EXCW->use_buffer)
	  XDrawRectangle(EXInfo.display,EXCW->buffer,EXInfo.gc,
			 x0,y0,x1,y1);
     else
	  XDrawRectangle(EXInfo.display,EXCW->window,EXInfo.gc,
			 x0,y0,x1,y1);
}

void EXClear()
{
     if (!EXCW->use_buffer)
	  XClearWindow(EXInfo.display,EXCW->window);
     else
     {
	  XWindowAttributes wa;
	  XGetWindowAttributes(EXInfo.display, EXCW->window, &wa);
	  XSetForeground(EXInfo.display, EXInfo.gc, 
			 BlackPixelOfScreen(DefaultScreenOfDisplay(EXInfo.display)));
	  XFillRectangle(EXInfo.display, EXCW->buffer, EXInfo.gc, 
			 0, 0, wa.width, wa.height);
	  XSetForeground(EXInfo.display, EXInfo.gc, 
			 WhitePixelOfScreen(DefaultScreenOfDisplay(EXInfo.display)));
     }
}

void EXSetLineWidth(int lw)
{
     XSetLineAttributes(EXInfo.display,EXInfo.gc,lw,LineSolid,CapNotLast,JoinMiter);
}

// Maybe the interface is not as efficient as it should...
void EXFillPolygon(const List &X, const List &Y)
{
     XPoint *P=(XPoint*)malloc(X.N*sizeof(XPoint));
     for (long i=0;i<X.N;i++)
     {
	  P[i].x=X(i+1);
	  P[i].y=Y(i+1);
     }
     if (EXCW->use_buffer)
	  XFillPolygon(EXInfo.display,EXCW->buffer,EXInfo.gc,
		       P,X.N,Convex,CoordModeOrigin);
     else
	  XFillPolygon(EXInfo.display,EXCW->window,EXInfo.gc,
		        P,X.N,Convex,CoordModeOrigin);
     free(P);
}

void EXGetEvent(XEvent *evento)
{
     XNextEvent(EXInfo.display,evento);
}

unsigned long EXAllocNamedColor(const char *colorname)
{
     XColor hardwarecolor, exactcolor;
     unsigned long color=0;
     int status;
     
     status=XAllocNamedColor(EXInfo.display,EXInfo.colormap,
			     colorname, &hardwarecolor,&exactcolor);
     if (status!=0) color=hardwarecolor.pixel;
     else printf("Error allocating color %s\n",colorname);
     return(color);
}

unsigned long EXAllocRGBColor(double red, double green, double blue)
{
     XColor search_color;
     unsigned long color=0;
     int status;
     
     search_color.red=(int)(65535*red);
     search_color.green=(int)(65535*green);
     search_color.blue=(int)(65535*blue);
     status=XAllocColor(EXInfo.display,EXInfo.colormap,&search_color);
     if (status==0) printf("Error allocating RGB color\n");
     else color=search_color.pixel;
     return(color);
}

void EXSetColor(unsigned long int colornum)
{
     XSetForeground(EXInfo.display,EXInfo.gc,colornum);
}

int EXKeyPressed()
{
     XEvent evento;
     int pressed;
     pressed=XCheckWindowEvent(EXInfo.display,EXCW->window,
			       KeyPressMask,&evento);
     XPutBackEvent(EXInfo.display,&evento);
     return pressed;    
}    

char EXReadKey()
{
     XEvent evento;
     KeySym tecla;
     
     do
     {
	  XNextEvent(EXInfo.display,&evento);
     }while(evento.type!=KeyPress);
     
     tecla=EXKey2Keysym(&evento);
     return (char)tecla;
}

int EXPointer()
{
     XEvent evento;
     int pressed;
     pressed=XCheckWindowEvent(EXInfo.display,EXCW->window, 
			       ButtonPressMask | 
			       ButtonReleaseMask | ButtonMotionMask, &evento);
     XPutBackEvent(EXInfo.display,&evento);
     return pressed;
}

int EXPointerPressed()
{
     XEvent evento;
     int pressed=XCheckWindowEvent(EXInfo.display,EXCW->window,
				   ButtonPressMask,&evento);
     XPutBackEvent(EXInfo.display, &evento);
     return pressed;
}
     	
EXPointerState EXReadPointer()
{
     XEvent evento;
     EXPointerState C;
//     unsigned int estado;
     static int button=0;
     
     do
     {	
	  XNextEvent(EXInfo.display,&evento);
     }while(evento.type!=MotionNotify && evento.type!=ButtonPress && 
	    evento.type!=ButtonRelease && evento.xany.window!=EXCW->window);
     if (evento.type==MotionNotify)
     {
	  C.x=evento.xmotion.x;
	  C.y=evento.xmotion.y;
	  C.button=button;
     }
     if (evento.type==ButtonPress)
     {
	  C.x=evento.xbutton.x;
	  C.y=evento.xbutton.y;
	  C.button=evento.xbutton.button;
	  button=C.button;
     }
     if (evento.type==ButtonRelease)
     {
	  C.x=evento.xbutton.x;
	  C.y=evento.xbutton.y;
	  button=0;
	  C.button=0;
     }
     return C;
}

void EXPrepareFont(const char *font_id)
{
//     XFontStruct *font_struct;
     EXInfo.font=XLoadQueryFont(EXInfo.display,font_id);     
     if (EXInfo.font == (XFontStruct *)NULL)
     {
	  /* if this fails, go to nofont, which is the default */
	  fprintf(stderr,"Error loading font [%s]\n",font_id);
	  EXInfo.font=XLoadQueryFont(EXInfo.display,"fixed");
	  if (EXInfo.font ==(XFontStruct *)NULL)
	  {
	       fprintf(stderr,"Error loading font fixed\n");
	       XCloseDisplay(EXInfo.display);
	       exit(1);
	  }
     }
//     EXInfo.font=EXLoadFont(EXInfo.display,font_id,"fixed");
     XSetFont(EXInfo.display,EXInfo.gc,EXInfo.font->fid);
}

void EXDrawString(int x, int y, const char *cadena)
{
     if (EXCW->use_buffer)
	  XDrawString(EXInfo.display,EXCW->buffer,EXInfo.gc,
		      x,y,cadena,strlen(cadena));
     else
	  XDrawString(EXInfo.display,EXCW->window,EXInfo.gc,
		      x,y,cadena,strlen(cadena));
}

long EXTextAscent()
{
     return EXInfo.font->ascent;
}

long EXTextDescent()
{
     return EXInfo.font->descent;
}

long EXTextWidth(const char *cadena, long n)
{
     long nr=n;
     if (!n) nr=strlen(cadena);
     return XTextWidth(EXInfo.font, cadena, nr);
}

XImage* EXGetImage(int x0, int y0, int xsize, int ysize)
{
     if (EXCW->use_buffer)
	  return XGetImage(EXInfo.display,EXCW->buffer,
			   x0,y0,xsize,ysize,AllPlanes,XYPixmap);
     else
	  return XGetImage(EXInfo.display,EXCW->window,
			   x0,y0,xsize,ysize,AllPlanes,XYPixmap);
}

void EXPutImage(XImage* imagen, int x0, int y0, int xsize, int ysize)
{
     if (EXCW->use_buffer)
	  XPutImage(EXInfo.display,EXCW->buffer,EXInfo.gc,
		    imagen,0,0,x0,y0,xsize,ysize);
     else
	  XPutImage(EXInfo.display,EXCW->window,EXInfo.gc,
		    imagen,0,0,x0,y0,xsize,ysize);
}

palette* EXPalette(double R1, double G1, double B1,
		   double R2, double G2, double B2, int Ncolors)
{
     unsigned long int *palette;
     int i;
     double s, rojo, verde, azul;
     palette=(unsigned long int*)malloc((Ncolors+1)*sizeof(long int));
     palette[0]=Ncolors;
     for (i=1;i<=Ncolors;i++)
     {
	 s=(double)i/(double)Ncolors;
	 rojo=R1+s*(R2-R1);
	 verde=G1+s*(G2-G1);
	 azul=B1+s*(B2-B1);
	 palette[i]=EXAllocRGBColor(rojo,verde,azul);
     }
     return palette;
}     

/////////////////////

KeySym EXKey2Keysym (XEvent *event)
{
     XComposeStatus compose;
     KeySym keysym;
     XKeyEvent *keyevent;
     char cad[20];
     
     keyevent=(XKeyEvent *) event;
     XLookupString(keyevent,cad,19,&keysym,&compose);
     return(keysym);
}



// XFontStruct *EXLoadFont(const char font_name[], 
// 		      const char nofont_name[])
// {
//      XFontStruct *font_struct;
//      font_struct=XLoadQueryFont(display,font_name);                
//      if (font_struct == (XFontStruct *)NULL)
//      {
// 	  /* if this fails, go to nofont, which is the default */
// 	  fprintf(stderr,"Error loading font [%s]\n",font_name);
// 	  font_struct=XLoadQueryFont(display,nofont_name);
// 	  if (font_struct ==(XFontStruct *)NULL)
// 	  {
// 	       fprintf(stderr,"Error loading font [%s]\n", nofont_name);
// 	       XCloseDisplay(display);
//                     exit(1);
// 	  }
//      }
//      return (font_struct);
// }

	 
#endif
