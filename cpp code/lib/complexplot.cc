// Complex plot,
// Routines in order to show complex functions, with a color code
// 110912

#include"complexplot.h"

// returns a different "top" color, periodic-continuous in alpha
// alpha in [0,1] !!!
Vector Get_Color_From_Phase(double alpha)
{
     Vector R(3), G(3), B(3);
     R(1)=1.0; G(2)=1.0; B(3)=1.0;
     if (alpha<0.33)
     {
	  double lambda=3.0*alpha;
	  return (1.0-lambda)*R+lambda*G;
     }
     if (alpha<0.66)
     {
	  double lambda=3.0*(alpha-0.33);
	  return (1.0-lambda)*G+lambda*B;
     }
     double lambda=3.0*(alpha-0.66);
     return (1.0-lambda)*B+lambda*R;
}

long Find_Color_Index(cmplx z, long Nint, long Nphase, double Rmax)
{
     double R=abs(z);
     double argument=arg(z); //+M_PI;
     if (argument<0.0) argument+=2.0*M_PI; // should be reall a "clamp"
     if (R>Rmax) R=Rmax;
     long ii=(long)floor(R/Rmax*(Nint-1));
     long ip=(long)floor(argument/(2.0*M_PI)*(Nphase-1));
     // printf("ii: %ld, ip: %ld, Index: %ld\n",ii,ip,ip*Nint+ii);
     return ip*Nint+ii;
}

// white=true if we're using a white background; otherwise, it will be black
palette *Build_Palette_With_Phases(long Nc_int, long Nc_phase, bool white)
{
     palette *P=(palette*)malloc(Nc_int*Nc_phase*sizeof(palette));
//     printf("Palette built with %ld entries\n",Nc_int*Nc_phase);
     for (long ip=0;ip<Nc_phase;ip++)
     {
	  Vector PC=Get_Color_From_Phase((double)ip/(double)Nc_phase);
	  if (!white) // Black background
	  {
	       for (long ii=0;ii<Nc_int;ii++)
	       {
		    double fac=(double)ii/(double)(Nc_int-1);
		    Vector C=fac*PC;
		    P[ip*Nc_int+ii]=EXAllocRGBColor(C(1),C(2),C(3));
	       }
	  }
	  else 
	  {
	       Vector White(3); White.Set_Value(1.0);
	       Vector Dir=PC-White;
	       for (long ii=0;ii<Nc_int;ii++)
	       {
		    double fac=(double)ii/(double)(Nc_int-1);
		    Vector C=White+fac*Dir;
		    P[ip*Nc_int+ii]=EXAllocRGBColor(C(1),C(2),C(3));
	       }
	  }

                    // the top color for this phase
     }
     return P;
}

// dotsize: size of each little square
// white=true -> white background
// Rmax: saturation value
void Complex_Plot(const CMatrix &V, long Nc_int, long Nc_phase, long xsize, 
		  bool white, double saturation)
{
     // Find the maximum absolute value
     double Rmax=0.0;
     for (long i=1;i<=V.N1;i++)
	  for (long j=1;j<=V.N2;j++)
	  {
	       double f=abs(V(i,j));
	       if (f>Rmax) Rmax=f;
	  }

     Rmax*=saturation; // heuristic: by allowing a bit of saturation, it usually
     // looks nicer... but not necessarily so!
     
     // First, generate the palette
     palette *P=Build_Palette_With_Phases(Nc_int,Nc_phase,white);

     double dx=(double)xsize/(double)V.N1;
     long idx=MAX(1,(long)ceil(dx)+1);
     // Now, do the plotting
     for (long i=1;i<=V.N1;i++)
	  for (long j=1;j<=V.N2;j++)
	  {
	       cmplx w=V(i,j);
	       palette Color=P[Find_Color_Index(w,Nc_int,Nc_phase,Rmax)];
	       EXSetColor(Color);
	       EXFillRectangle((long)round((i-1)*dx),(long)round((j-1)*dx),
			       idx,idx);
	  }
     EXFlush();
}
