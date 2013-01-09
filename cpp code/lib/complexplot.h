// Complex plot,
// Routines in order to show complex functions, with a color code
// 110912

#include"cmatrix.h"
#include"easyx.h"

// returns a different "top" color, periodic-continuous in alpha
// alpha in [0,1] !!!
Vector Get_Color_From_Phase(double alpha);

long Find_Color_Index(cmplx z, long Nint, long Nphase, double Rmax);

palette *Build_Palette_With_Phases(long Nc_int, long Nc_phase, bool white);

void Complex_Plot(const CMatrix &V, long Nc_int, long Nc_phase, long xsize, bool white, double saturation);
