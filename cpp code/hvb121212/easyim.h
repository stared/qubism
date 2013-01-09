// 100914, new version of easyim, uses Imlib2
// Older version, 050511, with Imlib11

// Works with a "current image"

#ifndef EASYIM_HEADER
#define EASYIM_HEADER

#include"easyx.h"
#include"Imlib2.h"

typedef Imlib_Image EI_Image;

void EI_Start();
EI_Image EI_Load(const char *name); // also sets current image
int EI_Get_Width();
int EI_Get_Height();
void EI_Render(int x, int y);
void EI_Render_Scaled(int x, int y, int w, int h);
void EI_Free();
EI_Image EI_Capture(int x, int y, int w, int h); // also sets current image
void EI_Save(const char *name);

#endif
