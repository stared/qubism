// Routines to deal with wavefunctions in a consistent way
// 110914

#include"cmatrix.h"
#include"parsing.h"

void Write_WF(const CVector &Psi, long L, long s, const char *filename);
void Read_WF(CVector &Psi, long &L, long &s, const char *filename);
