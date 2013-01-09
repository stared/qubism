// MANYBODY - manybody calculations 
// Javier Rodriguez Laguna
// 100614-110509
#include"cmatrix.h"
#include"graph.h"

CMatrix Spin_Op(int ms); // for a single 1/2; 0: z, 1: p, -1: m.
CMatrix Spin_Op(long mult, int ms); // multiplicity and ms (-1,0,1)

CMatrix Sz_Op(long mult); // given the multiplicity, i.e.: 2 for spin 1/2
CMatrix Sx_Op(long mult);
CMatrix Sy_Op(long mult);

CMatrix C_Op(int s); // hard boson operator; 0: number, 1: creator, -1: annih

CMatrix Heisenberg_Ham(const Graph &G, const Vector &J, long s);
CMatrix ITF_Ham(const Graph &G, const Vector &J, double Gamma, long s);

CMatrix Site_Op(const CMatrix &B, long i, long N);

CMatrix Site_Op(const CMatrix &B, long i, const List &Ldim);

//CMatrix Get_DM(const CVector &V, long s, long L, const List &Ll);

CMatrix Trace_On(const CMatrix &Rho, const List &Lsites);

CMatrix Trace_On(const CMatrix &Rho, const List &Ls1, const List &Ldim);

double Entropy(const CMatrix &Rho);
		  
double Entropy(const CVector &V, long sector, long i, long L);

CMatrix F_Site_Op(const CMatrix &c, long i, long L);
