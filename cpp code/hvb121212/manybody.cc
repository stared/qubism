// MANYBODY - manybody calculations, complex numbers, high spin, no bmatrices
// Javier Rodriuez-aguna
// 100615

#include"manybody.h"

// return spin 1/2 matrix, 0, +1 or -1
CMatrix Spin_Op(int s) // for a single 1/2
{
     return Spin_Op(2,s);
}

// give multiplicity and component (0,+1,-1)
CMatrix Spin_Op(long mult, int s)
{
     double spin=mult/2.0-0.5;
     CMatrix S(mult);
     if (s==0)
     {
	  for (long i=1;i<=mult;i++)
	       S(i,i)=i-spin-1.0;
	  return S;
     }
     for (long i=1;i<mult;i++)
     {
	  double ms=i-spin-1.0;
	  double factor=sqrt(spin*(spin+1.0) - ms*(ms+1.0));
	  S(i+1,i)=factor;     
     }
     if (s==-1) S.Herm();
     return S;
}

CMatrix Sz_Op(long mult)
{
     return Spin_Op(mult,0);
}

CMatrix Sx_Op(long mult)
{
     return 0.5*(Spin_Op(mult,1)+Spin_Op(mult,-1));
}

CMatrix Sy_Op(long mult)
{
     return -0.5*M_I*(Spin_Op(mult,1)-Spin_Op(mult,-1));
}

CMatrix C_Op(int s) //  particle operator; 0: number, 1: creator, -1: annih
{
     CMatrix B(2);
     switch(s)
     {
     case 0: B(2,2)=1.0; break;
     case -1: B(1,2)=1.0; break;
     case 1: B(2,1)=1.0; break;
     }
     return B;
}

// for 2-states per site
CMatrix Site_Op(const CMatrix &B, long i, long N)
{
     if (i==1) return Tens_Prod_Unit(B,1<<(N-1));
     if (i==N) return Tens_Prod_Unit(1<<(N-1),B);
     return Tens_Prod_Unit(1<<(i-1),Tens_Prod_Unit(B,1<<(N-i)));
}

// for general number of states per site
// Ldim is a list with the dimensionality of each site
// B is a matrix with the correct dimensionality, of course!
CMatrix Site_Op(const CMatrix &B, long i, const List &Ldim)
{
     long Ntot=Ldim.Prod();
     long Np=Ntot/Ldim(i);
     long N=Ldim.N;
     if (i==1) return Tens_Prod_Unit(B,Np);
     if (i==N) return Tens_Prod_Unit(Np,B);
     long Np1=Ldim.Prod(1,i-1);
     long Np2=Ldim.Prod(i+1,N);
     return Tens_Prod_Unit(Np1,Tens_Prod_Unit(B,Np2));
}

// ITF hamiltonian with graph G, couplings J and field Gamma
// s is the spin multiplicity
CMatrix ITF_Ham(const Graph &G, const Vector &J, double Gamma, long s)
{
     long N=G.N; long Nl=G.Nl;
     CMatrix Ham;
     CMatrix Sz=Spin_Op(s,0);
     CMatrix Sx=Sx_Op(s);

     List Ldim(N); Ldim.Set_Value(s);
     for (long k=1;k<=Nl;k++)
     {
	  long s1, s2;
	  G.Get_Link_Sites(s1,s2,k);
	  CMatrix Sz1=Site_Op(Sz,s1,Ldim);
	  CMatrix Sz2=Site_Op(Sz,s2,Ldim);
	  CMatrix Szz=Sz1*Sz2;
	  Ham-=J(k)*Szz;
     }
     for (long i=1;i<=N;i++)
     {
	  CMatrix SX=Site_Op(Sx,i,Ldim);
	  Ham-=Gamma*SX;
     }
     return Ham;
}

// Heisenberg hamiltonian with graph G and couplings J
// s is the spin multiplicity
CMatrix Heisenberg_Ham(const Graph &G, const Vector &J, long s)
{
     long N=G.N; long Nl=G.Nl;
     CMatrix Ham;
     CMatrix Sz=Spin_Op(s,0);
     CMatrix Sp=Spin_Op(s,+1);
     CMatrix Sm=Spin_Op(s,-1);
     List Ldim(N); Ldim.Set_Value(s);
     for (long k=1;k<=Nl;k++)
     {
	  long s1, s2;
	  G.Get_Link_Sites(s1,s2,k);
	  CMatrix Sz1=Site_Op(Sz,s1,Ldim);
	  CMatrix Sz2=Site_Op(Sz,s2,Ldim);
	  CMatrix Szz=Sz1*Sz2;
	  Sz1.Destroy(); Sz2.Destroy();
	  CMatrix Sp1=Site_Op(Sp,s1,Ldim);
	  CMatrix Sm2=Site_Op(Sm,s2,Ldim);
	  CMatrix Spm=Sp1*Sm2;
	  Sp1.Destroy(); Sm2.Destroy();
	  Spm+=Herm(Spm);
	  Ham+=J(k)*( Szz + 0.5*( Spm ));
     }
     return Ham;
}

// // return the integer made up with the bits from L
// long Dec(long i, const List &L)
// {
//      // printf("Dec %s\n",B2string(i));
//      // L.Write();
//      long x=0;
//      for (long k=1;k<=L.N;k++)
// 	  if (get_bit(i,L(k)-1))
// 	       x=put_bit(x,k-1,1);
//      // printf("x: %ld\n",x);
//      return x;
// }

// Spin 1/2 system, trace Rho on the sites in list L
CMatrix Trace_On(const CMatrix &Rho, const List &Lsites)
{
     long Lt=(long)round(log((double)Rho.N1)/log(2.0)); // easier way?
     List Ldim(Lt); Ldim.Set_Value(2);
     return Trace_On(Rho,Lsites,Ldim);
}

// // Opposite of Select / Combination
// // expand list A to have L elements, Ai will go to position Bi.
// List Expand(const List &A, const List &B, long L)
// {
//      List R(L);
//      for (long i=1;i<=B.N;i++)
// 	  R(B(i))=A(i);
//      return R;
// }

// General spin system, trace Rho on the sites in list L1, they have
// dimensions in Ldim
CMatrix Trace_On(const CMatrix &Rho, const List &Ls1, const List &Ldim)
{
     long l=Ldim.N, l1=Ls1.N; // l2=l-l1;
     // Find the dimension of the resulting Hilbert space
     List Ldim1=Combination(Ls1,Ldim);
     List Ls=List_Interval(1,l);
     List Ls2=Substraction(Ls,Ls1);
     List Ldim2=Combination(Ls2,Ldim);
     long N1=Ldim1.Prod();
     long N=Ldim.Prod(); long N2=N/N1;

     CMatrix R(N1);
     for (long i1=1;i1<=N1;i1++)
	  for (long j1=1;j1<=N1;j1++)
	       for (long k2=1;k2<=N2;k2++)
	       {
		    List I1=To_List(i1-1,Ldim1);
		    List J1=To_List(j1-1,Ldim1);
		    List K2=To_List(k2-1,Ldim2);
		    List I(l), J(l);
		    for (long p=1;p<=l;p++)
			 if (Ls1.Find(p)) 
			 {
			      I(p)=I1(Ls1.Find(p));
			      J(p)=J1(Ls1.Find(p));
			 }
			 else 
			 {
			      I(p)=K2(Ls2.Find(p));
			      J(p)=K2(Ls2.Find(p));
			 }
		    long i=To_Number(I,Ldim)+1;
		    long j=To_Number(J,Ldim)+1;
		    R(i1,j1)+=Rho(i,j);
	       }
     return R;
}

double Entropy(const CMatrix &Rho)
{
     CMatrix Basis; Vector Eigen;
     Rho.Diagonalize(Basis,Eigen);
     double S=0.0;
     for (long i=1;i<=Eigen.N;i++)
	  S+=(Eigen(i)<1e-16 ? 0.0 : -Eigen(i)*log(Eigen(i)));
     return S; 
}

// Multiply If*..*If*c*I*..*I (tensor!)
// Not very efficient, but it works
CMatrix F_Site_Op(const CMatrix &c, long i, long L)
{
     static CMatrix If;
     If.Create(2); If.Zero(); If(1,1)=1.0; If(2,2)=-1.0;
     CMatrix Left(1); Left(1,1)=1.0; 
     for (long j=1;j<i;j++)
	  Left=Tens_Prod(Left,If);
     CMatrix C=Tens_Prod(Left,c);
     return Tens_Prod_Unit(C,1<<(L-i));
}
