// Program to visualize a quantum manybody wf
// 100915-110503-110719-110912
#include"cmatrix.h"
#include"parsing.h"
#include"easyx.h"
#include"easyim.h"
#include"complexplot.h"
#include"wf.h"

typedef struct 
{
     long xsize; // window size
     long Ncint; // number of colors, for intensity
     long Ncphase; // number of colors, phase
     bool alt;
     bool probab, phase; 
     bool white;
     bool nolines;
     bool png, nokey, positive;
     long L, s; 
     double saturation;
}Plot_Parameters; 

////////////////////////////////////////////////
// THE WF-PLOT GENERATING PROCEDURES
////////////////////////////////////////////////

// Get the many-body vector Psi, representing a state with 2^L comps
// into a matrix, in the hyerarchical way.
CMatrix Standard_2D_Plot(const CVector &Psi, long L)
{
     List Ls(L/2); Ls.Set_Value(4);
     long Lt=1<<(L/2);
     CMatrix Im(Lt);
     for (long i=1;i<=Psi.N;i++)
     {
	  List D=To_List(i-1,Ls);; // number in base-4
	  long x=0, y=0, l=Lt/2;
	  for (long k=1;k<=D.N;k++)
	  {
	       if (D(k)==1) x+=l;
	       if (D(k)==2) y+=l;
	       if (D(k)==3) { x+=l; y+=l;}
	       l/=2;
	  }
	  Im(x+1,y+1)=Psi(i);
     }
     return Im;
} 

// Alternate version, "Piotr"-like
CMatrix Alt_2D_Plot(const CVector &Psi, long L)
{
     List Ls(L/2); Ls.Set_Value(4);
     long Lt=1<<(L/2);
     CMatrix Im(Lt);
     for (long i=1;i<=Psi.N;i++)
     {
	  List D=To_List(i-1,Ls);; // number in base-4
	  long x=0, y=0, l=Lt/2;
	  for (long k=1;k<=D.N;k++)
	  {
	       if (D(k)==1) x+=l;
	       if (D(k)==2) { x+=l; y+=l;}
	       if (D(k)==3) y+=l;
	       l/=2;
	  }
	  Im(x+1,y+1)=Psi(i);
     }
     return Im;
} 

CMatrix Spin1_2D_Plot(const CVector &P, long L)
{
     long ltot=(long)pow(3,L/2);
     long Ltot=ltot*ltot;
     CMatrix M(ltot);
     List Ls(L); Ls.Set_Value(3);
     for (long i=1;i<=Ltot;i++)
     {
	  List Li=To_List(i-1,Ls);
	  long x=0, y=0, l=ltot/3;
	  for (long k=1;k<=L;k+=2)
	  {
	       y+=Li(k)*l;
	       x+=Li(k+1)*l;
	       l/=3;
	  }
	  M(x+1,y+1)=P(i);
     }
     return M;
}


///////////////////////////////////////////////////////////////
// NOW, SHOW!
///////////////////////////////////////////////////////////////

// Draw a grid, dividing into pxp little squares
void Plot_Grid(long xsize,long p)
{
     double ps=(double)xsize/(double)p;
     for (long i=0;i<=p;i++)
     {
	  long ips=(long)round(i*ps);
	  if (ips<0) ips=0;
	  if (ips>=xsize) ips=xsize-1;
	  EXLine(ips,0,ips,xsize-1);
	  EXLine(0,ips,xsize-1,ips);
     }
}

void Plot_Lines(const Plot_Parameters &P)
{
     if (P.nolines) return;
     palette color1=(P.white ? EXAllocRGBColor(0.4,0.4,0.4) : 
		     EXAllocRGBColor(0.3,0.3,0.1));
     palette color2=(P.white ? EXAllocRGBColor(0.7,0.7,0.7) : 
		     EXAllocRGBColor(0.2,0.1,0.1));

     EXSetColor(color2);
     Plot_Grid(P.xsize,sqr(P.s));
     EXSetColor(color1);
     Plot_Grid(P.xsize,P.s);
     EXFlush();
}

// Make the phase of the maximum component real-positive
void Positivize(CVector &V)
{
     // First, locate the first non-zero comp
     long imax=V.Max_Index();
     cmplx value=V(imax);
     cmplx phase=value/abs(value);
     V*=conj(phase);
}

// get the effective wf when the first qubits are measured and the results are in B
CVector Effective_WF(const CVector &Psi, long L, long s, const List &B1)
{
     long L1=B1.N; // number of bits in given part
     long L2=L-L1; // number of bits in open part
     long Ntot2=(long)floor(pow(s,L2));
     CVector Psi2(Ntot2);
     List Ns(L); Ns.Set_Value(s);
     List Ns2(L2); Ns2.Set_Value(s);
     for (long i2=0;i2<Ntot2;i2++)
     {
	  List B2=To_List(i2,Ns2);
	  List B(B1); B.Append(B2);
	  long i=To_Number(B,Ns);
	  Psi2(i2+1)=Psi(i+1);
     }
     // double norm=Psi2.Norm();
     return Psi2;
}

void Shift(List &L, long p) // surely there is a faster way... to hvb????
{
     // printf("Shift by %ld\n",p);
     //L.Write();
     List L2(L);
     long N=L.N;
     for (long i=1;i<=N;i++)
     {
	  long ip=i+p;
	  if (ip>N) ip-=N;
	  if (ip<1) ip+=N;
	  L2(i)=L(ip);
     }
     L=L2;
     //L.Write();
}

void Shift(CVector &Psi, long L, long s, long p)
{
     printf("Shifting the wf by %ld\n",p);
     long Nt=Psi.N;
     List Ns(L); Ns.Set_Value(s);
     CVector Psi2(Psi);
     for (long i=1;i<=Nt;i++)
     {
	  List Li=To_List(i-1,Ns);
	  Shift(Li,p);
	  long ip=To_Number(Li,Ns)+1;
	  Psi(ip)=Psi2(i);
     }
}

// increase the even/odd qbits by p
void Inc_Bits(List &Bits, long s, bool even, long p)
{
     if (!Bits.N) return;
     List B2(Bits.N/2);
     for (long i=1;i<=B2.N;i++)
	  B2(i)=Bits(even?2*i:2*i-1);
     List Ns(B2.N); Ns.Set_Value(s);
     long I=To_Number(B2,Ns);
     I+=p;
     if (I<0) I+=(1<<B2.N);
     B2=To_List(I,Ns);
     for (long i=1;i<=B2.N;i++)
	  Bits(even?2*i:2*i-1)=B2(i);
}

// Show, and do the magnification magic
void Tweak_Show(CVector &Psi, const Plot_Parameters &P)
{
     CMatrix PsiM;
     if (P.s==3) PsiM=Spin1_2D_Plot(Psi,P.L);
     if (P.s==2) PsiM=(P.alt ? Alt_2D_Plot(Psi,P.L) : Standard_2D_Plot(Psi,P.L));
     Complex_Plot(PsiM,P.Ncint,P.Ncphase,P.xsize,P.white,P.saturation);
     Plot_Lines(P);
     if (P.nokey) return;

     List Bits; // list of "fixed bits", chosen by the user by pressing magnification keys
     bool done=false;
     do
     {
	  char q=EXReadKey();
	  if ((P.s==3) && (Bits.N<P.L))
	       switch(q)
	       {
	       case '1': Bits.Append(0); Bits.Append(0); break;
	       case '2': Bits.Append(0); Bits.Append(1); break;
	       case '3': Bits.Append(0); Bits.Append(2); break;
	       case '4': Bits.Append(1); Bits.Append(0); break;
	       case '5': Bits.Append(1); Bits.Append(1); break;
	       case '6': Bits.Append(1); Bits.Append(2); break;
	       case '7': Bits.Append(2); Bits.Append(0); break;
	       case '8': Bits.Append(2); Bits.Append(1); break;
	       case '9': Bits.Append(2); Bits.Append(2); break;
	       }
	  if ((P.s==2) && (Bits.N<P.L))
	  {
	       if (!P.alt) // quadrants are always numbered 1 2;3 4.
		    switch(q)
		    {
		    case '1': Bits.Append(0); Bits.Append(0); break;
		    case '2': Bits.Append(0); Bits.Append(1); break;
		    case '3': Bits.Append(1); Bits.Append(0); break;
		    case '4': Bits.Append(1); Bits.Append(1); break;
		    }
	       else
		    switch(q)
		    {
		    case '1': Bits.Append(0); Bits.Append(0); break;
		    case '2': Bits.Append(0); Bits.Append(1); break;
		    case '3': Bits.Append(1); Bits.Append(1); break;
		    case '4': Bits.Append(1); Bits.Append(0); break;
		    }
	  }
	  if (q=='0' && (Bits.N>0)) Bits.Part(1,Bits.N-2); // zoom out
	  if (q=='r' || q=='R') Shift(Psi,P.L,P.s,+1);
	  if (q=='l' || q=='L') Shift(Psi,P.L,P.s,-1);

	  if (q=='p' || q=='P') Inc_Bits(Bits,P.s,true,1);
	  if (q=='o' || q=='O') Inc_Bits(Bits,P.s,true,-1);
	  if (q=='q' || q=='Q') Inc_Bits(Bits,P.s,false,-1);
	  if (q=='a' || q=='A') Inc_Bits(Bits,P.s,false,1);

	  if (q=='z' || q=='Z') done=true;

	  Bits.Write();
	  CVector Psi_eff=Effective_WF(Psi,P.L,P.s,Bits);
	  long Leff=P.L-Bits.N;
	  
	  if (P.s==3) PsiM=Spin1_2D_Plot(Psi_eff,Leff);
	  if (P.s==2) PsiM=(P.alt ? Alt_2D_Plot(Psi_eff,Leff) : 
			   Standard_2D_Plot(Psi_eff,Leff));
	  
	  Complex_Plot(PsiM,P.Ncint,P.Ncphase,P.xsize,P.white,P.saturation);
	  Plot_Lines(P);
	  EXFlush();
     }while(!done);


}

int main(int argc, char *argv[])
{
     
     if (argc<2) 
     {
	  printf("Syntaxis: xqview filename [options]\n");
	  printf("Options:\n");
	  printf("alt: Alternative 2D plot\n");
	  printf("probab: Plot probability, not amplitude\n");
	  printf("phase: Plot only phases\n");
	  printf("png: Generate a PNG image\n");
	  printf("white: Use white background\n");
	  printf("nolines: Do not plot dividing lines\n");
	  printf("nokey: Do not wait for a keypress\n");
	  printf("SN: Saturation value (N in 1..10)\n");
	  printf("xsizeN: Window size in pixels (default N=500)\n");
	  printf("positive: Adapt global phase so that first non-zero comp is positive\n");
	  exit(1);
     }


     Plot_Parameters P;
     P.alt=P.probab=P.phase=P.white=P.nolines=P.png=P.nokey=P.positive=false;
     P.xsize=500;
     P.Ncint=40;
     P.Ncphase=40;
     P.saturation=0.7;

     CVector V;
     Read_WF(V,P.L,P.s,argv[1]);
  
     for (long i=1;i<=argc;i++)
     {
	  if (is_here(argv[i],"alt")) P.alt=true;
	  if (is_here(argv[i],"probab")) P.probab=true;
	  if (is_here(argv[i],"phase")) P.phase=true;
	  if (is_here(argv[i],"png")) P.png=true;
	  if (is_here(argv[i],"white")) P.white=true;
	  if (is_here(argv[i],"nokey")) P.nokey=true;
	  if (is_here(argv[i],"nolines")) P.nolines=true;
	  if (is_here(argv[i],"S"))
	  {
	       int p=get_int(argv[i],"S^",1);
	       P.saturation=p/10.0;
	  }
	  if (is_here(argv[i],"xsize"))
	       P.xsize=(long)get_int(argv[i],"e^",1);
	  if (is_here(argv[i],"positive")) P.positive=true;
     }

     if (P.positive) // find first non-zero component and make it real and positive
	  Positivize(V);

     EXStart(100,100,P.xsize,P.xsize);
     EXEnableBuffer();
     
     if (P.png) EI_Start();
     if (P.phase)
	  for (long i=1;i<=V.N;i++)
	  {
	       if (abs(V(i))>1e-13) // which number should we put here?
	       V(i)=V(i)/abs(V(i));
	  }
     if (P.probab)
	  for (long i=1;i<=V.N;i++)
	       V(i)=norm(V(i));

     Tweak_Show(V,P);
     
     if (P.png)
     {
	  char name[60];
	  sprintf(name,"%s.png",argv[1]);
	  EI_Capture(0,0,P.xsize,P.xsize);
	  EI_Save(name);
	  EI_Free();
     }
}

