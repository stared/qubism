////////////////////////////////////////////////////////
// hvb++ 1.0, CMATRIX
// Copyleft: Javier Rodriguez Laguna
// 080725-121211

#include "cmatrix.h"

#ifndef CMATRIX3
#define CMATRIX3

///////////////////////////////////////////////////////////////////////////
//  Real Implementations
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////
/// CLASS CVECTOR

CVector::CVector(): N(0), D(NULL) {} 

CVector::CVector(long n)
{
      Start(n);
      Zero();
}

CVector::CVector(const CVector &V)
{
     Start(V.N);
     if (N) memcpy(D,V.D,(N+1)*sizeof(cmplx));
}

CVector::CVector(const Vector &V)
{
     Start(V.N);
     if (N) 
	  for (long i=1;i<=N;i++)
	       D[i]=V(i);
}

CVector::CVector(cmplx *data, long n)
{
     Start(n);
     memcpy(D,data,(n+1)*sizeof(cmplx));
}

CVector::~CVector() { Destroy(); }

// Here we're sure that the CVector exists beforehand,
// So it is previously destroyed
void CVector::Create(long n)
{
     Destroy();
     Start(n);
}

void CVector::Create()
{
     Create(0);
}

// Use when you are sure that the CVector is not created
// Useful when initializing an array of Matrices
void CVector::Start(long n)
{
     N=n;
     if (!N) { D=NULL; return;}
     D=(cmplx*)malloc((n+1)*sizeof(cmplx));    
     
#ifdef DEBUG 
     if (!D) merror("Error allocating vector."); 
#endif
#ifdef MEM_CONTROL
     total_memory+=N;
     if (total_memory>max_total_memory) 
	  max_total_memory=total_memory;
     num_vec++;
#endif
}

void CVector::Load(cmplx *d, long n)
{
     Destroy();
     N=n;
     D=d;
#ifdef MEM_CONTROL
     num_vec++;
     total_memory+=N;
     if (total_memory>max_total_memory) 
	  max_total_memory=total_memory;
#endif
}

void CVector::Load_Copy(cmplx *d1, long n)
{
     cmplx *d2=(cmplx*)malloc((n+1)*sizeof(cmplx));
     memcpy(d2,d1,(n+1)*sizeof(cmplx));
     Load(d2,n);
}

// Acts as if V goes into our CVector
void CVector::Transfer(CVector &V)
{
     Destroy();
     D=V.D;
     N=V.N;
     V.D=(cmplx*)NULL;
     V.N=0; 
}

void CVector::Destroy()
{
     if (N) 
     { 	
	  free(D);
#ifdef MEM_CONTROL
	  total_memory-=N;
	  num_vec--;
#endif
	  N=0;
	  D=NULL;
     }
}

void CVector::Zero()
{
     if (N) memset(D,0,(N+1)*sizeof(cmplx));
}

bool CVector::Is_Zero(double tolerance)
{
     for (long i=1;i<=N;i++)
	  if (abs(D[i])>tolerance) return false;
     return true;
}

void CVector::Set_Value(cmplx x)
{
     if (N) for (long i=1;i<=N;i++) D[i]=x;
}

double CVector::Norm() const
{
     long ix=1, n=N;
     cmplx norm=zdotc_(&n,D+1,&ix,D+1,&ix);
     return sqrt(real(norm));
}

void CVector::Write(FILE *fich) const
{
     for (long i=1;i<=N;i++)
	  fprintf(fich,"(%16.12g,%16.12g)",real(D[i]),imag(D[i]));
     fprintf(fich,"\n\n");
}

void CVector::Write() const
{
     Write(stdout);
}

void CVector::Write_Col(FILE *fich) const
{
     for (long i=1;i<=N;i++)
	  fprintf(fich," %16.12g %16.12g\n",real(D[i]),imag(D[i]));
     fprintf(fich,"\n");
}

void CVector::Write_Col() const
{
     Write_Col(stdout);
}

void CVector::Write(int prec, FILE *f) const
{
     char form[10];
     sprintf(form,"(%%%d.%dg,%%%d.%dg) ",prec+4,prec,prec+4,prec);
     for (long i=1;i<=N;i++)
	  fprintf(f,form,real(D[i]),imag(D[i]));
     fprintf(f,"\n\n");
}

void CVector::Write(int prec) const
{
     Write(prec,stdout);
}

void CVector::Append(const cmplx x) 
{
     if (!N) Create(1);
     else
     {
	  N++;
	  D=(cmplx*)realloc(D,(N+1)*sizeof(cmplx));
#ifdef MEM_CONTROL
	  total_memory++;
	  if (total_memory>max_total_memory)
	       max_total_memory=total_memory;
#endif
     }
     D[N]=x;     
}

void CVector::Append(const CVector &V)
{
     long nold=N;
     if (!N) Create (V.N);
     else
     {
          N+=V.N;
          D=(cmplx*)realloc(D,(N+1)*sizeof(cmplx));
#ifdef MEM_CONTROL
          total_memory+=V.N;
          if (total_memory>max_total_memory)
               max_total_memory=total_memory;
#endif
     }
     for (long i=1;i<=V.N;i++)
          D[nold+i]=V(i);
}

// insert x in position i
void CVector::Insert(const cmplx x, long i)
{
     N++;
     D=(cmplx*)realloc(D,(N+1)*sizeof(cmplx));
#ifdef MEM_CONTROL
     total_memory++;
     if (total_memory>max_total_memory)
	  max_total_memory=total_memory;
#endif
     memmove(D+i+1,D+i,(N-i)*sizeof(cmplx));
     D[i]=x;
}

void CVector::Remove(long i)
{
     memmove(D+i,D+i+1,(N-i-1)*sizeof(cmplx));
     N--;
     D=(cmplx*)realloc(D,(N+1)*sizeof(cmplx));
#ifdef MEM_CONTROL
     total_memory--;
#endif
}

// Absolute lowest value
cmplx CVector::Min() const
{
     double E0=abs(D[1]);
     for (long i=2;i<=N;i++)
	  if (abs(D[i])<E0) E0=abs(D[i]);
     return E0;
}

// Get the k-th lowest value... not very fast, but effective
cmplx CVector::Min(long k) const
{
     CVector Acum;
     Copy(Acum,*this);
     Acum.Order();
     return Acum(k);
}

long CVector::Min_Index() const
{
     double E0=abs(D[1]);
     long imin=1;
     for (long i=2;i<=N;i++)
	  if (abs(D[i])<E0) { E0=abs(D[i]); imin=i; }
     return imin;
}

cmplx CVector::Max() const
{
     double E0=abs(D[1]);
     for (long i=2;i<=N;i++)
	  if (abs(D[i])>E0) E0=abs(D[i]);
     return E0;
}

cmplx CVector::Max(long k) const
{
     CVector Acum;
     Copy(Acum,*this);
     Acum.Order(-1);
     return Acum(k);
}

long CVector::Max_Index() const
{
     double E0=abs(D[1]);
     long imax=1;
     for (long i=2;i<=N;i++)
	  if (abs(D[i])>E0) { E0=abs(D[i]); imax=i; }
     return imax;
}

cmplx CVector::Sum() const
{
     if (!N) return (cmplx)0.0;
     cmplx sum=D[1];
     for (long i=2;i<=N;i++)
	  sum+=D[i];
     return sum;
}

cmplx CVector::Sum(long i1, long i2) const
{
     cmplx sum=D[i1];
     for (long i=i1+1;i<=i2;i++)
	  sum+=D[i];
     return sum;
}

cmplx CVector::Average() const
{
     return Sum()/(double)N;
}

cmplx CVector::Variance() const
{
     if (N==1) return 0.0;
     cmplx sumsq=0.0;
     for (long i=1;i<=N;i++)
          sumsq+=D[i]*D[i];
     sumsq/=(double)N;
     cmplx aver=Average();
     return sumsq-aver*aver;
}

cmplx CVector::Deviation() const
{
     return sqrt(Variance());
}

int complex_order_12(const void* x, const void* y) // ascending order
{
     const double dx=abs(*((const cmplx *)x));
     const double dy=abs(*((const cmplx *)y));
     return (dx > dy) - (dx < dy);
}

int complex_order_21(const void* x, const void* y) // descending order
{
     const double dx=abs(*(const cmplx *)x);
     const double dy=abs(*(const cmplx *)y);
     return (dx < dy) - (dx > dy);
}

void CVector::Order(int p)
{     
     if (p<0)
	  qsort(D+1,N,sizeof(cmplx),complex_order_21);
     else
	  qsort(D+1,N,sizeof(cmplx),complex_order_12);
}

cmplx CVector::operator() (long n) const
{
#ifdef DEBUG
     if (n<0 || n>N) merror("Error getting Cvector comp."); 
#endif
    return D[n];
}

cmplx& CVector::operator() (long n)
{
#ifdef DEBUG
    if (n<0 || n>N) merror("Error putting vector comp."); 
#endif
    return D[n];
}

cmplx CVector::Last() const
{
     return D[N];
}

CVector& CVector::operator=(const CVector& W)
{
     if (!W.N) { Destroy(); return (*this); }
     if (this==&W) return *this;
     Copy(*this,W);
     return (*this);
}

CVector& CVector::operator=(const Vector& W)
{
     if (!W.N) { Destroy(); return (*this); }
     Copy(*this,W);
     return (*this);
}


void CVector::operator+=(const CVector& W)
{
     if (!N) { Create(W.N); Zero(); }
#ifdef DEBUG
      if (N!=W.N) merror("Incompatible sizes in vector +=");
#endif
      cmplx alpha=1.0;
      long ix=1;
      zaxpy_(&N,&alpha,W.D+1,&ix,D+1,&ix);
}

void CVector::operator-=(const CVector& W)
{
     if (!N) { Create(W.N); Zero(); }
#ifdef DEBUG
     if (N!=W.N) merror("Incompatible sizes in vector -=");
#endif
     cmplx alpha=-1.0;
     long ix=1;
     zaxpy_(&N,&alpha,W.D+1,&ix,D+1,&ix);
}

void CVector::operator*=(const cmplx x)
{
#ifdef DEBUG
    if (!N) return;
#endif
    long ix=1;
    cmplx cosa=x;
    zscal_(&N,&cosa,D+1,&ix);
}

void CVector::operator*=(const double x)
{
     *this*=(cmplx)x;
}


void CVector::operator/=(const cmplx x)
{
#ifdef DEBUG
     if (!N) return;
#endif
     long ix=1;
     cmplx cosa=1.0/x;
     zscal_(&N,&cosa,D+1,&ix);
}

// return false if normalization was not possible!
bool CVector::Normalize()
{
     cmplx norm=Norm();
     if (norm==0.0) return false;
     (*this)/=norm;
     return true;
}

CVector Normalize(const CVector &V)
{
     CVector R(V);
     R.Normalize();
     return R;
}

void CVector::Part(long n1, long n2)
{
#ifdef DEBUG
     if (n1<1 || n2>N) 
	  merror("The part can't be larger than the whole!");
#endif
     CVector R(n2-n1+1);
     memcpy(R.D+1,D+n1,(n2-n1+1)*sizeof(cmplx));
     Transfer(R);
}

CVector Part(const CVector &V, long n1, long n2)
{
     CVector R(V);
     R.Part(n1,n2);
     return R;
}

void CVector::Reverse()
{
     for (long i=1;i<=N/2;i++)
     {
	  cmplx acum=D[i];
	  D[i]=D[N+1-i];
	  D[N+1-i]=acum;
     }
}

CVector Reverse(const CVector &V)
{
     CVector R(V);
     R.Reverse();
     return R;
}

void CVector::Insert(const CVector& V, long n)
{
#ifdef DEBUG
     if (n+V.N-1>N)
	  merror("Inserting a vector which is too big.");
#endif
     memcpy(D+n,V.D+1,V.N*sizeof(cmplx));
}

CVector Insert(const CVector &V, const CVector &W, long n)
{
     CVector R(V);
     R.Insert(W,n);
     return R;
}

void CVector::Sqr()
{
     for (long i=1;i<=N;i++)
     {
	  cmplx x=D[i];
	  D[i]=x*x;
     }
}

CVector Sqr(const CVector &V) 
{
     CVector R(V);
     R.Sqr();
     return R;
}

void CVector::Real()
{
     for (long i=1;i<=N;i++)
	  D[i]=real(D[i]);
}

CVector Real(const CVector &V)
{
     CVector R(V);
     R.Real();
     return R;
}

void CVector::Imag()
{
     for (long i=1;i<=N;i++)
	  D[i]=imag(D[i]);
}

CVector Imag(const CVector &V)
{
     CVector R(V);
     R.Imag();
     return R;
}

void CVector::Abs()
{
     for (long i=1;i<=N;i++)
	  D[i]=abs(D[i]);
}

CVector Abs(const CVector &V)
{
     CVector R(V);
     R.Abs();
     return R;
}

void CVector::Conj()
{
     for (long i=1;i<=N;i++)
	  D[i]=conj(D[i]);
}

CVector Conj(const CVector &V)
{
     CVector R(V);
     R.Conj();
     return R;
}

Vector To_Real(const CVector &V)
{
     Vector R(V.N);
     for (long i=1;i<=V.N;i++)
	  R(i)=(double)real(V(i));
     return R;
}


// CVector CVector::operator+(const CVector& V) const
// {
//      CVector R(V);
//      R+=(*this);
//      return R;
// }

// CVector CVector::operator-(const CVector& V) const
// {
//      CVector R(V);
//      R-=(*this);
//      return R;
// }

// CVector CVector::operator-() const
// {
//      CVector R(*this);
//      R*=-1.0;
//      return R;
// }

int CVector::Save_Binary(FILE *fich) const
{
     int ausgang;
     ausgang=fwrite(&N,sizeof(long),1,fich);
     if (ausgang!=1) return 0;
     if (!N) return 1;
     ausgang=fwrite(D,sizeof(cmplx),N+1,fich);
     if (ausgang!=N+1) return 0;
     return 1;
}

int CVector::Save_Binary(const char *name) const
{
     FILE *fich=fopen(name,"wb");
     if (!fich) return 0;
     int status=Save_Binary(fich);
     fclose(fich);
     return status;
}

int CVector::Load_Binary(FILE *fich)
{
     int ausgang;
     ausgang=fwrite(&N,sizeof(long),1,fich);
     if (ausgang!=1) return 0;
     if (!N) { D=(cmplx*)NULL; return 1; }
     Create(N);
     ausgang=fwrite(D,sizeof(cmplx),N+1,fich);
     if (ausgang!=N+1) return 0;
     return 1;
}

int CVector::Load_Binary(const char *name)
{
     FILE *fich=fopen(name,"rb");
     if (!fich) return 0;
     int status=Load_Binary(fich);
     fclose(fich);
     return status;
}

int CVector::Save(const char *name) const
{
     FILE *fich=fopen(name,"wt");
     if (!fich) return 0;
     int ausgang=Save(fich);
     fclose(fich);
     return ausgang;
}

int CVector::Save(FILE *fich) const
{
     for (long i=1;i<=N;i++)
	  fprintf(fich,"(%14.10g,%14.10g)\n",real(D[i]),imag(D[i]));
     return 1;
}

int CVector::Load(const char *name) 
{
     FILE *fich=fopen(name,"rt");
     if (!fich) return 0;
     int ausgang=Load(fich);
     fclose(fich);
     return ausgang;
}

int CVector::Load(FILE *fich)
{
     Create(0);
     float x,y;
     while(!feof(fich))
     {
	  if (!fscanf(fich,"%g %g",&x,&y))
	       merror("Error in CVector::Load");
	  cmplx z(x,y);
	  Append(z);
     }
     return 1;
}

///////////////////////////////////////////////////////
// EXTERNAL VECTOR FUNCTIONS

void Copy(CVector& B, const CVector& A)
{
     if (B.N) B.Destroy(); 
     if (!A.N) return;
     B.Start(A.N);
     memcpy(B.D,A.D,(A.N+1)*sizeof(cmplx));
}

void Copy(CVector& B, const Vector& A)
{
     if (B.N) B.Destroy(); 
     if (!A.N) return;
     B.Start(A.N);
     for (long i=1;i<=A.N;i++)
	  B(i)=A(i);
}

CVector operator+(const CVector &V, const CVector &W)
{
#ifdef DEBUG
     if (V.N!=W.N) merror("Inconsistent dimensions in vector +\n");
#endif
     CVector R(V.N);
     for (long i=1;i<=V.N;i++)
	  R(i)=V(i)+W(i);
     return R;
}

CVector operator-(const CVector &V, const CVector &W)
{
#ifdef DEBUG
     if (V.N!=W.N) merror("Inconsistent dimensions in vector +\n");
#endif
     CVector R(V.N);
     for (long i=1;i<=V.N;i++)
	  R(i)=V(i)-W(i);
     return R;
}

CVector operator-(const CVector &V)
{
     CVector R(V.N);
     return (-1.0)*R;
}

CVector operator*(cmplx x, const CVector& V)
{
     CVector R(V);
     R*=x;
     return R;
}

CVector operator*(double x, const CVector& V)
{
     CVector R(V);
     R*=x;
     return R;
}


CVector operator*(const CVector& V, cmplx x)
{
     CVector R(V);
     R*=x;
     return R;
}

CVector operator*(const CVector& V, double x)
{
     CVector R(V);
     R*=x;
     return R;
}


CVector operator/(const CVector& V, cmplx x)
{
     CVector R(V);
     R/=x;
     return R;
}

cmplx Dot(const CVector& V1, const CVector& V2)
{
#ifdef DEBUG
     if (V1.N!=V2.N) merror("Dot product of vectors of diferent dim.");
#endif
     long ix=1, n=V1.N;
     return zdotc_(&n,V1.D+1,&ix,V2.D+1,&ix);
}

CVector Cross(const CVector& V, const CVector& W)
{
#ifdef DEBUG
     if (V.N!=3 || W.N!=3) merror ("Cross product: vector dim must be 3.");
#endif
     CVector R(3);
     R(1)=V(2)*W(3)-V(3)*W(2);
     R(2)=V(3)*W(1)-V(1)*W(3);
     R(3)=V(1)*W(2)-V(2)*W(1);
     return R;
}

CVector Tens_Prod(const CVector &V, const CVector &W)
{
     long N=V.N*W.N;
     CVector R(N);
     long k=1;
     for (long i=1;i<=V.N;i++)
	  for (long j=1;j<=W.N;j++)
	  {
	       R.D[k]=V.D[i]*W.D[j];
	       k++;
	  }
     return R;
}

void Tens_Prod(CVector &R, const CVector &V, const CVector &W)
{
     long N=V.N*W.N;
     R.Create(N);
     long k=1;
     for (long i=1;i<=V.N;i++)
	  for (long j=1;j<=W.N;j++)
	  {
	       R.D[k]=V.D[i]*W.D[j];
	       k++;
	  }      
}

CVector Elem_Mult(const CVector &A, const CVector &B)
{
     CVector V(A.N);
     Elem_Mult(V,A,B);
     return V;
}

void Elem_Mult(CVector &V, const CVector &A, const CVector &B)
{
#ifdef DEBUG
     if (A.N!=B.N) merror("Different dimensions in Elem_Mult\n");
#endif
     long N=A.N;
     if (V.N!=N) V.Create(N);
     for (long i=1;i<=N;i++)
	  V(i)=A(i)*B(i);
}

// R <- alpha * V + R
void Zaxpy(CVector &R, const CVector &V, cmplx alpha)
{
     long ix=1, n=R.N;
     zaxpy_(&n,&alpha,V.D+1,&ix,R.D+1,&ix);     
}

///////////////////////////////////////////////////////////////
// MATRICIAL METHODS
///////////////////////////////////////////////////////////////

CMatrix::CMatrix(long n1, long n2) : N1(n1), N2(n2)
{
     N1=N2=0;
     Start(n1,n2);
     Zero();
}

CMatrix::CMatrix(const CMatrix & M)
{
     N1=N2=0; D=(cmplx*)NULL;
     Copy(*this,M);
}

CMatrix::CMatrix(const Matrix &M)
{
     N1=N2=0; D=(cmplx*)NULL;
     Copy(*this,M);
}

CMatrix::CMatrix()
{
     N1=N2=0;
     D=NULL;
}

CMatrix::~CMatrix()
{
     Destroy();
}

void CMatrix::Create(long n1, long n2)
{
     Destroy();
     Start(n1, n2);
}

void CMatrix::Create()
{
     Create(0,0);
}

void CMatrix::Start(long n1, long n2) // n2=0
{
     N1=n1; N2=n2;
     if (!N2) N2=N1;
     if (!N1)
     {
	  D=NULL;
	  return;
     }
     D=(cmplx*)malloc((N1*N2+1)*sizeof(cmplx));
#ifdef DEBUG
     if (!D) merror ("Error allocating matrix.");
#endif
#ifdef MEM_CONTROL
     total_memory+=N1*N2;
     if (total_memory>max_total_memory)
	  max_total_memory=total_memory;
     num_mat++;
#endif
}

// CAUTION: the data must be stored in columns
void CMatrix::Load(cmplx* d, long n1, long n2)
{
     Destroy();
     N1=n1;
     N2=n2; if (!N2) N2=N1;
     D=d;
#ifdef MEM_CONTROL
     total_memory+=N1*N2;
     if (total_memory>max_total_memory)
	  max_total_memory=total_memory;
     num_mat++;
#endif
}

// CAUTION: the data must be stored in columns
void CMatrix::Load_Copy(cmplx *d1, long n1, long n2)
{
     if (!n2) n2=n1;
     cmplx *d2=(cmplx*)malloc((n1*n2+1)*sizeof(cmplx));
     memcpy(d2,d1,(n1*n2+1)*sizeof(cmplx));
     Load(d2,n1,n2);
}

void CMatrix::Transfer(CMatrix &M)
{
     Destroy();
     D=M.D;
     N1=M.N1;
     N2=M.N2;
     M.N1=M.N2=0;
     M.D=(cmplx*)NULL;
}

void CMatrix::Destroy()
{
     if (N1) free(D);
     D=NULL;
#ifdef MEM_CONTROL
     total_memory-=N1*N2;
     if (N1*N2) num_mat--;
#endif
     N1=N2=0;
}

void CMatrix::Zero()
{
     memset(D+1,0,N1*N2*sizeof(cmplx));
}

bool CMatrix::Is_Zero(double tolerance)
{
     for (long i=1;i<=N1*N2;i++)
	  if (abs(D[i])>tolerance) return false;
     return true;
}

void CMatrix::Unit()
{
     Zero();
     long n=MIN(N1,N2);
     for (long i=1;i<=n;i++)
	  Elem(i,i)=1.0;
}

void CMatrix::Set_Value(cmplx x)
{
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=x;
}

void CMatrix::Resize(long n1, long n2) // m=0
{
     if (!n2) n2=n1;
     if (n1==N1 && n2==N2) return; 
     long nr1=MIN(N1,n1), nr2=MIN(N2,n2);
     long oldN1=N1, oldN2=N2;
     
     cmplx *D2=(cmplx*)malloc((N1*N2+1)*sizeof(cmplx));
     memcpy(D2+1,D+1,N1*N2*sizeof(cmplx));
     
     Create(n1,n2);
     if (n1>oldN1 || n2>oldN2)
	  Zero(); // In case the new matrix is bigger
     for (long i=1;i<=nr2;i++)
	  memcpy(D+N1*(i-1)+1,D2+oldN1*(i-1)+1,nr1*sizeof(cmplx));
     free(D2);
}

cmplx& CMatrix::Elem(long i, long j)
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  merror("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

cmplx CMatrix::Elem(long i, long j) const
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  merror("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

cmplx& CMatrix::operator()(long i, long j)
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  merror("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

cmplx CMatrix::operator() (long i, long j) const
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  merror("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

void CMatrix::Write(FILE *f) const
{
     for (long i=1;i<=N1;i++)
     {
         for (long j=1;j<=N2;j++)
	 {
	      cmplx x=Elem(i,j);
	      fprintf(f,"(%14.10g,%14.10g) ",real(x),imag(x));
	 }
         fprintf(f,"\n");
     }
     fprintf(f,"\n");
}

void CMatrix::Write() const
{
     Write(stdout);
}

void CMatrix::Write(int prec, FILE *f) const
{
     cmplx x;
     char form[40];
     sprintf(form,"(%%%d.%dg,%%%d.%dg) ",prec+6,prec,prec+6,prec);
     for (long i=1;i<=N1;i++)
     {
         for (long j=1;j<=N2;j++)
	 {
	      x=Elem(i,j);
	      fprintf(f,form,real(x),imag(x));
	 }
         fprintf(f,"\n");
     }
     fprintf(f,"\n");
}

void CMatrix::Write(int prec) const
{
     Write(prec,stdout);
}

CMatrix& CMatrix::operator=(const CMatrix& M)
{
     if (!(M.N1*M.N2)) 
     { 
	  Destroy(); 
	  return(*this); 
     }
     if (this==&M) return (*this);     
     Destroy();
     Copy(*this,M);
     return (*this);
}

CMatrix& CMatrix::operator=(const Matrix& M)
{
     if (!(M.N1*M.N2)) 
     { 
	  Destroy(); 
	  return(*this); 
     }
     Destroy();
     Copy(*this,M);
     return (*this);
}


void CMatrix::operator+=(const CMatrix& M)
{
     if (!N1) { Create(M.N1,M.N2); Zero(); }
#ifdef DEBUG
     if (N1!=M.N1 || N2!=M.N2) 
	  merror("Adding matrices with different dims.");
#endif
     cmplx alpha=1.0;
     long ix=1, n=N1*N2;
     zaxpy_(&n,&alpha,M.D+1,&ix,D+1,&ix);
}

void CMatrix::operator-=(const CMatrix& M)
{
     if (!N1) { Create(M.N1,M.N2); Zero(); }
#ifdef DEBUG
     if (N1!=M.N1 || N2!=M.N2) 
	  merror("Adding matrices with different dims.");
#endif
     cmplx alpha=-1.0;
     long ix=1, n=N1*N2;
     zaxpy_(&n,&alpha,M.D+1,&ix,D+1,&ix);
}

void CMatrix::operator*=(cmplx x)  
{
#ifdef DEBUG
    if (!N1) return;
#endif
    long ix=1, n=N1*N2;
    cmplx cosa=x;
    zscal_(&n,&cosa,D+1,&ix);
}

void CMatrix::operator*=(double x)
{
     (*this)*=(cmplx)x;
}

void CMatrix::operator*=(const CMatrix& M)
{
#ifdef DEBUG
     if (N2!=M.N1)
	  merror("Incompatible dimensions in *=");
     if (!N1 || !N2 || !M.N1 || !M.N2)
	  merror("Zero matrix in *=\n");
#endif
     cmplx alpha=1.0, beta=0.0;
     char ca='N', cb='N';
     long n1=N1, n2=M.N2, k=M.N1;
     cmplx *D2=(cmplx*)malloc((n1*n2+1)*sizeof(cmplx));
     if (n1<0 || n2<0) merror("Nonexistent matrix A!\n");
     zgemm_(&ca,&cb,&n1,&n2,&k,&alpha,D+1,&n1,M.D+1,&k,&beta,D2+1,&n1);
     Load(D2,n1,n2);     
}

CVector CMatrix::Col(long n) const
{
     CVector R;
     Col(R,n);
     return R;
}

CVector CMatrix::Row(long n) const
{
     CVector R;
     Row(R,n);
     return R;
}

void CMatrix::Col(CVector &R, long n) const
{
     cmplx *d=(cmplx*)malloc((N1+1)*sizeof(cmplx));
     memcpy(d+1,D+1+(n-1)*N1,N1*sizeof(cmplx));
     R.Load(d,N1);
}

void CMatrix::Put_Col(const CVector &V, long n)
{
#ifdef DEBUG
     if ((V.N!=N1) || (n>N2)) merror("Put_Col is impossible.");
#endif
     memcpy(D+(n-1)*N1+1,V.D+1,V.N*sizeof(cmplx));
}

void CMatrix::Row(CVector &R, long n) const
{
     R.Create(N2);
     for (long i=1;i<=N2;i++)
	  R(i)=Elem(n,i);
}

void CMatrix::Put_Row(const CVector &V, long n)
{
     for (long i=1;i<=N2;i++)
	  Elem(n,i)=V(i);
}

void CMatrix::Append_Col(const CVector &V)
{
     if (N1)
     	  Resize(N1,N2+1);
     else
     	  Create(V.N,1);
     Put_Col(V,N2);
}

void CMatrix::Append_Row(const CVector &V)
{
     if (N2)
     	  Resize(N1+1,N2);
     else
     	  Create(1,V.N);
     Put_Row(V,N1);
}

void CMatrix::Swap_Cols(long k1,long k2)
{
     cmplx acum;
     for(long i=1;i<=N1;i++)
     {
	  acum=Elem(i,k1);
	  Elem(i,k1)=Elem(i,k2);
	  Elem(i,k2)=acum;
     }
}

void CMatrix::Swap_Rows(long k1,long k2)
{
    cmplx acum;
    for(long i=1;i<=N2;i++)
    {
	 acum=Elem(k1,i);
	 Elem(k1,i)=Elem(k2,i);
	 Elem(k2,i)=acum;
    }
}

// Sort the columns of the matrix wrt the first element of each column
void CMatrix::Sort_Cols()
{
     qsort(D+1,N2,N1*sizeof(cmplx),complex_order_12);
}

// Isn't there a BLAS routine to do this?
cmplx CMatrix::Elem(const CVector &V1, const CVector &V2) const
{
     CVector W;
     Multiply(W,(*this),V2);
     return Dot(V1,W);
}

// Take the matrix element between 2 columns of matrix M
cmplx CMatrix::Elem(const CMatrix &M, long c1, long c2) const
{
#ifdef DEBUG
     if (N1!=N2) merror("Can't take matrix elem of non-square matrix.");
     if (N2!=M.N1) merror("Incompatible dimensions in Elem.");
#endif
     cmplx alpha=1.0, beta=0.0;
     long incx=1, n1=N1, n2=N2;
     char c='N';
     cmplx *d=(cmplx*)malloc(N1*sizeof(cmplx));
     zgemv_(&c, &n1, &n2, &alpha, D+1, &n1, M.D+N1*(c2-1)+1, &incx, 
	    &beta, d, &incx);
     long ix=1, n=M.N1;
     cmplx value=zdotc_(&n,M.D+N1*(c1-1)+1,&ix,d,&ix);     
     free(d);
     return value;
}

void CMatrix::Change_Basis(const CMatrix &B)
{
     (*this)=::Change_Basis(*this,B);
}

// Unitary basis change
CMatrix Change_Basis(const CMatrix &M, const CMatrix &B)
{
     long ncolB=B.N2; // number of columns in B
     CMatrix R(ncolB);
     for (long i=1;i<=ncolB;i++)
          for (long j=1;j<=ncolB;j++)
               R(i,j)=M.Elem(B,i,j);
     return R;
}

bool CMatrix::Gram_Schmidt()
{
     CVector V, W;
     cmplx dotprod;
     
     for (long k=1;k<=N2;k++)
     {
	  Col(V,k);
	  for (long j=1;j<=k-1;j++)
	  {
	       Col(W,j);
	       dotprod=-Dot(W,V);
	       Zaxpy(V,W,dotprod);
	  }
	  bool success=V.Normalize();
	  if (!success) return false;
	  Put_Col(V,k);
     }
     return true;
}

void CMatrix::T()
{
     (*this)=::T(*this);
}

CMatrix T(const CMatrix &M)
{
     CMatrix R(M.N2,M.N1);
     for (long i=1;i<=M.N2;i++)
          for (long j=1;j<=M.N1;j++)
               R(i,j)=M(j,i);
     return R;
}

void CMatrix::Herm()
{
     (*this)=::Herm(*this);
}

CMatrix Herm(const CMatrix &M)
{
     CMatrix R(M.N2,M.N1);
     for (long i=1;i<=M.N2;i++)
          for (long j=1;j<=M.N1;j++)
               R(i,j)=conj(M(j,i));
     return R;
}

void CMatrix::Part(long n10, long n20, long n1f, long n2f)
{
     CMatrix R;
     ::Part(R,*this,n10,n20,n1f,n2f);
     Transfer(R);
}

CMatrix Part(const CMatrix &M, long i0,long j0,long i1,long j1)
{
     CMatrix R;
     Part(R,M,i0,j0,i1,j1);
     return R;
}

void Part(CMatrix &R, const CMatrix &M, long n10, long n20, long n1f, long n2f)
{
     long m1=n1f-n10+1;
     long m2=n2f-n20+1;
     R.Create(m1,m2);
     for (long i=1;i<=m2;i++)
          memcpy(R.D+m1*(i-1)+1,M.D+M.N1*(n20+i-2)+n10,m1*sizeof(cmplx));
}


void CMatrix::Add(const CMatrix &M, long n1, long n2)
{
#ifdef DEBUG
     if (n1+M.N1>N1) merror("Incompatible dimensions in Add.");
     if (n2+M.N2>N2) merror("Incompatible dimensions in Add.");
#endif
     long m=n2-n1+1, n=M.N1, ix=1;
     cmplx alpha=1.0;
     for (long i=1;i<=m;i++)
	  zaxpy_(&n,&alpha,D+N1*(n1+i-1)+1,&ix,M.D+M.N1*(i-1)+1,&ix);
}

void CMatrix::Insert(const CMatrix &M, long i, long j)
{
#ifdef DEBUG
     if (i+M.N1-1>N1) merror("Incompatible dimensions in Insert.");
     if (j+M.N2-1>N2) merror("Incompatible dimensions in Insert.");
#endif
     long m=M.N2;
     for (long k=1;k<=m;k++)
	  memcpy(D+N1*(j+k-2)+i,M.D+M.N1*(k-1)+1,M.N1*sizeof(cmplx));
}

void CMatrix::Invert()
{
#ifdef DEBUG
     if (N1!=N2) merror("Can't invert non-square matrix!");
#endif
     CMatrix R(N1);
     R.Unit();
     Solve(R);
     Transfer(R);
}

CMatrix Invert(const CMatrix &M)
{
#ifdef DEBUG
     if (M.N1!=M.N2) merror("Can't invert non-square matrix!");
#endif     
     CMatrix R(M.N1);
     R.Unit();
     M.Solve(R);
     return R;
}

// OPTIMAL ROUTINE!!!
void CMatrix::Solve(CVector &R) const
{
#ifdef DEBUG
     if (R.N!=N1) merror("Incompatible dimensions in Solve routine.");
#endif
     long N=N1; // number of linear equations
     long M=1; // number of RHS's
     long info;
     cmplx *d=(cmplx*)malloc((N1*N2+1)*sizeof(cmplx));
     memcpy(d,D,(N1*N2+1)*sizeof(cmplx));
     long *ip=(long*)malloc(N*sizeof(long));
     zgesv_(&N, &M, d+1, &N, ip, R.D+1, &N, &info);
     free(d); free(ip);
     if (info) printf("Singular matrix! Result of Solve is worthless!\n");
}

CVector Solve(const CMatrix &M, const CVector &b)
{
     CVector B(b);
     M.Solve(B);
     return B;
}

// OPTIMAL ROUTINE!!!
void CMatrix::Solve(CMatrix &R) const
{
#ifdef DEBUG
     if (N1!=R.N1) merror("Incompatible dimensions in solver.");
#endif
     long N=N1; // number of linear equations
     long M=R.N2; // number of RHS's
     long info;
     cmplx *d=(cmplx*)malloc((N1*N2+1)*sizeof(cmplx));
     memcpy(d,D,(N1*N2+1)*sizeof(cmplx));
     long *ip=(long*)malloc(N*sizeof(long));
     zgesv_(&N, &M, d+1, &N, ip, R.D+1, &N, &info);
     free(ip); free(d);
     if (info) printf("Singular matrix! Result of Solve is worthless!\n");
}

// what is the return value in this function?
int CMatrix::LU_Decomp(int *I) 
{
     long M=N1;
     long N=N2;
     cmplx *A=D+1;
     long lda=M;
     long info;
     zgetrf_(&M,&N,A,&lda,I,&info);
     if (info<0) merror("Problems with LU decomp.");
     return info;
}

cmplx CMatrix::Det() const
{
#ifdef DEBUG
     if (N1!=N2) merror("Trying to get det of non-square matrix.");
#endif
     int *I=(int*)malloc(N1*sizeof(int));
     CMatrix B(*this);
     B.LU_Decomp(I);
     int signo=1;
     for (long i=1;i<=N1;i++)
       	  if (i!=I[i-1]) signo*=-1;
     cmplx det=1.0;
     for (long i=1;i<=N1;i++)
	  det*=B(i,i);
     det*=(cmplx)signo;
     free(I);
     return det;
}

cmplx Det(const CMatrix &M)
{
     return M.Det();
}

// so as not to make checks for "squarity", assumes that
// and uses N1
cmplx CMatrix::Trace() const
{
     cmplx result=0.0;
     for (long i=1;i<=N1;i++)
	  result+=Elem(i,i);
     return result;
}

cmplx Trace(const CMatrix &M)
{
     return M.Trace();
}

// returns dimension of Kernel of matrix
// In K, columns span the kernel, although they need not be orthogonal
long CMatrix::Kernel(CMatrix &K) const
{
     // Totally arbitrary... CAUTION
     double VERY_TINY=1e-6;

     int *indx=(int*)malloc((N1+1)*sizeof(int)); 
     CMatrix A(*this);
     A.LU_Decomp(indx);
     free(indx);

     long dim=0; // Kernel dimension
     // Positions of the zeroes in the diagonal of the LU decomposed matrix
     long* zero_pos=(long*)malloc((N1+1)*sizeof(long));  
     // Find dimension and fill in the positions
     for (long i=1;i<=N1;i++)
	  if (abs(A(i,i))<VERY_TINY) 
	  {
	       dim++;
	       zero_pos[dim]=i;
	  }
     
     // If dim=0, there is no kernel. 
     if (!dim) { free(zero_pos); return 0;}
              
     K.Create(N2,dim);
     K.Zero();
     for (long i=1;i<=dim;i++)
     {	  
	  //CVector V(N2);
	  K(zero_pos[i],i)=1.0; // Make up
	  for (long j=N2;j>=1;j--)
	       if (abs(A(j,j))>VERY_TINY)
	       {
		    cmplx sum=0.0;
		    for (long k=N1;k>j;k--)
			 sum+=A(j,k)*K(k,i);
		    K(j,i)=-sum/A(j,j);
	       }
     }    
     free(zero_pos);
     return dim;
}

long CMatrix::Eigenvectors(CMatrix &V, cmplx E) const
{
     CMatrix A(*this);
     for (long i=1;i<=N1;i++)
	  A(i,i)-=E;
     return A.Kernel(V);
}

// All I/O routines return 1 if ok, 0 if there was any error
int CMatrix::Save_Binary(FILE *fich) const
{
     int ausgang;
     ausgang=fwrite(&N1,sizeof(long),1,fich);
     if (ausgang!=1) return 0;
     ausgang=fwrite(&N2,sizeof(long),1,fich);
     if (ausgang!=1) return 0;
     if (!(N1*N2)) return 1;
     ausgang=fwrite(D,sizeof(cmplx),N1*N2+1,fich);
     if (ausgang!=N1*N2+1) return 0;
     return 1;
}

int CMatrix::Save_Binary(const char *s) const
{
     FILE *fich=fopen(s,"wb");
     if (!fich) return 0;
     int ausgang=Save_Binary(fich);
     fclose(fich);
     return ausgang;
}

int CMatrix::Load_Binary(FILE *fich)
{
     int ausgang;
     ausgang=fread(&N1,sizeof(long),1,fich);
     if (ausgang!=1) return 0;
     ausgang=fread(&N2,sizeof(long),1,fich);
     if (ausgang!=1) return 0;
     if (!(N1*N2)) { D=(cmplx*)NULL; return 1;}
     Create(N1,N2);
     ausgang=fread(D,sizeof(cmplx),N1*N2+1,fich);
     if (ausgang!=N1*N2+1) return 0;
     return 1;
}

int CMatrix::Load_Binary(const char *s) 
{
     FILE *fich=fopen(s,"rb");
     if (!fich) return 0;
     int ausgang=Load_Binary(fich);
     fclose(fich);
     return ausgang;
}

int CMatrix::Save(const char *name) const
{
     FILE *fich=fopen(name,"wt");
     if (!fich) return 0;
     int ausgang=Save(fich);
     fclose(fich);
     return ausgang;
}

int CMatrix::Load(const char *name)
{
     FILE *fich=fopen(name,"rt");
     if (!fich) return 0;
     int ausgang=Load(fich);
     fclose(fich);
     return ausgang;
}

int CMatrix::Save(FILE *fich) const
{
     fprintf(fich,"%ld %ld\n",N1,N2);
     for (long i=1;i<=N1;i++)
     {
	  for (long j=1;j<=N2;j++)
	  {
	       cmplx x=Elem(i,j);
	       fprintf(fich,"(%14.10g,%14.10g) ",real(x),imag(x));
	  }
	  fprintf(fich,"\n");
     }
     if (ferror(fich)) return 0;
     return 1;
}

int CMatrix::Load(FILE *fich)
{
     float x, y;
     long n1, n2;
     if (!fscanf(fich,"%ld %ld\n",&n1,&n2))
	  merror("Error in CMatrix::Load");
     Create(n1,n2);
     for (long i=1;i<=n1;i++)
	  for (long j=1;j<=n2;j++)
	  {
	       if (feof(fich)) return 0; 
	       if (!fscanf(fich,"(%g,%g)",&x,&y))
		    merror("Error in CMatrix::Load");
	       cmplx z(x,y);
	       Elem(i,j)=z;
	  }
     if (ferror(fich)) return 0;
     return 1;
}

////////////////////////////////////////////////////////////

// Full diagonalization: all eigenvalues, all eigenvectors
// NOT OPTIMAL: matrix is preserved 
// CAUTION: Hermitian matrix is assumed!!!!!!!!!!!!
void CMatrix::Diagonalize(CMatrix &B, Vector &E) const
{
#ifdef DEBUG
     if (N1!=N2) merror("Trying to diagonalize non-square matrix.");
#endif
     char jobz='V', range='A', uplo='U';
     long N=N1;
     cmplx *A=(cmplx*)malloc(N*N*sizeof(cmplx));
     memcpy(A,D+1,N*N*sizeof(cmplx));
     long lda=N;
     double vl, vu;
     long il, iu;
     char askmachine='S';
     double abstol=2.0*dlamch_(&askmachine); // This is the optimal absolute
     // tolerance to give, dlamch('S') is the minimum number whose inverse
     // does not overflow
     long m; // total number of eigenvalues found
     E.Create(N); double *W=E.D+1;
     B.Create(N); cmplx *Z=B.D+1;
     long ldz=N;
     cmplx *work=(cmplx*)malloc(8*N*sizeof(cmplx));
     long lwork=8*N;
     double *rwork=(double*)malloc(7*N*sizeof(double));
     long *iwork=(long*)malloc(5*N*sizeof(long));
     long *ifail=(long*)malloc(N*sizeof(long));
     long info;
     zheevx_(&jobz,&range,&uplo,&N,A,&lda,&vl,&vu,&il,&iu,&abstol,
	     &m,W,Z,&ldz,work,&lwork,rwork,iwork,ifail,&info);
#ifdef DEBUG
     if (info) printf("Error during diagonalization.");
#endif
     free(iwork);
     free(rwork);
     free(ifail);
     free(work);
     free(A);
//     printf("%ld eigenvalues obtained\n",m);
}


void CMatrix::SVD(CMatrix &BU, CMatrix &BVt, Vector &SV) const
{
     char jobu='S', jobvt='S';
     long M=N1, N=N2;
     cmplx *A=(cmplx*)malloc(M*N*sizeof(cmplx));
     memcpy(A,D+1,M*N*sizeof(cmplx));
     long lda=M;
     SV.Create(min(M,N)); double *S=SV.D+1;
     BU.Create(M,M); cmplx *U=BU.D+1;
     BVt.Create(N,N); cmplx *Vt=BVt.D+1;
     BU.Zero(); BVt.Zero();
     long ldu=M; long ldvt=N;
     long lwork=10*max(M,N);
     cmplx *work=(cmplx*)malloc(lwork*sizeof(cmplx));
     double *rwork=(double*)malloc(lwork*sizeof(double));
     long info;
     zgesvd_(&jobu,&jobvt,&M,&N,A,&lda,S,U,&ldu,Vt,&ldvt,work,&lwork, 
	     rwork,&info);
     free(A);
     free(work);
     free(rwork);
}

void CMatrix::Spectrum(Vector &E) const
{
#ifdef DEBUG
     if (N1!=N2) merror("Trying to diagonalize non-square matrix.");
#endif
     char jobz='N', range='A', uplo='U';
     long N=N1;
     cmplx *A=(cmplx*)malloc(N*N*sizeof(cmplx));
     memcpy(A,D+1,N*N*sizeof(cmplx));
     long lda=N;
     double vl, vu;
     long il, iu;
     char askmachine='S';
     double abstol=2.0*dlamch_(&askmachine); // This is the optimal abstol
     long m; // total number of eigenvalues found
     E.Create(N); double *W=E.D+1;
     //B.Create(N); real *Z=B.D+1;
     cmplx *Z=NULL;
     long ldz=N;
     cmplx *work=(cmplx*)malloc(8*N*sizeof(cmplx));
     long lwork=8*N;
     long *iwork=(long*)malloc(5*N*sizeof(long));
     double *rwork=(double*)malloc(7*N*sizeof(double));     
     long *ifail=(long*)malloc(N*sizeof(long));
     long info;
     zheevx_(&jobz,&range,&uplo,&N,A,&lda,&vl,&vu,&il,&iu,&abstol,
	     &m,W,Z,&ldz,work,&lwork,rwork,iwork,ifail,&info);
     free(A); free(work); free(iwork); free(rwork); free(ifail);
//     printf("info=%ld\n",info);
//     if (info) merror("Error during diagonalization.");
}
// void zheevx_(char* jobz, char* range, char* uplo, long* n, cmplx* A,
// 		  long* lda, double* vl, double* vu,
// 		  long* il,long* iu, double* abstol, long* M, double* W,
// 		  cmplx* Z, long* ldz, cmplx* work,
// 		  long* lwork,double* rwork,long* iwork,long*);


void CMatrix::Tridiagonalize(CMatrix &B, CVector &Diag, CVector &S) const
{
     // Needs dsytrd y dorgtr
     // Reduce to tri-diagonal form, basis is in "strange form"
     char uplo='U';
     long N=N1;
     cmplx *A=(cmplx*)malloc((N*N+1)*sizeof(cmplx));
     memcpy(A+1,D+1,N*N*sizeof(cmplx));
     long lda=N;
     Diag.Create(N); cmplx* d=Diag.D+1;
     S.Create(N-1); cmplx* e=S.D+1;
     cmplx* tau=(cmplx*)malloc(N*sizeof(cmplx));
     long lwork=N*N;
     cmplx* work=(cmplx*)malloc(lwork*sizeof(cmplx));
     long info;
     zhetrd_(&uplo, &N, A+1, &lda, d, e,
	     tau, work, &lwork, &info);
     if (info) merror("Error during tridiagonalization, stage 1.");
     // Compute the basis in "normal form" for tri-diagonal reduction
     zungtr_(&uplo, &N, A+1, &lda, tau, work, &lwork, &info);
     if (info) merror("Error during tridiagonalization, stage 2.");
     B.Create(N);
     B.Load(A,N,N);
     free(work); free(tau);

}

void CMatrix::Ns_Diagonalize(CMatrix &B, CVector &E) const
{
#ifdef DEBUG
     if (N1!=N2) merror("Can't Ns_Diagonalize non-square matrix.");
#endif
     
     char jobl='V';
     char jobr='N';
     long N=N1;
     cmplx *A=(cmplx*)malloc(N*N*sizeof(cmplx));
     memcpy(A,D+1,N*N*sizeof(cmplx));
     long lda=N;
     E.Create(N); 
     B.Create(N); 
     cmplx* w=E.D+1;
     cmplx* vl=B.D+1;
     long ldvl=N;
     cmplx* vr=(cmplx*)NULL;
     long ldvr=1;
     long lwork=6*N;
     cmplx* work=(cmplx*)malloc(lwork*sizeof(cmplx));
     double *rwork=(double*)malloc(2*N*sizeof(double));
     long info;
 
     zgeev_(&jobl,&jobr,&N,A,&lda,
	    w, vl, &ldvl, vr, &ldvr, 
	    work, &lwork, rwork, &info);
     free(A); free(work); free(rwork);
     if (info) merror("Error diagonalizing non-symmetric matrix.");
}


/////////////////////////////////////////////////////
// EXTERNAL MATRICIAL FUNCTIONS

void Copy(CMatrix& B, const CMatrix& A) // B <- A raw and strict copy. 
{
     B.Destroy();
     if (!A.N1) return;
     B.Start(A.N1,A.N2);
     memcpy(B.D+1,A.D+1,A.N1*A.N2*sizeof(cmplx));
}

void Copy(CMatrix &B, const Matrix &A)
{
     B.Destroy();
     if (!A.N1) return;
     B.Start(A.N1,A.N2);
     for (long i=1;i<=A.N1;i++)
	  for (long j=1;j<=A.N2;j++)
	       B(i,j)=A(i,j); // sorry, no other way...
}

CMatrix operator-(const CMatrix &M)
{
     CMatrix R(M);
     R*=(cmplx)-1.0;
     return R;
}

CMatrix operator+(const CMatrix &A, const CMatrix &B)
{
#ifdef DEBUG
     if (A.N1!=B.N1 || A.N2!=B.N2) 
	  merror("Adding matrices with different dimensions.");
#endif
     CMatrix R(A);
     R+=B;
     return R;
}

// It is correct like this, although it doesn't seem to be!!
// This thing takes its arguments opposite: B-A!!
CMatrix operator-(const CMatrix &A, const CMatrix &B)
{
#ifdef DEBUG
     if (A.N1!=B.N1 || A.N2!=B.N2) 
	  merror("Substracting matrices with different dimensions.");
#endif
     CMatrix R(A);
     R-=B;
     return R;
}

CMatrix operator*(cmplx K, const CMatrix &A)
{
     CMatrix R(A);
     R*=K;
     return R;
}

CMatrix operator*(double K, const CMatrix &A)
{
     return ((cmplx)K)*A;
}

CMatrix operator*(const CMatrix &A, cmplx K)
{
     CMatrix R(A);
     R*=K;
     return R;
}

CMatrix operator*(const CMatrix &A, double K)
{
     return ((cmplx)K)*A;
}


CVector operator*(const CMatrix &A, const CVector &V)
{
#ifdef DEBUG
     if (A.N2!=V.N) merror("Wrong dims in CMatrix-CVector product.");
#endif
     CVector R;
     Multiply(R,A,V);
     return R;
}

CMatrix operator*(const CMatrix &A, const CMatrix &B)
{
#ifdef DEBUG
     if (A.N2!=B.N1) merror("Wrong dims in CMatrix-CMatrix product.");
#endif
     CMatrix R;
     Multiply(R,A,B);
     return R;
}

void Multiply(CVector &R, const CMatrix &M, const CVector &V)
{
     cmplx alpha=1.0, beta=0.0;
     long incx=1, n1=M.N1, n2=M.N2;
     char c='N';
     R.Create(n1);
     zgemv_(&c, &n1, &n2, &alpha, M.D+1, &n1, V.D+1, &incx, 
	    &beta, R.D+1, &incx);
}
void Multiply(CMatrix &R, const CMatrix &M1, const CMatrix &M2)
{
     cmplx alpha=1.0, beta=0.0;
#ifdef DEBUG
     if (M1.N2!=M2.N1) 
	  merror("Incompatible dimensions in Multiply\n");
//     if (M1.N1==0 || M1.N2==0 || M2.N1==0 || M2.N2==0)
//	  merror("Zero matrix in Multiply\n");
#endif
     if (M1.N1==0 || M1.N2==0 || M2.N1==0 || M2.N2==0)
     {
	  R.Destroy();
	  return;
     }
	 
     char ca='N', cb='N';
     long n1=M1.N1, n2=M2.N2, k=M2.N1;
     cmplx *D2=(cmplx*)malloc((n1*n2+1)*sizeof(cmplx));
     if (n1<0 || n2<0) merror("Nonexistent matrix B!\n");
     zgemm_(&ca,&cb,&n1,&n2,&k,&alpha,M1.D+1,&n1,M2.D+1,&k,&beta,D2+1,&n1);
     R.Load(D2,n1,n2);     
}

// The most general CMatrix-CMatrix product routine
// R <- alpha*M1*M2 + beta*R
// M1 is transposed if T1=true
// M2 is transposed if T2=true
void Multiply_Add(CMatrix &R, const CMatrix &M1, const CMatrix &M2,
		  cmplx alpha, cmplx beta, bool T1, bool T2)
{
     long n1=(T1 ? M1.N2 : M1.N1), 
	  n2=(T2 ? M2.N1 : M2.N2), 
	  k=(T1 ? M1.N1 : M1.N2);
     if ((n1<0) || (n2<0)) 
	  printf("n1: %ld, n2: %ld\n",n1,n2);

#ifdef DEBUG
     long kk=(T2 ? M2.N2 : M2.N1);
     if (k!=kk) merror("Wrong internal dimensions in Multiply_Add");
     if (R.N1!=n1 || R.N2!=n2) merror("Wrong external dimensions in Multiply_Add");
     if (M1.N1==0 || M1.N2==0 || M2.N1==0 || M2.N2==0)
	  merror("Zero matrix in Multiply_Add\n");
#endif
     char ca=(T1==true ? 'C' : 'N');
     char cb=(T2==true ? 'C' : 'N');
     long N1=M1.N1;
     long N2=M2.N1;
     zgemm_(&ca,&cb,&n1,&n2,&k,&alpha,M1.D+1,&N1,M2.D+1,&N2,&beta,R.D+1,&n1);
}

// Diagonalize a tridiagonal matrix (diagonal D, subdiagonal S)
// and return the spectrum in D
// void Td_Spectrum(Vector &D, Vector &S)
// {
//      char compz='N';
//      long N=D.N;
//      cmplx *z=(cmplx*)NULL;
//      long ldz=N;
//      cmplx *work=(cmplx*)NULL;
//      long info;
//      zsteqr_(&compz,&N,D.D+1,S.D+1,z,&ldz,work,&info);
// }

// B should be already initialized, a Unit matrix
void Td_Diagonalize(CMatrix &B, Vector &D, Vector &S)
{
     char compz='V';
     long N=D.N;
     cmplx *z=B.D+1;
     long ldz=N;
     cmplx *work=(cmplx*)malloc(2*(N-1)*sizeof(cmplx)); 
     long info;
     zsteqr_(&compz,&N,D.D+1,S.D+1,z,&ldz,work,&info);
     free(work);
     if (info) merror("Error during Td_Diagonalize.");
}


CMatrix Tens_Prod(const CMatrix &A, const CMatrix &B)
{
     CMatrix R;
     Tens_Prod(R,A,B);
     return R;
}

void   Tens_Prod(CMatrix &R, const CMatrix &A, const CMatrix &B)
{
     long nA1=A.N1, nA2=A.N2, nB1=B.N1, nB2=B.N2;
     R.Create(nA1*nB1,nA2*nB2);
     for (long i1=1;i1<=nA1;i1++)
	for (long i2=1;i2<=nA2;i2++)
	   for (long j1=1;j1<=nB1;j1++)
	      for (long j2=1;j2<=nB2;j2++)
		 R( (i1-1)*nB1 + j1, (i2-1)*nB2 + j2 )=
		      A(i1,i2)*B(j1,j2);
}

CMatrix Tens_Prod_Unit(const CMatrix & A, long m)
{
     CMatrix R;
     Tens_Prod_Unit(R,A,m);
     return R;
}

void   Tens_Prod_Unit(CMatrix& R, const CMatrix& A, long m)
{
     long n1=A.N1, n2=A.N2;
     R.Create(n1*m,n2*m);
     R.Zero();
     for (long i1=1;i1<=n1;i1++)
	  for (long i2=1;i2<=n2;i2++)
	       for (long j=1;j<=m;j++)
		    R((i1-1)*m+j,(i2-1)*m+j)=A(i1,i2);
}

CMatrix Tens_Prod_Unit(long n, const CMatrix &A)
{
     CMatrix R;
     Tens_Prod_Unit(R,n,A);
     return R;
}

void   Tens_Prod_Unit(CMatrix& R, long n, const CMatrix& A)
{
     long m1=A.N1, m2=A.N2;
     R.Create(n*m1,n*m2);
     R.Zero();
     for (long i=1;i<=n;i++)
	  for (long j1=1;j1<=m1;j1++)
	       for (long j2=1;j2<=m2;j2++)
		    R((i-1)*m1+j1,(i-1)*m2+j2)=A(j1,j2);
}

// Unit is at side "s" (left or right)
void   Tens_Prod_Unit(CMatrix &R, const CMatrix &A, long n, Side s)
{
     if (s==Left)
	  Tens_Prod_Unit(R,n,A);
     else
	  Tens_Prod_Unit(R,A,n);
}

CMatrix Tens_Prod_Unit(const CMatrix &A, long n, Side s)
{
     if (s==Left)
	  return Tens_Prod_Unit(n,A);
     else
	  return Tens_Prod_Unit(A,n);
}

void Elem_Mult(CMatrix &R, const CMatrix &A, const CMatrix &B)
{
#ifdef DEBUG
     if (A.N1!=B.N1 || A.N2!=B.N2)
	  merror("Different dimensions in Elem_Mult\n");
#endif
     long N1=A.N1, N2=A.N2;
     if (R.N1!=N1 || R.N2!=N2)
	  R.Create(N1,N2);
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       R(i,j)=A(i,j)*B(i,j);
}

CMatrix Elem_Mult(const CMatrix &A, const CMatrix &B)
{
     CMatrix R(A.N1,A.N2);
     Elem_Mult(R,A,B);
     return R;
}

void Ket_Bra(CMatrix &R, const CVector &V, const CVector &W)
{
     long N1=V.N, N2=W.N;
     R.Create(N1,N2);
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       R(i,j)=conj(V(i))*W(j);
}


CMatrix Diag(const CVector &E)
{
     long N=E.N;
     CMatrix R(N);
     for (long i=1;i<=N;i++)
	  R(i,i)=E(i);
     return R;
}

// Order V and swap columns of M accordingly
// Hm... don't quite like the implementation...
void Order(CVector &V, CMatrix &M)
{
     // Uufff... we should make a documentation for this
     // Bubble algorithm, I'm afraid...
     cmplx acum;
     long nswaps;
     do
     {
	  nswaps=0;
	  for (long i=1;i<M.N2;i++)
	       if (abs(V(i))>abs(V(i+1)))
	       {
		    acum=V(i);
		    V(i)=V(i+1);
		    V(i+1)=acum;
		    M.Swap_Cols(i,i+1);
		    nswaps++;
	       }
     }while(nswaps);
}

void CMatrix::Real()
{
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=real(Elem(i,j));
}

CMatrix Real(const CMatrix &M)
{
     CMatrix R(M.N1,M.N2);
     for (long i=1;i<=M.N1;i++)
	  for (long j=1;j<=M.N2;j++)
	       R(i,j)=real(M(i,j));
     return R;
}

void CMatrix::Imag()
{
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=imag(Elem(i,j));
}

CMatrix Imag(const CMatrix &M)
{
     CMatrix R(M.N1,M.N2);
     for (long i=1;i<=M.N1;i++)
	  for (long j=1;j<=M.N2;j++)
	       R(i,j)=imag(M(i,j));
     return R;
}

void CMatrix::Abs()
{
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=abs(Elem(i,j));
}

CMatrix Abs(const CMatrix &M)
{
     CMatrix R(M.N1,M.N2);
     for (long i=1;i<=M.N1;i++)
	  for (long j=1;j<=M.N2;j++)
	       R(i,j)=abs(M(i,j));
     return R;
}

Matrix To_Real(const CMatrix &M)
{
     Matrix R(M.N1,M.N2);
     for (long i=1;i<=M.N1;i++)
	  for (long j=1;j<=M.N2;j++)
	       R(i,j)=(double)real(M(i,j));
     return R;
}

CMatrix To_Cmplx(const Matrix &M)
{
     CMatrix R(M.N1,M.N2);
     for (long i=1;i<=M.N1;i++)
	  for (long j=1;j<=M.N2;j++)
	       R(i,j)=(cmplx)M(i,j);
     return R;
}

double Frobenius_Norm(const CMatrix &M)
{
#ifdef DEBUG
     if (!M.N1 || !M.N2) merror("Frobenius norm of non-existent matrix\n");
#endif
     double norma=0.0;
     for (long i=1;i<=M.N1;i++)
	  for( long j=1;j<=M.N2;j++)
	       norma+=norm(M(i,j));
     return sqrt(norma);
}

#endif
