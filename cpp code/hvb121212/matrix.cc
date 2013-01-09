////////////////////////////////////////////////////////
// hvb++ 1.0
// Copyleft: Javier Rodríguez Laguna
// 060224-080725-101022-121004

#include "matrix.h"

#ifndef MATRIX
#define MATRIX


long num_vec=0;
long num_mat=0;
long total_memory=0;
long max_total_memory=0;
bool matrix_error=false;
void Mem_Control()
{
     // In order to catch any possible leaks
     printf("** Memory Control:\n");
     printf("Maximum number of doubles alloc'd: %ld\n",max_total_memory);
     printf("Total doubles allocated: %ld\n",total_memory);
     printf("Number of stored matrices: %ld\n", num_mat);          
     printf("Number of stored vectors: %ld\n\n",num_vec);
}

/////////////////////
// Vector creators
/////////////////////

Vector::Vector(): N(0), D(NULL) {} 

Vector::Vector(long n)
{
      Start(n);
      Zero();
}

Vector::Vector(const Vector &V)
{
     Start(V.N);
    if (N) memcpy(D,V.D,(N+1)*sizeof(double));
}

Vector::Vector(double *data, long n)
{
     Start(n);
     memcpy(D,data,(n+1)*sizeof(double));
}

Vector::~Vector() { Destroy(); }

// Call Create if it is sure that the Vector existed beforehand.
// It is destroyed, then allocated.
void Vector::Create(long n)
{
     Destroy();
     Start(n);
}

// Call Start if you assume that the Vector does not exist. 
// Useful when we're initializing an array of them.
void Vector::Start(long n)
{
     N=n;
     if (!N) { D=NULL; return;}
     D=(double*)malloc((n+1)*sizeof(double));    
     
#ifdef DEBUG 
     if (!D) merror("Error allocating vector."); 
     total_memory+=N;
     if (total_memory>max_total_memory) 
	  max_total_memory=total_memory;
     num_vec++;
#endif
}

// Load a double* into a Vector
void Vector::Load(double *d, long n)
{
     Destroy();
     N=n;
     D=d;
#ifdef DEBUG
     num_vec++;
     total_memory+=N;
     if (total_memory>max_total_memory) 
	  max_total_memory=total_memory;
#endif
}

// Copy a double* into a Vector
void Vector::Load_Copy(double *d1, long n)
{
     double *d2=(double*)malloc((n+1)*sizeof(double));
     memcpy(d2,d1,(n+1)*sizeof(double));
     Load(d2,n);
}

// The data of V go to our vector, V is left empty.
void Vector::Transfer(Vector &V)
{
     Destroy();
     D=V.D;
     N=V.N;
     V.D=(double*)NULL;
     V.N=0; 
}

void Vector::Destroy()
{
     if (N) 
     { 	
	  free(D);
#ifdef DEBUG
	  total_memory-=N;
	  num_vec--;
#endif
	  N=0;
	  D=NULL;
     }
}


/////////////////////////////////////////
// Vector methods: Element handling
/////////////////////////////////////////

double Vector::operator() (long n) const
{
#ifdef DEBUG
     if (n<0 || n>N) merror("Error getting vector comp."); 
#endif
    return D[n];
}

double& Vector::operator() (long n)
{
#ifdef DEBUG
    if (n<0 || n>N) merror("Error putting vector comp."); 
#endif
    return D[n];
}

void Vector::Zero()
{
     if (N) memset(D,0,(N+1)*sizeof(double));
}

// Find if the vector is zero within a given tolerance
bool Vector::Is_Zero(double tolerance)
{
     for (long i=1;i<=N;i++)
	  if (fabs(D[i])>tolerance) return false;
     return true;
}

// Set the value to a given one.
void Vector::Set_Value(double x)
{
     if (N) for (long i=1;i<=N;i++) D[i]=x;
}

void Vector::Append(const double x) 
{
     if (!N) Create(1);
     else
     {
	  N++;
	  D=(double*)realloc(D,(N+1)*sizeof(double));
#ifdef DEBUG
	  total_memory++;
	  if (total_memory>max_total_memory)
	       max_total_memory=total_memory;
#endif
     }
     D[N]=x;     
}

void Vector::Append(const Vector &V)
{
     long nold=N;
     if (!N) Create (V.N);
     else
     {
          N+=V.N;
          D=(double*)realloc(D,(N+1)*sizeof(double));
#ifdef DEBUG
          total_memory+=V.N;
          if (total_memory>max_total_memory)
               max_total_memory=total_memory;
#endif
     }
     for (long i=1;i<=V.N;i++)
          D[nold+i]=V(i);
}

// Insert x at position i, Vector size increases by 1!
void Vector::Insert(const double x, long i)
{
     N++;
     D=(double*)realloc(D,(N+1)*sizeof(double));
#ifdef DEBUG
     total_memory++;
     if (total_memory>max_total_memory)
	  max_total_memory=total_memory;
#endif
     memmove(D+i+1,D+i,(N-i)*sizeof(double));
     D[i]=x;
}

void Vector::Remove(long i)
{
     memmove(D+i,D+i+1,(N-i-1)*sizeof(double));
     N--;
     D=(double*)realloc(D,(N+1)*sizeof(double));
#ifdef DEBUG
     total_memory--;
#endif
}

void Copy(Vector& B, const Vector& A)
{
     if (B.N) B.Destroy(); 
     if (!A.N) return;
     B.Start(A.N);
     memcpy(B.D,A.D,(A.N+1)*sizeof(double));
}

/////////////////////////
// Vector I/O routines
/////////////////////////

void Vector::Write_Col(FILE *fich) const
{
     for (long i=1;i<=N;i++)
	  fprintf(fich,"%12.8g \n",D[i]);
     fprintf(fich,"\n");
}

void Vector::Write_Col() const
{
     Write_Col(stdout);
}

void Vector::Write(FILE *fich) const
{
     for (long i=1;i<=N;i++)
	  fprintf(fich,"%12.8g ",D[i]);
     fprintf(fich,"\n\n");
}

void Vector::Write() const
{
     Write(stdout);
}

int Vector::Save_Binary(FILE *fich) const
{
     int ausgang;
     ausgang=fwrite(&N,sizeof(long),1,fich);
     if (ausgang!=1) return 0;
     if (!N) return 1;
     ausgang=fwrite(D,sizeof(double),N+1,fich);
     if (ausgang!=N+1) return 0;
     return 1;
}

int Vector::Save_Binary(const char *name) const
{
     FILE *fich=fopen(name,"wb");
     if (!fich) return 0;
     int status=Save_Binary(fich);
     fclose(fich);
     return status;
}

int Vector::Load_Binary(FILE *fich)
{
     int ausgang;
     ausgang=fwrite(&N,sizeof(long),1,fich);
     if (ausgang!=1) return 0;
     if (!N) { D=(double*)NULL; return 1; }
     Create(N);
     ausgang=fwrite(D,sizeof(double),N+1,fich);
     if (ausgang!=N+1) return 0;
     return 1;
}

int Vector::Load_Binary(const char *name)
{
     FILE *fich=fopen(name,"rb");
     if (!fich) return 0;
     int status=Load_Binary(fich);
     fclose(fich);
     return status;
}

int Vector::Save(const char *name) const
{
     FILE *fich=fopen(name,"wt");
     if (!fich) return 0;
     int ausgang=Save(fich);
     fclose(fich);
     return ausgang;
}

int Vector::Save(FILE *fich) const
{
//     fprintf(fich,"%ld\n",N);
     for (long i=1;i<=N;i++)
	  fprintf(fich,"%14.10g\n",D[i]);
     return 1;
}

int Vector::Load(const char *name) 
{
     FILE *fich=fopen(name,"rt");
     if (!fich) return 0;
     int ausgang=Load(fich);
     fclose(fich);
     return ausgang;
}

int Vector::Load(FILE *fich)
{
     Create(0);
     float x;
     while(!feof(fich))
     {
	  if (!fscanf(fich,"%g",&x))
	       merror("File reading error in Vector::Load");
	  Append(x);
     }
     return 1;
}

/////////////////////////
// Vector operators
/////////////////////////

Vector& Vector::operator=(const Vector& W)
{
     if (!W.N) { Destroy(); return (*this); }
     if (this==&W) return *this;
     Copy(*this,W);
     return (*this);
}

void Vector::operator+=(const Vector& W)
{
     if (!N) { Create(W.N); Zero(); }
#ifdef DEBUG
      if (N!=W.N) merror("Incompatible sizes in vector +=");
#endif
      double alpha=1.0;
      long ix=1;
      daxpy_(&N,&alpha,W.D+1,&ix,D+1,&ix);
}

void Vector::operator-=(const Vector& W)
{
     if (!N) { Create(W.N); Zero(); }
#ifdef DEBUG
     if (N!=W.N) merror("Incompatible sizes in vector -=");
#endif
     double alpha=-1.0;
     long ix=1;
     daxpy_(&N,&alpha,W.D+1,&ix,D+1,&ix);
}

void Vector::operator*=(const double x)
{
#ifdef DEBUG
    if (!N) return;
#endif
    long ix=1;
    double cosa=x;
    dscal_(&N,&cosa,D+1,&ix);
}

void Vector::operator/=(const double x)
{
#ifdef DEBUG
     if (!N) return;
#endif
     long ix=1;
     double cosa=1.0/x;
     dscal_(&N,&cosa,D+1,&ix);
}

Vector operator-(const Vector &M)
{
     Vector R(M);
     R*=-1.0;
     return R;
}

Vector operator+(const Vector &A, const Vector &B)
{
#ifdef DEBUG
     if (A.N!=B.N) 
	  merror("Adding vectors with different dimensions.");
#endif
     Vector R(A);
     R+=B;
     return R;
}

// It is correct like this, although it doesn't seem to be!!
// This thing takes its arguments opposite: B-A!!
Vector operator-(const Vector &A, const Vector &B)
{
#ifdef DEBUG
     if (A.N!=B.N) 
	  merror("Substracting vectors with different dimensions.");
#endif
     Vector R(A);
     R-=B;
     return R;
}

Vector operator*(double x, const Vector& V)
{
     Vector R(V);
     R*=x;
     return R;
}

Vector operator*(const Vector& V, double x)
{
     Vector R(V);
     R*=x;
     return R;
}

Vector operator/(const Vector& V, double x)
{
     Vector R(V);
     R/=x;
     return R;
}

/////////////////////////////////
/// Vector functions
/////////////////////////////////

double Vector::Norm() const
{
     long ix=1, n=N;
     double norm=ddot_(&n,D+1,&ix,D+1,&ix);
     return sqrt(norm);
}

double Norm(const Vector &V)
{
     return V.Norm();
}

// Absolute lowest value
double Vector::Min() const
{
     double E0=D[1];
     for (long i=2;i<=N;i++)
	  if (D[i]<E0) E0=D[i];
     return E0;
}

// Get the k-th lowest value... not very fast, but effective
double Vector::Min(long k) const
{
     Vector Acum;
     Copy(Acum,*this);
     Acum.Order();
     return Acum(k);
}

long Vector::Min_Index() const
{
     double E0=D[1];
     long imin=1;
     for (long i=2;i<=N;i++)
	  if (D[i]<E0) { E0=D[i]; imin=i; }
     return imin;
}

double Vector::Max() const
{
     double E0=D[1];
     for (long i=2;i<=N;i++)
	  if (D[i]>E0) E0=D[i];
     return E0;
}

double Vector::Max(long k) const
{
     Vector Acum;
     Copy(Acum,*this);
     Acum.Order(-1);
     return Acum(k);
}

long Vector::Max_Index() const
{
     double E0=D[1];
     long imax=1;
     for (long i=2;i<=N;i++)
	  if (D[i]>E0) { E0=D[i]; imax=i; }
     return imax;
}

double Vector::Sum() const
{
     if (!N) return 0.0;
     double sum=D[1];
     for (long i=2;i<=N;i++)
	  sum+=D[i];
     return sum;
}

// Sum from i1 to i2
double Vector::Sum(long i1, long i2) const
{
     double sum=D[i1];
     for (long i=i1+1;i<=i2;i++)
	  sum+=D[i];
     return sum;
}

double Vector::Average() const
{
     return Sum()/(double)N;
}

double Vector::Variance() const
{
     if (N==1) return 0.0;
     double sumsq=0.0;
     for (long i=1;i<=N;i++)
          sumsq+=sqr(D[i]);
     sumsq/=(double)N;
     double aver=Average();
     return sumsq-sqr(aver);
}

double Vector::Deviation() const
{
     return sqrt(Variance());
}

double Min(const Vector &V)
{
     return V.Min();
}

double Max(const Vector &V)
{
     return V.Max();
}

double Sum(const Vector &V)
{
     return V.Sum();
}

double Average(const Vector &V)
{
     return V.Average();
}

double Deviation(const Vector &V)
{
     return V.Deviation();
}

double Variance(const Vector &V)
{
     return V.Variance();
}

int order_12(const void* x, const void* y) // ascending order
{
     const double *dx=(const double *)x;
     const double *dy=(const double *)y;
     return (*dx > *dy) - (*dx < *dy);
}

int order_21(const void* x, const void* y) // descending order
{
     const double *dx=(const double *)x;
     const double *dy=(const double *)y;
     return (*dx < *dy) - (*dx > *dy);
}

void Vector::Order(int p)
{     
     if (p<0)
	  qsort(D+1,N,sizeof(double),order_21);
     else
	  qsort(D+1,N,sizeof(double),order_12);
}

double Vector::Last() const
{
     return D[N];
}

// return false if normalization was not possible!
bool Vector::Normalize()
{
     double norm=Norm();
     if (norm==0.0) return false;
     (*this)/=norm;
     return true;
}

Vector Normalize(const Vector &V)
{
     Vector R(V);
     R.Normalize();
     return R;
}

void Vector::Part(long n1, long n2)
{
#ifdef DEBUG
     if (n1<1 || n2>N) 
	  merror("The part can't be larger than the whole!");
#endif
     Vector R(n2-n1+1);
     memcpy(R.D+1,D+n1,(n2-n1+1)*sizeof(double));
     Transfer(R);
}

Vector Part(const Vector &V, long n1, long n2) 
{
     Vector R(V);
     R.Part(n1,n2);
     return R;
}

void Vector::Reverse()
{
     for (long i=1;i<=N/2;i++)
	  SWAP(D[i],D[N+1-i]);
}

Vector Reverse(const Vector &V)
{
     Vector R(V);
     R.Reverse();
     return R;
}

void Vector::Insert(const Vector& V, long n)
{
#ifdef DEBUG
     if (n+V.N-1>N)
	  merror("Inserting a vector which is too big.");
#endif
     memcpy(D+n,V.D+1,V.N*sizeof(double));
}

Vector Insert(const Vector &V, const Vector &W, long n)
{
     Vector R(V);
     R.Insert(W,n);
     return R;
}

void Vector::Sqr()
{
     for (long i=1;i<=N;i++)
	  D[i]=sqr(D[i]);
}

Vector Sqr(const Vector &V)
{
     Vector R(V);
     R.Sqr();
     return R;
}

double Dot(const Vector& V1, const Vector& V2)
{
#ifdef DEBUG
     if (V1.N!=V2.N) merror("Dot product of vectors of diferent dim.");
#endif
     long ix=1, n=V1.N;
     return ddot_(&n,V1.D+1,&ix,V2.D+1,&ix);
}

Vector Cross(const Vector& V, const Vector& W)
{
#ifdef DEBUG
     if (V.N!=3 || W.N!=3) merror ("Cross product: vector dim must be 3.");
#endif
     Vector R(3);
     R(1)=V(2)*W(3)-V(3)*W(2);
     R(2)=V(3)*W(1)-V(1)*W(3);
     R(3)=V(1)*W(2)-V(2)*W(1);
     return R;
}

Vector Tens_Prod(const Vector &V, const Vector &W)
{
     long N=V.N*W.N;
     Vector R(N);
     long k=1;
     for (long i=1;i<=V.N;i++)
	  for (long j=1;j<=W.N;j++)
	  {
	       R.D[k]=V.D[i]*W.D[j];
	       k++;
	  }
     return R;
}

void Tens_Prod(Vector &R, const Vector &V, const Vector &W)
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

Vector Elem_Mult(const Vector &A, const Vector &B)
{
     Vector V(A.N);
     Elem_Mult(V,A,B);
     return V;
}

void Elem_Mult(Vector &V, const Vector &A, const Vector &B)
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
void Daxpy(Vector &R, const Vector &V, double alpha)
{
     long ix=1, n=R.N;
     daxpy_(&n,&alpha,V.D+1,&ix,R.D+1,&ix);     
}

///////////////////////////////////////////////////////////////

////////////////////////////////
// Matrix creators
////////////////////////////////

Matrix::Matrix(long n1, long n2) : N1(n1), N2(n2)
{
     N1=N2=0;
     Start(n1,n2);
     Zero();
}

Matrix::Matrix(const Matrix & M)
{
     N1=N2=0; D=(double*)NULL;
     Copy(*this,M);
}

Matrix::Matrix()
{
     N1=N2=0;
     D=NULL;
}

Matrix::~Matrix()
{
     Destroy();
}

void Matrix::Create(long n1, long n2)
{
     Destroy();
     Start(n1, n2);
}

void Matrix::Start(long n1, long n2) // n2=0
{
     N1=n1; N2=n2;
     if (!N2) N2=N1;
     if (!N1)
     {
	  D=NULL;
	  return;
     }
     D=(double*)malloc((N1*N2+1)*sizeof(double));
#ifdef DEBUG
     if (!D) merror ("Error allocating matrix.");
     total_memory+=N1*N2;
     if (total_memory>max_total_memory)
	  max_total_memory=total_memory;
     num_mat++;
#endif
}

// CAUTION: the data must be stored in columns
void Matrix::Load(double* d, long n1, long n2)
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
void Matrix::Load_Copy(double *d1, long n1, long n2)
{
     if (!n2) n2=n1;
     double *d2=(double*)malloc((n1*n2+1)*sizeof(double));
     memcpy(d2,d1,(n1*n2+1)*sizeof(double));
     Load(d2,n1,n2);
}

void Matrix::Transfer(Matrix &M)
{
     Destroy();
     D=M.D;
     N1=M.N1;
     N2=M.N2;
     M.N1=M.N2=0;
     M.D=(double*)NULL;
}

void Matrix::Destroy()
{
     if (N1) free(D);
     D=NULL;
#ifdef MEM_CONTROL
     total_memory-=N1*N2;
     if (N1*N2) num_mat--;
#endif
     N1=N2=0;
}


/////////////////////////////////////
// Matrix element handling
/////////////////////////////////////

void Matrix::Zero()
{
     memset(D+1,0,N1*N2*sizeof(double));
}

bool Matrix::Is_Zero(double tolerance)
{
     for (long i=1;i<=N1*N2;i++)
	  if (fabs(D[i])>tolerance) return false;
     return true;
}

void Matrix::Unit()
{
     Zero();
     long n=MIN(N1,N2);
     for (long i=1;i<=n;i++)
	  Elem(i,i)=1.0;
}

void Matrix::Set_Value(double x)
{
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       Elem(i,j)=x;
}

void Matrix::Resize(long n1, long n2) // m=0
{
     if (!n2) n2=n1;
     if (n1==N1 && n2==N2) return; 
     long nr1=MIN(N1,n1), nr2=MIN(N2,n2);
     long oldN1=N1, oldN2=N2;
     
     double *D2=(double*)malloc((N1*N2+1)*sizeof(double));
     memcpy(D2+1,D+1,N1*N2*sizeof(double));
     
     Create(n1,n2);
     if (n1>oldN1 || n2>oldN2)
	  Zero(); // In case the new matrix is bigger
     for (long i=1;i<=nr2;i++)
	  memcpy(D+N1*(i-1)+1,D2+oldN1*(i-1)+1,nr1*sizeof(double));
     free(D2);
}

double& Matrix::Elem(long i, long j)
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  merror("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

double Matrix::Elem(long i, long j) const
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  merror("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

double& Matrix::operator()(long i, long j)
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  merror("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

double Matrix::operator() (long i, long j) const
{
#ifdef DEBUG
     if (i<1 || i>N1 || j<1 || j>N2) 
	  merror("Error accessing matrix element.");
#endif
     return D[N1*(j-1)+i];
}

void Matrix::Write(FILE *f) const
{
     for (long i=1;i<=N1;i++)
     {
         for (long j=1;j<=N2;j++)
	 {
	      double x=Elem(i,j);
	      fprintf(f,"%12.8g ",x);
	 }
         fprintf(f,"\n");
     }
     fprintf(f,"\n");
}

void Matrix::Write() const
{
     Write(stdout);
}

Matrix& Matrix::operator=(const Matrix& M)
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

void Matrix::operator+=(const Matrix& M)
{
     if (!N1) { Create(M.N1,M.N2); Zero(); }
#ifdef DEBUG
     if (N1!=M.N1 || N2!=M.N2) 
	  merror("Adding matrices with different dims.");
#endif
     double alpha=1.0;
     long ix=1, n=N1*N2;
     daxpy_(&n,&alpha,M.D+1,&ix,D+1,&ix);
}

void Matrix::operator-=(const Matrix& M)
{
     if (!N1) { Create(M.N1,M.N2); Zero(); }
#ifdef DEBUG
     if (N1!=M.N1 || N2!=M.N2) 
	  merror("Adding matrices with different dims.");
#endif
     double alpha=-1.0;
     long ix=1, n=N1*N2;
     daxpy_(&n,&alpha,M.D+1,&ix,D+1,&ix);
}

void Matrix::operator*=(double x)  
{
#ifdef DEBUG
    if (!N1) return;
#endif
    long ix=1, n=N1*N2;
    double cosa=x;
    dscal_(&n,&cosa,D+1,&ix);
}

void Matrix::operator*=(const Matrix& M)
{
#ifdef DEBUG
     if (N2!=M.N1)
	  merror("Incompatible dimensions in *=");
     if (!N1 || !N2 || !M.N1 || !M.N2)
	  merror("Zero matrix in *=\n");
#endif
     double alpha=1.0, beta=0.0;
     char ca='N', cb='N';
     long n1=N1, n2=M.N2, k=M.N1;
     double *D2=(double*)malloc((n1*n2+1)*sizeof(double));
     dgemm_(&ca,&cb,&n1,&n2,&k,&alpha,D+1,&n1,M.D+1,&k,&beta,D2+1,&n1);
     Load(D2,n1,n2);     
}

Vector Matrix::Col(long n) const
{
     Vector R;
     Col(R,n);
     return R;
}

Vector Matrix::Row(long n) const
{
     Vector R;
     Row(R,n);
     return R;
}

void Matrix::Col(Vector &R, long n) const
{
     double *d=(double*)malloc((N1+1)*sizeof(double));
     memcpy(d+1,D+1+(n-1)*N1,N1*sizeof(double));
     R.Load(d,N1);
}

void Matrix::Put_Col(const Vector &V, long n)
{
#ifdef DEBUG
     if ((V.N!=N1) || (n>N2)) merror("Put_Col is impossible.");
#endif
     memcpy(D+(n-1)*N1+1,V.D+1,V.N*sizeof(double));
}

void Matrix::Row(Vector &R, long n) const
{
     R.Create(N2);
     for (long i=1;i<=N2;i++)
	  R(i)=Elem(n,i);
}

void Matrix::Put_Row(const Vector &V, long n)
{
     for (long i=1;i<=N2;i++)
	  Elem(n,i)=V(i);
}

void Matrix::Append_Col(const Vector &V)
{
     if (N1)
     	  Resize(N1,N2+1);
     else
     	  Create(V.N,1);
     Put_Col(V,N2);
}

void Matrix::Append_Row(const Vector &V)
{
     if (N2)
     	  Resize(N1+1,N2);
     else
     	  Create(1,V.N);
     Put_Row(V,N1);
}

void Matrix::Swap_Cols(long k1,long k2)
{
     double acum;
     for(long i=1;i<=N1;i++)
     {
	  acum=Elem(i,k1);
	  Elem(i,k1)=Elem(i,k2);
	  Elem(i,k2)=acum;
     }
}

void Matrix::Swap_Rows(long k1,long k2)
{
    double acum;
    for(long i=1;i<=N2;i++)
    {
	 acum=Elem(k1,i);
	 Elem(k1,i)=Elem(k2,i);
	 Elem(k2,i)=acum;
    }
}

// Sort the columns of the matrix wrt the first element of each column
void Matrix::Sort_Cols()
{
     qsort(D+1,N2,N1*sizeof(double),order_12);
}

// Isn't there a BLAS routine to do this?
double Matrix::Elem(const Vector &V1, const Vector &V2) const
{
     return Dot(V1,(*this)*V2);
}

// Take the matrix element between 2 columns of matrix M
double Matrix::Elem(const Matrix &M, long c1, long c2) const
{
#ifdef DEBUG
     if (N1!=N2) merror("Can't take matrix elem of non-square matrix.");
     if (N2!=M.N1) merror("Incompatible dimensions in Elem.");
#endif
     double alpha=1.0, beta=0.0;
     long incx=1, n1=N1, n2=N2;
     char c='N';
     double *d=(double*)malloc(N1*sizeof(double));
     dgemv_(&c, &n1, &n2, &alpha, D+1, &n1, M.D+N1*(c2-1)+1, &incx, 
	    &beta, d, &incx);
     long ix=1, n=M.N1;
     double value=ddot_(&n,M.D+N1*(c1-1)+1,&ix,d,&ix);     
     free(d);
     return value;
}

void Matrix::Change_Basis(const Matrix &B)
{
     (*this)=::Change_Basis(*this,B);
}

// Orthogonal basis change
Matrix Change_Basis(const Matrix &M, const Matrix &B)
{
     long ncolB=B.N2; // number of columns in B
     Matrix R(ncolB);
     for (long i=1;i<=ncolB;i++)
	  for (long j=1;j<=ncolB;j++)
	       R(i,j)=M.Elem(B,i,j);
     return R;
}

void Matrix::Gram_Schmidt()
{
     Vector V, W;
     double dotprod;
     
     for (long k=1;k<=N2;k++)
     {
	  Col(V,k);
	  for (long j=1;j<=k-1;j++)
	  {
	       Col(W,j);
	       dotprod=-Dot(V,W);
	       Daxpy(V,W,dotprod);
	  }
	  bool success=V.Normalize();
	  if (!success) 
	  {
	       matrix_error=true;
	       return;
	  }
	  Put_Col(V,k);
     }
}

void Matrix::T()
{
     (*this)=::T(*this);
}

Matrix T(const Matrix &M)
{
     Matrix R(M.N2,M.N1);
     for (long i=1;i<=M.N2;i++)
	  for (long j=1;j<=M.N1;j++)
	       R(i,j)=M(j,i);
     return R;
}


void Matrix::Part(long n10, long n20, long n1f, long n2f)
{
     Matrix R;
     ::Part(R,*this,n10,n20,n1f,n2f);
     Transfer(R);
}

Matrix Part(const Matrix &M, long i0,long j0,long i1,long j1)
{
     Matrix R;
     Part(R,M,i0,j0,i1,j1);
     return R;
}

void Part(Matrix &R, const Matrix &M, long n10, long n20, long n1f, long n2f)
{
#ifdef DEBUG
     if (n1f<n10 || n2f<n20) merror("End before start in Part.\n");
     if (n10<1 || n20<1) merror("Wrong boundaries in Part (begin).\n");
     if (n1f>M.N1 || n2f>M.N2) merror("Wrong boundaries in Part (end).\n");
#endif
     
     long m1=n1f-n10+1;
     long m2=n2f-n20+1;
     R.Create(m1,m2);
     for (long i=0;i<m2;i++)
	  memcpy(R.D+m1*i+1,
		 M.D+M.N1*(n20+i-1)+n10,
		 m1*sizeof(double));
}


void Matrix::Add(const Matrix &M, long n1, long n2)
{
#ifdef DEBUG
     if (n1+M.N1>N1) merror("Incompatible dimensions in Add.");
     if (n2+M.N2>N2) merror("Incompatible dimensions in Add.");
#endif
     long m=n2-n1+1, n=M.N1, ix=1;
     double alpha=1.0;
     for (long i=1;i<=m;i++)
	  daxpy_(&n,&alpha,D+N1*(n1+i-1)+1,&ix,M.D+M.N1*(i-1)+1,&ix);
}

void Matrix::Insert(const Matrix &M, long i, long j)
{
#ifdef DEBUG
     if (i+M.N1-1>N1) merror("Incompatible dimensions in Insert.");
     if (j+M.N2-1>N2) merror("Incompatible dimensions in Insert.");
#endif
     long m=M.N2;
     for (long k=1;k<=m;k++)
	  memcpy(D+N1*(j+k-2)+i,M.D+M.N1*(k-1)+1,M.N1*sizeof(double));
}

void Matrix::Invert()
{
#ifdef DEBUG
     if (N1!=N2) merror("Can't invert non-square matrix!");
#endif
     Matrix R(N1);
     R.Unit();
     Solve(R);
     Transfer(R);
}

// OPTIMAL ROUTINE
void Matrix::Solve(Vector &R) const
{
#ifdef DEBUG
     if (R.N!=N1) merror("Incompatible dimensions in Solve routine.");
#endif
     long N=N1; // number of linear equations
     long M=1; // number of RHS's
     long info;
     double *d=(double*)malloc((N1*N2+1)*sizeof(double));
     memcpy(d,D,(N1*N2+1)*sizeof(double));
     long *ip=(long*)malloc(N*sizeof(long));
     dgesv_(&N, &M, d+1, &N, ip, R.D+1, &N, &info);
     free(d); free(ip);
//     if (info) printf("Singular matrix! Result of Solve is worthless!\n");
}

Vector Solve(const Matrix &M, const Vector &b)
{
     Vector B(b);
     M.Solve(B);
     return B;
}

// OPTIMAL ROUTINE
// Solve Ax=b for all columns of matrix R, put results as columns of R.
void Matrix::Solve(Matrix &R) const
{
#ifdef DEBUG
     if (N1!=R.N1) merror("Incompatible dimensions in solver.");
#endif
     long N=N1; // number of linear equations
     long M=R.N2; // number of RHS's
     long info;
     double *d=(double*)malloc((N1*N2+1)*sizeof(double));
     memcpy(d,D,(N1*N2+1)*sizeof(double));
     long *ip=(long*)malloc(N*sizeof(long));
     dgesv_(&N, &M, d+1, &N, ip, R.D+1, &N, &info);
     free(ip); free(d);
     if (info) printf("Singular matrix! Result of Solve is worthless!\n");
}

// what is the return value in this function?
int Matrix::LU_Decomp(int *I) 
{
     long M=N1;
     long N=N2;
     double *A=D+1;
     long lda=M;
     long info;
     dgetrf_(&M,&N,A,&lda,I,&info);
     printf("Printing I:\n");
     for (long i=0;i<=N1;i++)
	  printf("I[%ld]=%d\n",i,I[i]);
     if (info<0) merror("Problems with LU decomp.");
     return info;
}

double Matrix::Det() const
{
#ifdef DEBUG
     if (N1!=N2) merror("Trying to get det of non-square matrix.");
#endif
     int *I=(int*)malloc(N1*sizeof(int));
     Matrix B(*this);
     B.LU_Decomp(I);
     int signo=1;
     for (long i=1;i<=N1;i++)
     	  if (i!=I[i-1]) signo*=-1;
     double det=1.0;
     for (long i=1;i<=N1;i++)
	  det*=B(i,i);
     det*=(double)signo;
     free(I);
     return det;
}

// so as not to make checks for "squarity", assumes that
// and uses N1
double Matrix::Trace() const
{
#ifdef DEBUG
     if (N1!=N2) merror("Trace of rectangular matrix");
#endif
     double result=0.0;
     for (long i=1;i<=N1;i++)
	  result+=Elem(i,i);
     return result;
}

// returns dimension of Kernel of matrix
// In K, columns span the kernel, although they need not be orthogonal
long Matrix::Kernel(Matrix &K) const
{
     // Totally arbitrary... CAUTION
     double VERY_TINY=1e-6;

     int *indx=(int*)malloc((N1+1)*sizeof(int)); 
     Matrix A(*this);
     A.LU_Decomp(indx);
     free(indx);

     long dim=0; // Kernel dimension
     // Positions of the zeroes in the diagonal of the LU decomposed matrix
     long* zero_pos=(long*)malloc((N1+1)*sizeof(long));  
     // Find dimension and fill in the positions
     for (long i=1;i<=N1;i++)
	  if (fabs(A(i,i))<VERY_TINY) 
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
	  //Vector V(N2);
	  K(zero_pos[i],i)=1.0; // Make up
	  for (long j=N2;j>=1;j--)
	       if (fabs(A(j,j))>VERY_TINY)
	       {
		    double sum=0.0;
		    for (long k=N1;k>j;k--)
			 sum+=A(j,k)*K(k,i);
		    K(j,i)=-sum/A(j,j);
	       }
     }    
     free(zero_pos);
     return dim;
}

long Matrix::Eigenvectors(Matrix &V, double E) const
{
     Matrix A(*this);
     for (long i=1;i<=N1;i++)
	  A(i,i)-=E;
     return A.Kernel(V);
}


// All I/O routines return 1 if ok, 0 if there was any error
int Matrix::Save_Binary(FILE *fich) const
{
     int ausgang;
     ausgang=fwrite(&N1,sizeof(long),1,fich);
     if (ausgang!=1) return 0;
     ausgang=fwrite(&N2,sizeof(long),1,fich);
     if (ausgang!=1) return 0;
     if (!(N1*N2)) return 1;
     ausgang=fwrite(D,sizeof(double),N1*N2+1,fich);
     if (ausgang!=N1*N2+1) return 0;
     return 1;
}

int Matrix::Save_Binary(const char *s) const
{
     FILE *fich=fopen(s,"wb");
     if (!fich) return 0;
     int ausgang=Save_Binary(fich);
     fclose(fich);
     return ausgang;
}

int Matrix::Load_Binary(FILE *fich)
{
     int ausgang;
     ausgang=fread(&N1,sizeof(long),1,fich);
     if (ausgang!=1) return 0;
     ausgang=fread(&N2,sizeof(long),1,fich);
     if (ausgang!=1) return 0;
     if (!(N1*N2)) { D=(double*)NULL; return 1;}
     Create(N1,N2);
     ausgang=fread(D,sizeof(double),N1*N2+1,fich);
     if (ausgang!=N1*N2+1) return 0;
     return 1;
}

int Matrix::Load_Binary(const char *s) 
{
     FILE *fich=fopen(s,"rb");
     if (!fich) return 0;
     int ausgang=Load_Binary(fich);
     fclose(fich);
     return ausgang;
}

int Matrix::Save(const char *name) const
{
     FILE *fich=fopen(name,"wt");
     if (!fich) return 0;
     int ausgang=Save(fich);
     fclose(fich);
     return ausgang;
}

int Matrix::Load(const char *name)
{
     FILE *fich=fopen(name,"rt");
     if (!fich) return 0;
     int ausgang=Load(fich);
     fclose(fich);
     return ausgang;
}

int Matrix::Save(FILE *fich) const
{
     fprintf(fich,"%ld %ld\n",N1,N2);
     for (long i=1;i<=N1;i++)
     {
	  for (long j=1;j<=N2;j++)
	       fprintf(fich,"%14.10g ",Elem(i,j));
	  fprintf(fich,"\n");
     }
     if (ferror(fich)) return 0;
     return 1;
}

int Matrix::Load(FILE *fich)
{
     float x;
     long n1, n2;
     if (!fscanf(fich,"%ld %ld\n",&n1,&n2))
	  merror("Error in Matrix::Load.");
     Create(n1,n2);
     for (long i=1;i<=n1;i++)
	  for (long j=1;j<=n2;j++)
	  {
	       if (feof(fich)) return 0; 
	       if (!fscanf(fich,"%g",&x))
		    merror("Error in Matrix::Load.");
	       Elem(i,j)=x;
	  }
     if (ferror(fich)) return 0;
     return 1;
}

////////////////////////////////////////////////////////////

// Full diagonalization: all eigenvalues, all eigenvectors
// Symmetric matrix.
// NOT OPTIMAL: matrix is preserved 
void Matrix::Diagonalize(Matrix &B, Vector &E) const
{
#ifdef DEBUG
     if (N1!=N2) merror("Trying to diagonalize non-square matrix.");
#endif
     char jobz='V', range='A', uplo='U';
     long N=N1;
     double *A=(double*)malloc(N*N*sizeof(double));
     memcpy(A,D+1,N*N*sizeof(double));
     long lda=N;
     double vl, vu;
     long il, iu;
     char askmachine='S';
     double abstol=2.0*dlamch_(&askmachine);
     long m; // total number of eigenvalues found
     E.Create(N); double *W=E.D+1;
     B.Create(N); double *Z=B.D+1;
     long ldz=N;
     double *work=(double*)malloc(8*N*sizeof(double));
     long lwork=8*N;
     long *iwork=(long*)malloc(5*N*sizeof(long));
     long *ifail=(long*)malloc(N*sizeof(long));
     long info;
     dsyevx_(&jobz,&range,&uplo,&N,A,&lda,&vl,&vu,&il,&iu,&abstol,
	     &m,W,Z,&ldz,work,&lwork,iwork,ifail,&info);
//     if (info) merror("Error during diagonalization.");
     free(iwork);
     free(ifail);
     free(work);
     free(A);
}

void Matrix::Spectrum(Vector &E) const
{
#ifdef DEBUG
     if (N1!=N2) merror("Trying to diagonalize non-square matrix.");
#endif
     char jobz='N', range='A', uplo='U';
     long N=N1;
     double *A=(double*)malloc(N*N*sizeof(double));
     memcpy(A,D+1,N*N*sizeof(double));
     long lda=N;
     double vl, vu;
     long il, iu;
     char askmachine='S';
     double abstol=2.0*dlamch_(&askmachine);
     long m; // total number of eigenvalues found
     E.Create(N); double *W=E.D+1;
     //B.Create(N); double *Z=B.D+1;
     double *Z=NULL;
     long ldz=N;
     double *work=(double*)malloc(8*N*sizeof(double));
     long lwork=8*N;
     long *iwork=(long*)malloc(5*N*sizeof(long));
     long *ifail=(long*)malloc(N*sizeof(long));
     long info;
     dsyevx_(&jobz,&range,&uplo,&N,A,&lda,&vl,&vu,&il,&iu,&abstol,
	     &m,W,Z,&ldz,work,&lwork,iwork,ifail,&info);
     free(A); free(work); free(iwork); free(ifail);
//     if (info) merror("Error during diagonalization.");
}

void Matrix::Tridiagonalize(Matrix &B, Vector &Diag, Vector &S) const
{
     // Needs dsytrd y dorgtr
     // Reduce to tri-diagonal form, basis is in "strange form"
     char uplo='U';
     long N=N1;
     double *A=(double*)malloc((N*N+1)*sizeof(double));
     memcpy(A+1,D+1,N*N*sizeof(double));
     long lda=N;
     Diag.Create(N); double* d=Diag.D+1;
     S.Create(N-1); double* e=S.D+1;
     double* tau=(double*)malloc(N*sizeof(double));
     long lwork=N*N;
     double* work=(double*)malloc(lwork*sizeof(double));
     long info;
     dsytrd_(&uplo, &N, A+1, &lda, d, e,
	     tau, work, &lwork, &info);
     if (info) merror("Error during tridiagonalization, stage 1.");
     // Compute the basis in "normal form" for tri-diagonal reduction
     dorgtr_(&uplo, &N, A+1, &lda, tau, work, &lwork, &info);
//     if (info) merror("Error during tridiagonalization, stage 2.");
     B.Create(N);
     B.Load(A,N,N);
     free(work); free(tau);

}

void Matrix::Ns_Diagonalize(Matrix &BL, Matrix &BR, 
			    Vector &ER, Vector &EI) const
{
#ifdef DEBUG
     if (N1!=N2) merror("Can't Ns_Diagonalize non-square matrix.");
#endif
     char jobl='V';
     char jobr='V';
     long N=N1;
     double *A=(double*)malloc(N*N*sizeof(double));
     memcpy(A,D+1,N*N*sizeof(double));
     long lda=N;
     ER.Create(N); EI.Create(N);
     BL.Create(N); BR.Create(N);
     double* wr=ER.D+1;
     double* wi=EI.D+1;
     double* vl=BL.D+1;
     long ldvl=N;
     double* vr=BR.D+1;
     long ldvr=N;
     long lwork=6*N;
     double* work=(double*)malloc(lwork*sizeof(double));
     long info;
 
     dgeev_(&jobl,&jobr,&N,A,&lda,
	    wr, wi, vl, &ldvl, vr, &ldvr, 
	    work, &lwork, &info);
     free(A); free(work); 
//     if (info) merror("Error diagonalizing non-symmetric matrix.");
}


/////////////////////////////////////////////////////
// EXTERNAL MATRICIAL FUNCTIONS

void Copy(Matrix& B, const Matrix& A) // B <- A raw and strict copy. 
{
     B.Destroy();
     if (!A.N1) return;
     B.Start(A.N1,A.N2);
     memcpy(B.D+1,A.D+1,A.N1*A.N2*sizeof(double));
}

Matrix operator-(const Matrix &M)
{
     Matrix R(M);
     R*=-1.0;
     return R;
}

Matrix operator+(const Matrix &A, const Matrix &B)
{
#ifdef DEBUG
     if (A.N1!=B.N1 || A.N2!=B.N2) 
	  merror("Adding matrices with different dimensions.");
#endif
     Matrix R(A);
     R+=B;
     return R;
}

// It is correct like this, although it doesn't seem to be!!
// This thing takes its arguments opposite: B-A!!
Matrix operator-(const Matrix &A, const Matrix &B)
{
#ifdef DEBUG
     if (A.N1!=B.N1 || A.N2!=B.N2) 
	  merror("Substracting matrices with different dimensions.");
#endif
     Matrix R(A);
     R-=B;
     return R;
}

Matrix operator*(double K, Matrix A)
{
     Matrix R(A);
     R*=K;
     return R;
}

Matrix operator*(Matrix A, double K)
{
     Matrix R(A);
     return R*K;
}

Vector operator*(const Matrix &A, const Vector &V)
{
#ifdef DEBUG
     if (A.N2!=V.N) merror("Wrong dims in Matrix-Vector product.");
#endif
     Vector R;
     Multiply(R,A,V);
     return R;
}

Matrix operator*(const Matrix &A, const Matrix &B)
{
#ifdef DEBUG
     if (A.N2!=B.N1) merror("Wrong dims in Matrix-Matrix product.");
#endif
     Matrix R;
     Multiply(R,A,B);
     return R;
}

void Multiply(Vector &R, const Matrix &M, const Vector &V)
{
     double alpha=1.0, beta=0.0;
     long incx=1, n1=M.N1, n2=M.N2;
     char c='N';
     R.Create(n1);
     dgemv_(&c, &n1, &n2, &alpha, M.D+1, &n1, V.D+1, &incx, 
	    &beta, R.D+1, &incx);
}
void Multiply(Matrix &R, const Matrix &M1, const Matrix &M2)
{
     double alpha=1.0, beta=0.0;
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
     double *D2=(double*)malloc((n1*n2+1)*sizeof(double));
     dgemm_(&ca,&cb,&n1,&n2,&k,&alpha,M1.D+1,&n1,M2.D+1,&k,&beta,D2+1,&n1);
     R.Load(D2,n1,n2);     
}

// The most general Matrix-Matrix product routine
// R <- alpha*M1*M2 + beta*R
// M1 is transposed if T1=true
// M2 is transposed if T2=true
void Multiply_Add(Matrix &R, const Matrix &M1, const Matrix &M2,
		  double alpha, double beta, bool T1, bool T2)
{
     long n1=(T1 ? M1.N2 : M1.N1), 
	  n2=(T2 ? M2.N1 : M2.N2), 
	  k=(T1 ? M1.N1 : M1.N2);
#ifdef DEBUG
     long kk=(T2 ? M2.N2 : M2.N1);
     if (k!=kk) merror("Wrong dimensions in Multiply_Add");
     if (M1.N1==0 || M1.N2==0 || M2.N1==0 || M2.N2==0)
	  merror("Zero matrix in Multiply_Add\n");
#endif
     char ca=(T1==true ? 'T' : 'N');
     char cb=(T2==true ? 'T' : 'N');
     long N1=M1.N1;
     long N2=M2.N1;
     dgemm_(&ca,&cb,&n1,&n2,&k,&alpha,M1.D+1,&N1,M2.D+1,&N2,&beta,R.D+1,&n1);
}

// Diagonalize a tridiagonal matrix (diagonal D, subdiagonal S)
// and return the spectrum in D
void Td_Spectrum(Vector &D, Vector &S)
{
     char compz='N';
     long N=D.N;
     double *z=(double*)NULL;
     long ldz=N;
     double *work=(double*)NULL;
     long info;
     dsteqr_(&compz,&N,D.D+1,S.D+1,z,&ldz,work,&info);
}

// B should be already initialized, a Unit matrix
void Td_Diagonalize(Matrix &B, Vector &D, Vector &S)
{
     char compz='V';
     long N=D.N;
     double *z=B.D+1;
     long ldz=N;
     double *work=(double*)malloc(2*(N-1)*sizeof(double)); 
     long info;
     dsteqr_(&compz,&N,D.D+1,S.D+1,z,&ldz,work,&info);
     free(work);
//     if (info) merror("Error during Td_Diagonalize.");
}

Matrix Laplacian_1D(long N, int b)
{
     Matrix L(N);
     for (long i=1;i<=N;i++)
     {
         if (i<N) L(i,i+1)=L(i+1,i)=-1.0;
         L(i,i)=2.0;
     }
     if (b) L(1,1)=L(N,N)=2.0;
     return L;
}

Matrix Tens_Prod(const Matrix &A, const Matrix &B)
{
     Matrix R;
     Tens_Prod(R,A,B);
     return R;
}

void   Tens_Prod(Matrix &R, const Matrix &A, const Matrix &B)
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

Matrix Tens_Prod_Unit(const Matrix & A, long m)
{
     Matrix R;
     Tens_Prod_Unit(R,A,m);
     return R;
}

void   Tens_Prod_Unit(Matrix& R, const Matrix& A, long m)
{
     long n1=A.N1, n2=A.N2;
     R.Create(n1*m,n2*m);
     R.Zero();
     for (long i1=1;i1<=n1;i1++)
	  for (long i2=1;i2<=n2;i2++)
	       for (long j=1;j<=m;j++)
		    R((i1-1)*m+j,(i2-1)*m+j)=A(i1,i2);
}

Matrix Tens_Prod_Unit(long n, const Matrix &A)
{
     Matrix R;
     Tens_Prod_Unit(R,n,A);
     return R;
}

void   Tens_Prod_Unit(Matrix& R, long n, const Matrix& A)
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
void   Tens_Prod_Unit(Matrix &R, const Matrix &A, long n, Side s)
{
     if (s==Left)
	  Tens_Prod_Unit(R,n,A);
     else
	  Tens_Prod_Unit(R,A,n);
}

Matrix Tens_Prod_Unit(const Matrix &A, long n, Side s)
{
     if (s==Left)
	  return Tens_Prod_Unit(n,A);
     else
	  return Tens_Prod_Unit(A,n);
}

void Elem_Mult(Matrix &R, const Matrix &A, const Matrix &B)
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

Matrix Elem_Mult(const Matrix &A, const Matrix &B)
{
     Matrix R(A.N1,A.N2);
     Elem_Mult(R,A,B);
     return R;
}

void Ket_Bra(Matrix &R, const Vector &V, const Vector &W)
{
     long N1=V.N, N2=W.N;
     R.Create(N1,N2);
     for (long i=1;i<=N1;i++)
	  for (long j=1;j<=N2;j++)
	       R(i,j)=V(i)*W(j);
}

Matrix Unit(long N1, long N2)
{
     Matrix U(N1,N2);
     U.Unit();
     return U;
}

Matrix Zero(long N1, long N2)
{
     Matrix Z(N1,N2);
     return Z;
}

Matrix Diag(const Vector &E)
{
     long N=E.N;
     Matrix R(N);
     for (long i=1;i<=N;i++)
	  R(i,i)=E(i);
     return R;
}

// Order V and swap columns of M accordingly
// Hm... don't quite like the implementation...
void Order(Vector &V, Matrix &M)
{
     // Uufff... we should make a documentation for this
     // Bubble algorithm, I'm afraid...
     double acum;
     long nswaps;
     do
     {
	  nswaps=0;
	  for (long i=1;i<M.N2;i++)
	       if (V(i)>V(i+1))
	       {
		    acum=V(i);
		    V(i)=V(i+1);
		    V(i+1)=acum;
		    M.Swap_Cols(i,i+1);
		    nswaps++;
	       }
     }while(nswaps);
}


double Elem(const Matrix &M, long i, long j)
{
     return M.Elem(i,j);
}

double Elem(const Matrix &M, const Vector &V1, const Vector &V2)
{
     return M.Elem(V1,V2);
}

Vector Col(const Matrix &M, long i)
{
     return M.Col(i);
}

Vector Row(const Matrix &M, long i)
{
     return M.Row(i);
}

// Matrix Change_Basis(const Matrix &M, const Matrix &B)
// {
//      Matrix A(M);
//      A.Change_Basis(B);
//      return A;
// }

// Matrix Gram_Schmidt(const Matrix &M)
// {
//      Matrix A(M);
//      A.Gram_Schmidt();
//      return A;
// }

// Matrix T(const Matrix &M)
// {
//      Matrix A(M);
//      A.T();
//      return A;
// }

Matrix Insert(const Matrix &M, const Matrix &Q, long i, long j)
{
     Matrix A(M);
     A.Insert(Q,i,j);
     return A;
}

Matrix Invert(const Matrix &M)
{
     Matrix A(M);
     A.Invert();
     return A;
}

double Trace(const Matrix &M)
{
     return M.Trace();
}

double Det(const Matrix &M)
{
     return M.Det();
}

double Frobenius_Norm(const Matrix &M)
{
#ifdef DEBUG
     if (!M.N1 || !M.N2) merror("Frobenius norm of non-existent matrix\n");
#endif
     double norm=0.0;
     for (long i=1;i<=M.N1;i++)
	  for( long j=1;j<=M.N2;j++)
	       norm+=sqr(M(i,j));
     return sqrt(norm);
}

#endif
