////////////////////////////////////////////////////////
// cmatrix, complex vector and matrix classes
// 080912-101022
// JaviRL
// Strongly based on matrix, from hvb++

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "matrix.h"

#ifndef CMATRIX_HEADER
#define CMATRIX_HEADER

// If you want to trace your errors easily, write this in your program
// #define DEBUG

///////////////////////////////////////////////////////////////////////////
/// HEADERS:
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
// Class cvector:
///////////////////////////////////////////////////////////////////////////

class CVector
{
public:
     long N; // length of the vector 
     cmplx *D; // the pointer itself
     
     CVector(long);
     CVector(const Vector &);
     CVector(const CVector &);
     CVector(cmplx *, long);
     CVector();
     ~CVector();
     void Start(long n=0);
     void Create(long n);
     void Create();
     void Load(cmplx *, long);
     void Load_Copy(cmplx *, long);
     void Transfer(CVector &);
     void Destroy();
     void Zero();
     void Set_Value(cmplx);
     bool Is_Zero(double tolerance);
     double Norm() const;
     void Write() const;
     void Write_Col() const;
     void Write(FILE *f) const;
     void Write_Col(FILE *f) const;
     void Write(int prec) const;
     void Write(int prec, FILE *f) const;
     void Append(const cmplx); 
     void Append(const CVector &);
     void Insert(const cmplx, long i);
     void Remove(long i);
     cmplx Min() const;
     cmplx Min(long) const;
     cmplx Max() const;
     cmplx Max(long) const;
     long Min_Index() const;
     long Min_Index(long) const;
     long Max_Index() const;
     long Max_Index(long) const;
     cmplx Sum() const;
     cmplx Sum(long i1, long i2) const;
     cmplx Average() const;
     cmplx Deviation() const;
     cmplx Variance() const;
     void Order(int p=1);
     cmplx operator() (long) const;
     cmplx& operator() (long);
     cmplx Last() const;
     CVector& operator=(const CVector&);
     CVector& operator=(const Vector&);
     void operator+=(const CVector&);
     void operator-=(const CVector&);
     void operator*=(const cmplx);
     void operator/=(const cmplx);
     void operator*=(const double);
     void operator/=(const double);

     bool Normalize();
     void Part(long, long);
     void Reverse();
     void Insert(const CVector&, long);
     void Sqr();
     void Real();
     void Imag();
     void Abs();
     void Conj();
 
     /* CVector operator+(const CVector&) const;  */
     /* CVector operator-(const CVector&) const; */
     /* CVector operator-() const; */

     int Save_Binary(const char *name) const;
     int Load_Binary(const char *name);
     int Save_Binary(FILE *fich) const;
     int Load_Binary(FILE *fich);
     int Save(const char *name) const;
     int Load(const char *name);
     int Save(FILE *fich) const;
     int Load(FILE *fich);
}; 

void Copy(CVector& B, const CVector& A);
void Copy(CVector& B, const Vector& A);


CVector Normalize(const CVector &V);
CVector Part(const CVector &V,long,long);
CVector Reverse(const CVector &V);
CVector Insert(const CVector &V, const CVector &W, long);
CVector Sqr(const CVector &V); 
CVector Abs(const CVector &V);
CVector Real(const CVector &V);
CVector Imag(const CVector &V);
CVector Conj(const CVector &V);
Vector To_Real(const CVector &V);

CVector operator+(const CVector &V, const CVector &W);
CVector operator-(const CVector &V);
CVector operator*(cmplx, const CVector &);
CVector operator*(const CVector &, cmplx);
CVector operator*(double, const CVector &);
CVector operator*(const CVector &, double);
CVector operator-(const CVector &, const CVector &);
CVector operator/(const CVector &, cmplx);
cmplx Dot(const CVector &, const CVector &);
CVector Cross(const CVector &, const CVector &);
CVector Tens_Prod(const CVector &, const CVector &);
void Tens_Prod(CVector &, const CVector &, const CVector &);
void Daxpy(CVector &, const CVector &, cmplx);
CVector Elem_Mult(const CVector &, const CVector &);
void Elem_Mult(CVector &, const CVector &, const CVector &);

/////////////////////////////////////////////////////////////////////////
// Class matrix:
/////////////////////////////////////////////////////////////////////////

class CMatrix
{
public:
     cmplx *D;
     long N1,N2;

     CMatrix(long n, long m=0);
     CMatrix(const CMatrix &);
     CMatrix(const Matrix &);
     CMatrix();
     ~CMatrix();
     void Start(long n, long=0);
     void Create(long n, long=0);
     void Create();
     void Load(cmplx*, long, long=0);
     void Load_Copy(cmplx*, long, long=0);
     void Transfer(CMatrix &M);
     void Destroy();
     void Zero();
     bool Is_Zero(double tolerance);
     void Unit();
     void Set_Value(cmplx);
     void Resize(long,long=0);
     cmplx& Elem(long, long);
     cmplx Elem(long, long) const;
     cmplx& operator()(long, long);
     cmplx operator() (long, long) const;
     void Write() const;
     void Write(FILE *f) const;
     void Write(int prec) const;
     void Write(int prec, FILE *f) const;
     CMatrix& operator=(const CMatrix&);
     CMatrix& operator=(const Matrix&);
     void operator+=(const CMatrix&);
     void operator-=(const CMatrix&);
     void operator*=(cmplx);  
     void operator*=(double);
     void operator*=(const CMatrix&);
     CVector Col(long n) const;
     CVector Row(long n) const;
     void Col(CVector &, long) const;
     void Put_Col(const CVector &,long);
     void Row(CVector &, long) const;
     void Put_Row(const CVector &, long);      
     void Append_Col(const CVector &);
     void Append_Row(const CVector &);
     void Swap_Cols(long,long);
     void Swap_Rows(long,long);
     cmplx Elem(const CVector &, const CVector &) const;
     cmplx Elem(const CMatrix &, long, long) const;
     void Sort_Cols();

     void Change_Basis(const CMatrix &);
     bool Gram_Schmidt();
     void T();
     void Herm(); // hermitian conjugate
     void Part(long, long, long, long);
     void Real();
     void Imag();
     void Abs();

     void Add(const CMatrix &, long, long);
     void Insert(const CMatrix &, long, long);

     void Solve(CVector &) const;
     void Solve(CMatrix &) const;
     void Invert();
     int LU_Decomp(int *);          

     cmplx Det() const; // determinant
     cmplx Trace() const; // trace
     long Kernel(CMatrix &) const;
     long Eigenvectors(CMatrix &, cmplx) const;
     
     void Diagonalize(CMatrix &,Vector &) const; 
     void Spectrum(Vector &) const;
     void Tridiagonalize(CMatrix &, CVector &, CVector &) const;
     
     void Ns_Diagonalize(CMatrix&, CVector &) const; 
     void SVD(CMatrix &, CMatrix &, Vector &) const;

     int Save_Binary(const char *name) const;
     int Load_Binary(const char *name);
     int Save_Binary(FILE *fich) const;
     int Load_Binary(FILE *fich);
     int Save(const char *name) const;
     int Load(const char *name);
     int Save(FILE *fich) const;
     int Load(FILE *fich);
};

void Copy(CMatrix& B, const CMatrix& A);
void Copy(CMatrix& B, const Matrix& A); 
//CMatrix Unit(long,long=0);
//CMatrix Zero(long,long=0);

CMatrix operator-(const CMatrix &);
CMatrix operator+(const CMatrix &, const CMatrix &);
CMatrix operator-(const CMatrix &, const CMatrix &);
CMatrix operator*(cmplx K, const CMatrix &A);
CMatrix operator*(double K, const CMatrix &A);
CVector operator*(const CMatrix &A, const CVector &V);
CMatrix operator*(const CMatrix &A, const CMatrix &B);

void Multiply(CVector &, const CMatrix &M, const CVector &V);
void Multiply(CMatrix &, const CMatrix &M1, const CMatrix &M2);
void Multiply_Add(CMatrix &, const CMatrix &M1, const CMatrix &M2,
		  cmplx alpha, cmplx beta, bool T1, bool T2);

CMatrix Tens_Prod(const CMatrix &, const CMatrix &);
void   Tens_Prod(CMatrix &, const CMatrix &, const CMatrix &);
CMatrix Tens_Prod_Unit(const CMatrix &, long);
void   Tens_Prod_Unit(CMatrix &, const CMatrix&, long);
CMatrix Tens_Prod_Unit(long, const CMatrix &);
void   Tens_Prod_Unit(CMatrix &, long, const CMatrix&);

void   Tens_Prod_Unit(CMatrix &, const CMatrix &, long, Side);
CMatrix Tens_Prod_Unit(const CMatrix &, long, Side);

void Elem_Mult(CMatrix &, const CMatrix &, const CMatrix &);
CMatrix Elem_Mult(const CMatrix &, const CMatrix &);

void Ket_Bra(CMatrix &, const CVector &, const CVector &);

CMatrix Diag(const CVector &E);

//CMatrix Laplacian_1D(long N, int b=0);

//void Td_Spectrum(CVector &D, CVector &S);
void Td_Diagonalize(CMatrix &B, CVector &D, CVector &S);

// Order V and swap columns of M accordingly
void Order(CVector &V, CMatrix &M);

CMatrix Change_Basis(const CMatrix &M, const CMatrix &B);
CMatrix T(const CMatrix &M);
CMatrix Herm(const CMatrix &M);
CMatrix Part(const CMatrix &M, long, long, long, long);
void Part(CMatrix &P, const CMatrix &M, long, long, long, long);
CMatrix Invert(const CMatrix &M);
CVector Solve(const CMatrix &A, const CVector &b);
cmplx Det(const CMatrix &M);
cmplx Trace(const CMatrix &M);
double Frobenius_Norm(const CMatrix &M);

CMatrix Real(const CMatrix &M);
CMatrix Imag(const CMatrix &M);
CMatrix Abs(const CMatrix &M);
Matrix To_Real(const CMatrix &M);
CMatrix To_Cmplx(const Matrix &M);


/////////////////////////////////////////////////////////////
/// BLAS-LAPACK HEADERS
extern "C"{
     // V <- alpha *V
     void zscal_(long *N, cmplx *alpha, cmplx* V, long *ix);
     // V <- V + W
     void zaxpy_(long *N, cmplx *alpha, cmplx *X, long *ix, cmplx *Y, long *iy);
     // Dot product
     cmplx zdotc_(long *N, cmplx *X, long *ix, cmplx *Y, long *iy);
     // Matrix-vector product
     void zgemv_(char *, long* n1, long* n2, cmplx* alpha, cmplx* A, long *lda,
		 cmplx* X, long *incx, cmplx* beta, cmplx* Y, long* incy);
     // Matrix-matrix product
     void zgemm_(char*, char*, long* n1, long* n2, long*k, cmplx* alpha, 
		 cmplx* A, long* lda, cmplx* B, long *ldb,
		 cmplx* beta, cmplx* C, long *ldc);
     // LU decomposition
     void zgetrf_(long*, long*, cmplx*, long*, int*, long*);
     // Linear equations solving
     void zgesv_(long*, long*, cmplx*, long*, long*, cmplx*, long*, long*);
     // Tridiagonal Matrix diagonalization
     void zsteqr_(char* compz, long* n, double* d, double* e,
		  cmplx* z, long* ldz, cmplx* work,long* info);
     // Full Diagonalization, expert driver
     void zheevx_(char* jobz, char* range, char* uplo, long* n, cmplx* A,
		  long* lda, double* vl, double* vu,
		  long* il, long* iu, double* abstol, long* M, double* W,
		  cmplx* Z, long* ldz, cmplx* work,
		  long* lwork, double* rwork, long* iwork,
		  long* ifail, long* info);
     // Non-symmetric full diagonalization, non-expert driver
     void zgeev_(char* jobl, char* jobr, long* N, cmplx* A, long* lda,
		 cmplx* w, cmplx* vl, long* ldvl, cmplx* vr, long* ldvr, 
		 cmplx* work, long* lwork, double* rwork, long* info);
     // Reduce to tri-diagonal form, basis is in "strange form"
     void zhetrd_(char* uplo, long* N, cmplx* A, long* lda, cmplx* D, cmplx* E,
		  cmplx* tau, cmplx* work, long* lwork, long* info);
     // Compute the basis in "normal form" for tri-diagonal reduction
     void zungtr_(char* uplo, long* N, cmplx* A, long* lda, cmplx* tau,
		  cmplx* work, long* lwork, long* info);
     // SVD
     void zgesvd_(char* jobu, char* jobvt, long *M, long *N, cmplx *A,
		  long *lda, double *S, cmplx *U, long *ldu, cmplx *Vt,
		  long *ldvt, cmplx *work, long *lwork, double *rwork,
		  long *info);
}

#endif



