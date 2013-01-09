////////////////////////////////////////////////////////
// hvb++ 1.0
// Copyleft: Javier Rodríguez Laguna
// 080725
// 080912, we return to double, instead of real
// 101022, simplify interface

//////////////////////////////////////////////////////////////////////////
// C++ library for linear algebra, JaviRL, started Nov 1999.
//
// Version 2: Aug 2003, split into modules.
// Version 3: May 2005, long ints, some functions added
// Version 3: Sep 2005, rationalize interface, link to blas and lapack
// Minor update,  060224
// Adding Format, 060822
// Version 3.1: Oct 2006
// Version 034: Mar 2008
//////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"

#ifndef MATRIX_HEADER
#define MATRIX_HEADER

// If you want to trace your errors easily, write this in your program
// #define DEBUG

///////////////////////////////////////////////////////////////////////////
/// HEADERS:
///////////////////////////////////////////////////////////////////////////

#ifdef DEBUG
extern long num_vec;
extern long num_mat;
extern long total_memory;
extern long max_total_memory;
extern bool matrix_error;
void Mem_Control();
#endif

///////////////////////////////////////////////////////////////////////////
// Class vector:
///////////////////////////////////////////////////////////////////////////

class Vector
{
public:
     long N; // length of the vector 
     double *D; // the pointer itself
     
     Vector(long);
     Vector(const Vector &);
     Vector(double *, long);
     Vector();
     ~Vector();
     void Start(long n=0);
     void Create(long n);
     void Load(double *, long);
     void Load_Copy(double *, long);
     void Transfer(Vector &);
     void Destroy();
     void Zero();
     void Set_Value(double);
     bool Is_Zero(double tolerance);
     double Norm() const;
     void Write_Col() const;
     void Write_Col(FILE *f) const;
     void Write() const;
     void Write(FILE *f) const;
     void Write(int prec) const;
     void Write(int prec, FILE *f) const;
     void Append(const double); 
     void Append(const Vector &);
     void Insert(const double, long i);
     void Remove(long i);
     double Min() const;
     double Min(long) const;
     double Max() const;
     double Max(long) const;
     long Min_Index() const;
     long Min_Index(long) const;
     long Max_Index() const;
     long Max_Index(long) const;
     double Sum() const;
     double Sum(long i1, long i2) const;
     double Average() const;
     double Deviation() const;
     double Variance() const;
     void Order(int p=1);
     double operator() (long) const;
     double& operator() (long);
     double Last() const;
     Vector& operator=(const Vector&);
     void operator+=(const Vector&);
     void operator-=(const Vector&);
     void operator*=(const double);
     void operator/=(const double);

     bool Normalize();
     void Part(long, long);
     void Reverse();
     void Insert(const Vector&, long);
     void Sqr();
 
     int Save_Binary(const char *name) const;
     int Load_Binary(const char *name);
     int Save_Binary(FILE *fich) const;
     int Load_Binary(FILE *fich);
     int Save(const char *name) const;
     int Load(const char *name);
     int Save(FILE *fich) const;
     int Load(FILE *fich);
}; 

void Copy(Vector& B, const Vector& A);
Vector operator-(const Vector &);
Vector operator+(const Vector &, const Vector &);
Vector operator-(const Vector &, const Vector &);
Vector operator*(double, const Vector &);
Vector operator*(const Vector &, double);
Vector operator/(const Vector &, double);
double Dot(const Vector &, const Vector &);
Vector Cross(const Vector &, const Vector &);
Vector Tens_Prod(const Vector &, const Vector &);
void Tens_Prod(Vector &, const Vector &, const Vector &);
void Daxpy(Vector &, const Vector &, double);
Vector Elem_Mult(const Vector &, const Vector &);
void Elem_Mult(Vector &, const Vector &, const Vector &);

Vector Normalize(const Vector &);
Vector Part(const Vector &,long, long);
Vector Reverse(const Vector &);
Vector Insert(const Vector &,const Vector&, long);
Vector Sqr(const Vector &);

double Norm(const Vector &);
double Min(const Vector &);
double Max(const Vector &);
double Sum(const Vector &);
double Average(const Vector &);
double Deviation(const Vector &);
double Variance(const Vector &);


/////////////////////////////////////////////////////////////////////////
// Class matrix:
/////////////////////////////////////////////////////////////////////////

class Matrix
{
public:
     double *D;
     long N1,N2;

     Matrix(long n, long m=0);
     Matrix(const Matrix &);
     Matrix();
     ~Matrix();
     void Start(long n, long=0);
     void Create(long n, long=0);
     void Load(double *, long, long=0);
     void Load_Copy(double *, long, long=0);
     void Transfer(Matrix &M);
     void Destroy();
     void Zero();
     bool Is_Zero(double tolerance);
     void Unit();
     void Set_Value(double);
     void Resize(long,long=0);
     double& Elem(long, long);
     double Elem(long, long) const;
     double& operator()(long, long);
     double operator() (long, long) const;
     void Write() const;
     void Write_Col() const;
     void Write(FILE *f) const;
     void Write_Col(FILE *f) const;
     Matrix& operator=(const Matrix&);
     void operator+=(const Matrix&);
     void operator-=(const Matrix&);
     void operator*=(double);  
     void operator*=(const Matrix&);
     Vector Col(long n) const;
     Vector Row(long n) const;
     void Col(Vector &, long) const;
     void Put_Col(const Vector &,long);
     void Row(Vector &, long) const;
     void Put_Row(const Vector &, long);      
     void Append_Col(const Vector &);
     void Append_Row(const Vector &);
     void Swap_Cols(long,long);
     void Swap_Rows(long,long);
     double Elem(const Vector &, const Vector &) const;
     double Elem(const Matrix &, long, long) const;
     void Sort_Cols();

     void Change_Basis(const Matrix &);
     void Gram_Schmidt();
     void T();
     void Part(long, long, long, long);
     void Add(const Matrix &, long, long);
     void Insert(const Matrix &, long, long);
     void Invert();

     void Solve(Vector &) const;
     void Solve(Matrix &) const;

     int LU_Decomp(int *);          
     double Det() const; // determinant
     double Trace() const; // trace
     long Kernel(Matrix &) const;
     long Eigenvectors(Matrix &, double) const;
     
     void Diagonalize(Matrix &,Vector &) const; 
     void Spectrum(Vector &) const;
     void Tridiagonalize(Matrix &, Vector &, Vector &) const;
     
     void Ns_Diagonalize(Matrix&, Matrix&, Vector &, Vector &) const; 

     int Save_Binary(const char *name) const;
     int Load_Binary(const char *name);
     int Save_Binary(FILE *fich) const;
     int Load_Binary(FILE *fich);
     int Save(const char *name) const;
     int Load(const char *name);
     int Save(FILE *fich) const;
     int Load(FILE *fich);
};

void Copy(Matrix& B, const Matrix& A);
Matrix Unit(long,long=0);
Matrix Zero(long,long=0);

Matrix operator-(const Matrix &);
Matrix operator+(const Matrix &, const Matrix &);
Matrix operator-(const Matrix &, const Matrix &);
Matrix operator*(double K, Matrix A);
Vector operator*(const Matrix &A, const Vector &V);
Matrix operator*(const Matrix &A, const Matrix &B);

void Multiply(Vector &, const Matrix &M, const Vector &V);
void Multiply(Matrix &, const Matrix &M1, const Matrix &M2);
void Multiply_Add(Matrix &, const Matrix &M1, const Matrix &M2,
		  double alpha, double beta, bool T1, bool T2);

Matrix Tens_Prod(const Matrix &, const Matrix &);
void   Tens_Prod(Matrix &, const Matrix &, const Matrix &);
Matrix Tens_Prod_Unit(const Matrix &, long);
void   Tens_Prod_Unit(Matrix &, const Matrix&, long);
Matrix Tens_Prod_Unit(long, const Matrix &);
void   Tens_Prod_Unit(Matrix &, long, const Matrix&);

void   Tens_Prod_Unit(Matrix &, const Matrix &, long, Side);
Matrix Tens_Prod_Unit(const Matrix &, long, Side);

void Elem_Mult(Matrix &, const Matrix &, const Matrix &);
Matrix Elem_Mult(const Matrix &, const Matrix &);

void Ket_Bra(Matrix &, const Vector &, const Vector &);

Matrix Diag(const Vector &E);

Matrix Laplacian_1D(long N, int b=0);

void Td_Spectrum(Vector &D, Vector &S);
void Td_Diagonalize(Matrix &B, Vector &D, Vector &S);

// Order V and swap columns of M accordingly
void Order(Vector &V, Matrix &M);


double Elem(const Matrix &M, long, long);
double Elem(const Matrix &M, const Vector &V1, const Vector &V2);
Vector Col(const Matrix &N, long);
Vector Row(const Matrix &M, long);
Matrix Change_Basis(const Matrix &M, const Matrix &B);
Matrix T(const Matrix &M);
Matrix Part(const Matrix &M, long,long,long,long);
void Part(Matrix &R,const Matrix &M,long,long,long,long);
Matrix Insert(const Matrix &M, long, long);
Matrix Invert(const Matrix &M);
Vector Solve(const Matrix &A, const Vector &b);
double Trace(const Matrix &M);
double Det(const Matrix &M);
double Frobenius_Norm(const Matrix &M);

/////////////////////////////////////////////////////////////
/// BLAS-LAPACK HEADERS
extern "C"{
     // obtain machine parameters
     double dlamch_(char *c);
     // V <- alpha *V
     void dscal_(long *N, double *alpha, double* V, long *ix);
     // V <- V + W
     void daxpy_(long *N, double *alpha, double *X, long *ix, double *Y, 
		 long *iy);
     // Dot product
     double ddot_(long *N, double *X, long *ix, double *Y, long *iy);
     // Matrix-vector product
     void dgemv_(char *, long* n1, long* n2, double* alpha, double* A, 
		 long *lda, double* X, long *incx, double* beta, double* Y, 
		 long* incy);
     // Matrix-matrix product
     void dgemm_(char*, char*, long* n1, long* n2, long*k, double* alpha, 
		 double* A, long* lda, double* B, long *ldb,
		 double* beta, double* C, long *ldc);
     // LU decomposition
     void dgetrf_(long*, long*, double*, long*, int*, long*);
     // Linear equations solving
     void dgesv_(long*, long*, double*, long*, long*, double*, long*, long*);
     // Tridiagonal Matrix diagonalization
     void dsteqr_(char* compz, long* n, double* d, double* e,
		  double* z, long* ldz, double* work,long* info);
     // Full Diagonalization, expert driver
     void dsyevx_(char*,char*,char*,long*,double*,long*,double*,double*,
		  long*,long*,double*,long*,double*,double*,long*,double*,
		  long*,long*,long*,long*);
     // Non-symmetric full diagonalization, non-expert driver
     void dgeev_(char* jobl, char* jobr, long* N, double* A, long* lda,
		 double* wr, double* wi, double* vl, long* ldvl,
		 double* vr, long* ldvr, double* work, long* lwork, 
		 long* info);
     // Reduce to tri-diagonal form, basis is in "strange form"
     void dsytrd_(char* uplo, long* N, double* A, long* lda, double* D, 
		  double* E, double* tau, double* work, long* lwork, 
		  long* info);
     // Compute the basis in "normal form" for tri-diagonal reduction
     void dorgtr_(char* uplo, long* N, double* A, long* lda, double* tau,
		  double* work, long* lwork, long* info);
}

#endif



