////////////////////////////////////////////////////////
// hvb++ 1.0
// Copyleft: Javier Rodríguez Laguna
// 080725-110824

// Common routines library, for all (or almost all) my programs
// JaviRL, 060404

#ifndef COMMON_HEADER
#define COMMON_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>
#include <complex>
#include <iostream>
#include <ctype.h>

////////////////////////////////////////////////////////
// General definitions
using namespace std;
#define cmplx complex<double>
extern cmplx M_I;

double sqr(double a);

#define MIN(a,b)  ((a<b)?(a):(b))
#define MAX(a,b)  ((a>b)?(a):(b))
#define SIGN(a,b) ((b)>=0.0 ? fabs(a) : -fabs(a))

void SWAP(int &a, int &b);
void SWAP(long &a, long &b);
void SWAP(double &a, double &b);
void SWAP(cmplx &a, cmplx &b);

enum Side {Left, Right};

////////////////////////////////////////////////////////
// Error handling

void merror(const char *) __attribute__ ((noinline));

////////////////////////////////////////////////////////
// Time measurement

double uhr();
void delay(double Dt);

////////////////////////////////////////////////////////
// Random number generator

void open_rng(unsigned long s);
double rand_double();
double rand_double(double a, double b);
long rand_int(long a, long b);
void close_rng();

////////////////////////////////////////////////////////
// Binary numbers

long  pow2(int num);
bool  get_bit(long X, int pos);
int   get_max_bit(long X);
long  flip_bit(long X, int pos); 
long  put_bit(long X, int pos, int value);
long  swap_bits(long X, int pos1, int pos2);
long  reverse_bits(long X, int maxbit=0);
long  insert_bit(long X, int pos, int value);
int   count_ones(long X); 
long  next_in_sector(long X);
char* B2string(long X, int maxbit=0);
char* B2stringR(long X, int maxbit=0);
void  printB(long X);


////////////////////////////////////////////////////////
// Class list

class List
{
public:
     long N; // size of list
     long *D;

     List(); // forwarded initialization
     List(long); // initialize a list with a given size
     List(long*); // initialize a list copying from a given int*
     List(long*,long); // initialize a list copying a number of items from int*
     List(const List&); // initialize a list from another
     ~List(); // destructor
     void Create(); // void allocation
     void Create(long); // allocates
     void Start(); // initialize when you're sure that it is not created
     void Destroy(); // deallocates
     void Zero(); // puts all elements to zero

     long Find(long,long k=1) const; // finds starting from second argument
     long Find_Non(long,long k=1) const;
     long Find_Higher(long,long k=1) const;
     long Find_Lower(long,long k=1) const;

     void Remove_Index(long);
     long Remove_All(long); // eliminates all occurrences of (long)
     long Remove_Higher(long); // remove all values higher than long
     long Remove_Lower(long);  // remove all values lower than long

     void Append(long); // appends an element at the end of the list
     void Append(const List &); // appends a full list at the end of another

     long Count(long x) const; // counts the number of appearances of x
     long Min() const;
     long Max() const;
     long Max_Index() const;
     long Min_Index() const;
     void Order(long); // orders the list
     void Uniquify(); // removes all repeated values
     void Write() const; // writes list with standard pen
     void Write(FILE *) const; 
     void Insert(long, long); // insert a number at a given position
     void Insert(long, const List&); // insert a list at a given position
     void Reverse(); // reverse the order of a list
     void Swap(long,long); // swaps the values at two positions
     void Part(long,long); // returns a part of the list
     long operator() (long) const;
     long& operator() (long);
     List& operator=(const List&);
     void operator+=(const List&);
     void operator-=(const List&);
     void operator*=(const long);
     void operator/=(const long);
     void operator+=(const long);
     void operator-=(const long);

     int Save_Binary(FILE *) const;
     int Save_Binary(const char *name) const;
     int Load_Binary(FILE *);
     int Load_Binary(const char *name);
     int Load(const char *);
     int Load(FILE *);
     int Save(const char *) const;
     int Save(FILE *) const;

     long Sum() const;
     long Sum(long,long) const; // sum of elements between two given
     long Prod() const; // product of all elements
     long Prod(long,long) const; // takes the product of all elements

     void Set_Value(long); // initialize the list to a given value
};

List operator+(const List &L1, const List &L2); 
List operator-(const List &L1, const List &L2);
List operator+(const List &L, const long p); 
List operator-(const List &L, const long p);
List operator+(const long p, const List &L); 
List operator-(const long p, const List &L);
List operator-(const List &L);
List operator*(const List &L, const long p);
List operator/(const List &L, const long p);
List operator*(const long p, const List &L);

void Copy(List &L1, const List &L2);
bool Is_Equal(const List &L1, const List &L2);
List Intersection(const List &L1, const List &L2);
List List_Interval(long, long);
List Combination(const List &, const List &);
bool Is_Contained(const List &L1, const List &L2);
List Substraction(const List &L1, const List &L2);
// List Complement(const List &L1, long L);

List Reverse(const List &L);
List Part(const List &L, long i1, long i2);

List  B2List(long X, int maxbit=0); // conversion of a binary number to a list
long  List2B(const List &L);

// working with numbers in any numbering system
// X is the list of digits in that represeantation
// Ls is the maximum cardinality of each digit.
long To_Number(const List &X, const List &Ls);
List To_List(long ix, const List &Ls);

#endif













