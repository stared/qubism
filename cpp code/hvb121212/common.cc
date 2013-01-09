////////////////////////////////////////////////////////
// hvb++ 1.0
// Copyleft: Javier Rodríguez Laguna
// Common routines for hvb
// 060404-080725-110930

#include"common.h"

cmplx M_I(0.0,1.0);

void SWAP(int &a, int &b)
{
     int c; c=a; a=b; b=c;
}

void SWAP(long &a, long &b)
{
     long c; c=a; a=b; b=c;
}

void SWAP(double &a, double &b)
{
     double c; c=a; a=b; b=c;
}

void SWAP(cmplx &a, cmplx &b)
{
     cmplx c; c=a; a=b; b=c;
}


////////////////////////////////////////////////////////
// ERROR HANDLING
void merror(const char *s) 
{
     printf("%s\n",s);
     fflush(stdout);
     exit(1);
}

////////////////////////////////////////////////////////
// CPU TIME
double uhr()
{
    static const long clock_interval = sysconf(_SC_CLK_TCK);
    struct tms t;
    times(&t);
    return (t.tms_utime + t.tms_stime) * (1.0/clock_interval);
}

void delay(double Dt)
{
     double start=uhr();
     do {} while(uhr()-start<Dt);
}

////////////////////////////////////////////////////////
// RANDOM NUMBERS
// We'll use Mersenne Twister 19937 (See GSL-documentation, code is from there)
#define MTRNG_N 624
#define MTRNG_M 397
#define MTRNG_RAND_MAX 0xffffffffUL
static const unsigned long MTRNG_UPPER_MASK = 0x80000000UL;
static const unsigned long MTRNG_LOWER_MASK = 0x7fffffffUL;

typedef struct
{
     unsigned long mt[MTRNG_N];
     int mti;
}mtrng_state_t;

mtrng_state_t *mtrng_state;

void open_rng(unsigned long int s)
{
     mtrng_state=(mtrng_state_t*)malloc(sizeof(mtrng_state_t));
     if (s == 0)
	  s = 4357;   /* the default seed is 4357 */
     mtrng_state->mt[0]= s & 0xffffffffUL;
     int i;
     for (i=1; i<MTRNG_N; i++)
     {
	 /* See Knuth's "Art of Computer Programming" Vol. 2, 3rd
	    Ed. p.106 for multiplier. */
	  mtrng_state->mt[i] =
	       (1812433253UL * 
		(mtrng_state->mt[i-1] ^ (mtrng_state->mt[i-1] >> 30)) + i);
	  mtrng_state->mt[i] &= 0xffffffffUL;
     }
     mtrng_state->mti = i;
}

void close_rng()
{
     free(mtrng_state);
}

unsigned long int get_deviate()
{
     unsigned long k;
     unsigned long int *const mt = mtrng_state->mt;
#define MTRNG_MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)
     if (mtrng_state->mti >= MTRNG_N)
     {   // generate N words at one time 
	  int kk;
	  for (kk=0; kk<MTRNG_N-MTRNG_M; kk++)
	  {
	       unsigned long y = (mt[kk] & MTRNG_UPPER_MASK) | 
		    (mt[kk + 1] & MTRNG_LOWER_MASK);
	       mt[kk] = mt[kk + MTRNG_M] ^ (y >> 1) ^ MTRNG_MAGIC(y);
	  }
	  for (; kk<MTRNG_N-1; kk++)
	  {
	       unsigned long y = (mt[kk] & MTRNG_UPPER_MASK) | 
		    (mt[kk + 1] & MTRNG_LOWER_MASK);
	       mt[kk] = mt[kk + (MTRNG_M - MTRNG_N)] ^ (y >> 1) ^ 
		    MTRNG_MAGIC(y);
	  }
	  {
	       unsigned long y = (mt[MTRNG_N - 1] & MTRNG_UPPER_MASK) | 
		    (mt[0] & MTRNG_LOWER_MASK);
	       mt[MTRNG_N - 1] = mt[MTRNG_M - 1] ^ (y >> 1) ^ MTRNG_MAGIC(y);
	  }
	  mtrng_state->mti = 0;
     }
     // Tempering
     k = mt[mtrng_state->mti];
     k ^= (k >> 11);
     k ^= (k << 7) & 0x9d2c5680UL;
     k ^= (k << 15) & 0xefc60000UL;
     k ^= (k >> 18);
     mtrng_state->mti++;
     return k;
}

double rand_double()
{
     return get_deviate() / 4294967296.0 ;
}

double rand_double(double a, double b)
{
     double x=rand_double();
     return (b-a)*x+a;
}

// Get a uniform integer deviate between a and b,
// a and b included!
long rand_int(long a, long b)
{
     long range=b-a+1;
     return (get_deviate() % range)+a;
}

double sqr(double a)
{
     return a*a;
}

////////////////////////////////////////////////////////
// BINARY NUMBERS

long pow2(int num)
// fast method to return a desired power of 2
{
     return(1<<num);
}

bool get_bit(long X, int pos)
{
     return (X & pow2(pos));
}

int get_max_bit(long X)
{ 
     int i=-1;
     if (!X) return (0);
     do
     {
	  i++;
	  X=X>>1;
     }while(X);
     return i;
}

long flip_bit(long X, int pos) 
{
     long power=pow2(pos);
     if (X & power) return(X-power);
     else return(X+power);
}

long put_bit(long X, int pos, int value)
{
     if (get_bit(X,pos)==value) return X;
     else return flip_bit(X,pos);
}

long swap_bits(long X, int pos1, int pos2)
{
     bool temp=get_bit(X,pos1);
     X=put_bit(X,pos1,get_bit(X,pos2));
     X=put_bit(X,pos2,temp);
     return X;
}

long reverse_bits(long X, int maxbit)
{
     if (!maxbit) maxbit=get_max_bit(X);
     int x=0;
     for (int i=0;i<=maxbit;i++)
	  if (get_bit(X,i)) 
	       x=put_bit(x,maxbit-i,1);
     return x;
}

// insert bit "b" in position "j" of number "X"
// e.g.: insert(ABCD,2,0) -> AB0CD
long insert_bit(long X, int pos, int value)
{
     // Let X=ABCD (15), j=2
     long I=1<<pos; // I=100 (4)
     long I2=I-1; // I2=11 (3)
     long low=X & I2; // low=CD  (3)
     long x=X-low; // x=AB00      (12)
     x*=2;        // x=AB000      (24)
     x+=low; // AB0CD             (27)
     if (value) x+=I; // AB1CD        (31)
     return x;
}

int count_ones(long X) // returns the number of ones in the expansion of X
{
     int number=0;
     long C=X;
     do
     {
	  if (C%2) number++;
	  C>>=1;
     }while(C);
     return number;
}

long next_in_sector(long X)
// returns the smallest number which is bigger than X and has the same
// number of ones in its expansion. 
{
     long i=0, n_ones=0, Y=X, R=X;
     bool bit=Y%2, nextbit;
     do
     {
	  Y>>=1;
	  nextbit=Y%2;
	  if (bit && !nextbit) break; // if no 01 pattern is found yet
	  if (bit) // if current bit is switched
	  {
	       n_ones++;
	       R-=pow2(i);
	  }
	  i++;
	  bit=nextbit;
     }while(Y); 
     X=swap_bits(R,i+1,i);
     X+=pow2(n_ones)-1;
     return X;	  
}

char* B2string(long X, int maxbit)
{
     if (!maxbit) maxbit=get_max_bit(X);
     char *s=new char[maxbit+1];
     for (int i=maxbit;i>=0;i--)
	  s[maxbit-i]=(get_bit(X,i) ? '1' : '0');
     s[maxbit+1]='\0';
     return s;
}

char* B2stringR(long X, int maxbit)
{
     if (!maxbit) maxbit=get_max_bit(X);
     char *s=new char[maxbit+1];
     for (int i=0;i<=maxbit;i++)
	  s[i]=(get_bit(X,i) ? '1' : '0');
     s[maxbit]='\0';
     return s;
} 

void printB(long X)
{
     char *s=B2string(X);
     printf("%s",s);
     delete s;
}

////////////////////////////////////////////////////////
// Class List

List::List()
{
     N=0; // non-created
     D=NULL;
}

List::List(long n)
{
     N=0;
     Create(n);
     Zero();
}

List::List(long* P)
{
     N=0;
     Create(P[0]); // zero element is assumed to contain the number of items
     memcpy(D,P,(N+1)*sizeof(long));
}

List::List(long* P, long n)
{
     N=0;
     Create(n);
     memcpy(D,P,(N+1)*sizeof(long));
}

List::List(const List& L)
{
     N=0;
     Create(L.N);
     if (L.N>0) memcpy(D,L.D,(N+1)*sizeof(long));
     else D=NULL;
}

void List::Create()
{
    N=0;
    D=(long*)NULL;
}

void List::Create(long n)
{
     if (N!=0) Destroy();
     D=(long*)malloc((n+1)*sizeof(long));
     D[0]=n;
     N=n;
}

void List::Start() // n=0
{
     D=(long*)NULL;
     N=0;
     Create();
}

void List::Destroy()
{
     if (N) free(D);
     D=NULL;
     N=0;
}

void List::Zero()
{
     if (N==0) return;
     for (long i=1;i<=N;i++)
	  D[i]=0;
}

List::~List()
{
     Destroy();
}

long List::Count(long x) const // counts the number of appearances of x
{
     long c=0;
     for (long i=1;i<=N;i++)
	  if (D[i]==x) c++;
     return c;
}

long List::Find(long n, long k) const // finds first appearance of "n"
{
     long i=k-1;
     if (N==0) return 0;
     do
     {
	 i++; 
     }while(i<N && D[i]!=n);
     if (i<=N && D[i]==n) return i;
     else return 0;
}



long List::Remove_All(long x) // removes all appearances of "n"
{
     if (D==(long*)NULL) return 0;
     long i,j,n;
     long *D2;
     n=0;
     for (i=1;i<=N;i++)
	  if (D[i]==x) n++;  
     D2=(long*)malloc((N-n+1)*sizeof(long));
     j=1;
     for (i=1;i<=N;i++)
	 if (D[i]!=x) 
	 {
	      D2[j]=D[i];
	      j++;
	 }
     free(D);
     D=D2;
     D[0]=N-n;
     N-=n;
     return n;
}

long List::Min() const
{
     long nmin=D[1];
     for (long k=2;k<=N;k++)
	  if (D[k]<nmin) nmin=D[k];
     return nmin;
}

long List::Max() const
{
     long nmax=D[1];
     for (long k=2;k<=N;k++)
	  if (D[k]>nmax) nmax=D[k];
     return nmax;
}

long List::Max_Index() const
{
     long nmax=D[1];
     long kmax=1;
     for (long k=2;k<=N;k++)
	  if (D[k]>nmax) { nmax=D[k]; kmax=k; }
     return kmax;
}

long List::Min_Index() const
{
     long nmin=D[1];
     long kmin=1;
     for (long k=2;k<=N;k++)
	  if (D[k]<nmin) { nmin=D[k]; kmin=k; }
     return kmin;
}

int order12 (const void *a, const void *b)
{
    const long *da = (const long *) a;
    const long *db = (const long *) b;
    return (*da > *db) - (*da < *db);
}

int order21(const void* x, const void* y) // descending order
{
    const long *xx=(const long *) x;
    const long *yy=(const long *) y;
    return (*yy > *xx) - (*xx > *yy);
}

void List::Order(long tipo) 
// if positive: ascending order; if negative, descending
{
     if (tipo>0)
	  qsort(D+1,N,sizeof(long),order12);
     else qsort(D+1,N,sizeof(long),order21);
}
     
void List::Append(long n)
{ 
     D=(long*)realloc(D,(N+2)*sizeof(long));
     D[N+1]=n;
     N++;
}

void List::Append(const List &L)
{
     D=(long*)realloc(D,(N+L.N+1)*sizeof(long));
     memcpy(D+N+1,L.D+1,L.N*sizeof(long));
     N+=L.N;
}

void List::Insert(long p, long x)
{
     D=(long*)realloc(D,(N+2)*sizeof(long));
     long *F=(long*)malloc((N-p+1)*sizeof(long));
     memcpy(F,D+p,(N-p+1)*sizeof(long));
     D[p]=x;
     memcpy(D+p+1,F,(N-p+1)*sizeof(long));
     N++;
     free(F);
}

void List::Insert(long p, const List &K)
{
     // Save the data beyond "p"
     long *F=(long*)malloc((N-p+1)*sizeof(long));
     memcpy(F,D+p,(N-p+1)*sizeof(long));
     // Now, extend D and copy K on D
     long n=K.N;
     D=(long*)realloc(D,(N+n+1)*sizeof(long));
     memcpy(D+p,K.D+1,n*sizeof(long));
     // Now, copy again the data beyond "p"
     memcpy(D+p+n,F,(N-p+1)*sizeof(long));
     N+=n;
     free(F);
}

void List::Write() const // l means length of each line of print
{
     Write(stdout);
}

void List::Write(FILE *arch) const
{
     if (!N) printf("Empty list\n");
     int l=10, i;
     for (i=1;i<=N;i++)
     {
	  fprintf(arch,"%5ld ",D[i]);
	  if (i%l==0) fprintf(arch,"\n");
     }
     if ((i-1)%l!=0) fprintf(arch,"\n");

}

long List::Find_Non(long x, long k) const
{
     long i=k-1;
     do
	  i++;
     while(i<N && D[i]==x);
     if (i<=N && D[i]!=x) return i;
     else return 0;
//     if (i>N) return 0;
//     else return i;     
}

long List::Find_Higher(long x, long k) const
{
     long i=k-1;
     do
	  i++;
     while(i<N && D[i]<=x);
     if (i<=N && D[i]>x) return i;
     else return 0;
//     if (i>N) return 0;
//     else return i;     
}

long List::Find_Lower(long x, long k) const
{
     long i=k-1;
     do
	  i++;
     while(i<N && D[i]>=x);
//     if (i>N) return 0;
//     else return i;     
     if (i<=N && D[i]<x) return i;
     else return 0;
}

void List::Uniquify()
{
     long i,j,x,counter=0;
     long *D2, *D3;
     D2=(long*)malloc((N+1)*sizeof(long));
     for (i=1;i<=N;i++)
     {
	  D2[i]=1;
	  for (j=1;j<i;j++)
	       if (D[i]==D[j]) D2[i]=0;
	  if (D2[i]==1) counter++;
     } // now D2[i] contains "1" at element "i" only if "i" is original
     D3=(long*)malloc((counter+1)*sizeof(long));
     x=1;
     for (i=1;i<=N;i++)
	  if (D2[i]==1) 
	  {
	       D3[x]=D[i];
	       x++;
	  }
     D3[0]=counter;
     free(D);
     free(D2);
     D=D3;
     N=counter;
}

void List::Swap(long i, long j)
{
     if (i>N || j>N) merror("Swap error!\n");
     long tmp=D[i];
     D[i]=D[j];
     D[j]=tmp;
}

void List::Part(long n0, long n1)
{
     N=n1-n0+1;
     if (N<1)
     {
	  N=0;
	  free(D); D=NULL;
	  return;
     }
     long* D2=(long*)malloc((N+1)*sizeof(long));
     D2[0]=N;
     memcpy(D2+1,D+n0,N*sizeof(long));
     free(D);
     D=D2;
}

List& List::operator=(const List& L0) 
{
     if (this==&L0) return *this;
     Copy(*this,L0);
     return *this;
}

long List::operator() (long n1) const
{
    if (n1<0 || n1>N) 
	 merror("Error getting list component.\n");
    return D[n1];
}

long& List::operator() (long n1)
{
    if (n1<0 || n1>N) 
	 merror("Error putting list component.\n");
    return D[n1];
}

void List::operator+=(const List &L)
{
     if (!N) { Create(L.N); Zero(); }
#ifdef DEBUG
      if (N!=L.N) merror("Incompatible sizes in List +=");
#endif
      for (long i=1;i<=N;i++)
	   D[i]+=L(i);
}

void List::operator-=(const List &L)
{
     if (!N) { Create(L.N); Zero(); }
#ifdef DEBUG
      if (N!=L.N) merror("Incompatible sizes in List +=");
#endif
      for (long i=1;i<=N;i++)
	   D[i]-=L(i);
}

void List::operator*=(const long p)
{
     if (!N) return;
     for (long i=1;i<=N;i++)
	  D[i]*=p;
}

void List::operator/=(const long p)
{
     if (!N) return;
     for (long i=1;i<=N;i++)
	  D[i]/=p;
}

void List::operator+=(const long p)
{
     if (!N) return;
     for (long i=1;i<=N;i++)
	  D[i]+=p;
}

void List::operator-=(const long p)
{
     if (!N) return;
     for (long i=1;i<=N;i++)
	  D[i]-=p;
}


void Copy(List& L1, const List &L2)
{
     L1.Destroy();
     L1.N=L2.N;
     if (L2.D!=NULL)
     {
	  L1.D=(long*)malloc((L1.N+1)*sizeof(long));
	  memcpy(L1.D,L2.D,(L1.N+1)*sizeof(long));
     }
}

// Returns true if lists are completely equal
bool Is_Equal(const List &L1, const List &L2)
{
     if (L1.N!=L2.N) return false;
     for (long i=1;i<=L1.N;i++)
	  if (L1(i)!=L2(i)) return false;
     return true;
}

// returns true if L1 is contained in L2
// more precisely, if all elements of L1 are also in L2
bool Is_Contained(const List &L1, const List &L2)
{
     for (int i=1;i<=L1.N;i++)
	  if (!L2.Find(L1(i))) return false;
     return true;
}

List Intersection(const List &L1, const List &L2)
{
     List L;
     if (L1.N*L2.N==0) return L;
     for (long i=1;i<=L1.N;i++)
	  if (L2.Find(L1(i))>0) L.Append(L1(i));
     return L;
}

void List::Remove_Index(long k)
{
     if (N==1) { Destroy(); return; }
     memmove(D+k,D+k+1,(N-k)*sizeof(long));
     D=(long*)realloc(D,N*sizeof(long));
     N--;
}

long List::Remove_Higher(long x) // remove all values higher than int
{
     long i,j,n;
     long *D2;
     n=0;
     for (i=1;i<=N;i++)
	  if (D[i]>x) n++;  
     D2=(long*)malloc((N-n+1)*sizeof(long));
     j=1;
     for (i=1;i<=N;i++)
     {
	 if (D[i]<=x) 
	 {
	      D2[j]=D[i];
	      j++;
	 }
     }
     free(D);
     D=D2;
     D[0]=N-n;
     N-=n;
     return n;
}


long List::Remove_Lower(long x)
{
     long i,j,n;
     long *D2;
     n=0;
     for (i=1;i<=N;i++)
	  if (D[i]<x) n++;  
     D2=(long*)malloc((N-n+1)*sizeof(long));
     j=1;
     for (i=1;i<=N;i++)
     {
	 if (D[i]>=x) 
	 {
	      D2[j]=D[i];
	      j++;
	 }
     }
     free(D);
     D=D2;
     D[0]=N-n;
     N-=n;
     return(n);
}

int List::Save_Binary(FILE *arch) const
{
     int ausgang=fwrite(&N,sizeof(long),1,arch);
     if (ausgang!=1) return 0;
     if (!N) return 1;
     ausgang=fwrite(D,sizeof(long),N+1,arch);
     if (ausgang!=N+1) return 0;
     return 1;
}

int List::Save_Binary(const char *name) const
{
     FILE* fich=fopen(name,"wb");
     int status=Save_Binary(fich);
     fclose(fich);
     return status;
}

int List::Load_Binary(FILE *arch)
{
     int ausgang=fread(&N,sizeof(long),1,arch);
     if (ausgang!=1) return 0;
     if (!N) return 1;
     Create(N);
     ausgang=fread(D,sizeof(long),N+1,arch);
     if (ausgang!=N+1) return 0;
     return 1;
}

int List::Load_Binary(const char *name)
{
     FILE* fich=fopen(name,"rb");
     int status=Load_Binary(fich);
     fclose(fich);
     return status;
}

int List::Load(const char *s)
{
     FILE *fich=fopen(s,"rt");
     int status=Load(fich);
     fclose(fich);
     return status;
}

int List::Load(FILE *fich)
{
     if (!fich) return 0;
     Create();
     long x, n;
     if (!fscanf(fich,"%ld\n",&n))
	  merror("Error in List::Load");
     while(!feof(fich))
     {
	  if (!fscanf(fich,"%ld\n",&x))
	       merror("Error in List::Load");
	  Append(x);
     }
     if (N!=n) return 0; // error
     return 1;
}

int List::Save(const char *name) const
{
     FILE *fich=fopen(name,"wt");
     int status=Save(fich);
     fclose(fich);
     return status;
}

int List::Save(FILE *fich) const
{
     if (!fich) return 0;
     fprintf(fich,"%ld\n",N);
     for (long i=1;i<=N;i++)
	  fprintf(fich,"%ld\n",D[i]);
     return 1;
}

void List::Reverse()
{
     for (long i=1;i<=N/2;i++)
	  SWAP(D[i],D[N+1-i]); 
}

long List::Prod(long i0, long i1) const
{
     long r=1;
     for (long k=i0;k<=i1;k++)
	  r*=D[k];
     return r;
}

long List::Prod() const
{
     return Prod(1,N);
}

long List::Sum(long i0, long i1) const
{
     long r=0;
     for (long k=i0;k<=i1;k++)
	  r+=D[k];
     return r;
}

long List::Sum() const
{
     return Sum(1,N);
}

void List::Set_Value(long S) 
{
     for (long k=1;k<=N;k++)
	  D[k]=S;
}


List operator+(const List &L1, const List &L2) 
{
     List R(L1);
     R+=L2;
     return R;
}

List operator-(const List &L1, const List &L2)
{
     List R(L1);
     R-=L2;
     return R;
}

List operator+(const List &L, const long p) 
{
     List R(L);
     R+=p;
     return R;
}

List operator-(const List &L, const long p)
{
     List R(L);
     R-=p;
     return R;
}

List operator+(const long p, const List &L) 
{
     List R(L);
     R+=p;
     return R;
}

List operator-(const long p, const List &L)
{
     List R(L);
     R-=p;
     return R;
}

List operator-(const List &L)
{
     List R(L);
     R*=-1;
     return L;
}

List operator*(const List &L, const long p)
{
     List R(L);
     R*=p;
     return R;
}

List operator/(const List &L, const long p)
{
     List R(L);
     R/=p;
     return R;
}

List operator*(const long p, const List &L)
{
     List R(L);
     R*=p;
     return R;
}

List List_Interval(long i1, long i2)
{
     List L(labs(i2-i1)+1);
     if (i2>i1)
	  for (long i=i1;i<=i2;i++) L(i-i1+1)=i;
     else
	  for (long i=i1;i>=i2;i--) L(i1-i+1)=i;
     return L;
}

// The resulting list is as large as L1
// For each element of L1, it looks up its value in L2
List Combination(const List &L1, const List &L2)
{
     List L(L1);
     for (long i=1;i<=L1.N;i++)
	  L(i)=L2(L1(i));
     return L;
}

// remove from list L1 all elements in L2
List Substraction(const List &L1, const List &L2)
{
     List L(L1);
     for (long i=1;i<=L2.N;i++)
     {
	  long x=L2(i);
	  if (L.Find(x)) L.Remove_All(x);
     }
     return L;
}

List Reverse(const List &L)
{
     List L2(L);
     L2.Reverse();
     return L2;
}

List Part(const List &L, long i1, long i2)
{
     List L2(L);
     L2.Part(i1,i2);
     return L2;
}


// get the list of all elements from 1 to L not in Ll
// List Complement(const List &Ll, long L)
// {
//      List Lr;
//      for (long i=1;i<=L;i++)
// 	  if (!Ll.Find(i)) Lr.Append(i);
//      return Lr;
// }


List  B2List(long X, int maxbit) // conversion of a binary number to a list
{
     List L;
     long Y=X;
     while(Y)
     {
	  if (Y%2) L.Append(1);
	  else L.Append(0);
	  Y/=2;
     }
     // now if size of list is nbits (or nbits is 0) done, otherwise, complete
     while (L.N<maxbit+1)
	    L.Append(0);
     return L;
}

long List2B(const List &L)
{
     long X=0;
     for (long i=1;i<=L.N;i++)
	  X+=L(i)*pow2(i-1);
     return X;
}

// Numbering system routines.
// Takes number given by list X in numbering system Ls and returns a long
long To_Number(const List &X, const List &Ls)
{
     long N=Ls.N;
     long ix=X(N);
     if (Ls.N<=1) return ix;
     long prod=1;
     for (long k=Ls.N-1;k>=1;k--)
     {
	  prod*=Ls(k+1);
	  ix+=X(k)*prod;
     }
     return ix;
}

// Opposite, takes a long and returns a list representing it in number system Ls
List To_List(long ix, const List &Ls)
{
     long i=ix;
     List X(Ls.N);
     for (long k=Ls.N;k>=1;k--)
     {
	  X(k)=i % Ls(k);
	  i/=Ls(k);
     }
     return X;
}
