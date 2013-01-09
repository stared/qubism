// Parsing routines
// 111130

#include"text.h"

String::String()
{
     N=0; D=(char*)NULL;
}

String::String(const char *z)
{
     Start();
     Create(z);
}

String::String(const String &S)
{
     Start();
     Create(S.D);
}

void String::Start()
{
     N=0; D=(char*)NULL;
}

String::~String()
{
     Destroy();
}

void String::Destroy()
{
     if (N) free(D);
     N=0;
}

void String::Create(const char *z)
{
     Destroy();
     if (!z) return;
     N=strlen(z);
     D=(char*)malloc(N+1);
     strncpy(D,z,N);
     D[N]=0;
}

// long is errorcode: chars if read smth (can be 0); -1 if smth went wrong
long String::Get_Line(FILE *fich)
{
     Destroy();
     char *z=NULL;
     ssize_t Nr; size_t N=0;
     Nr=getline(&z,&N,fich);
     if (Nr==-1) return -1;
     if (Nr) Copy(*this,z,(long)Nr-1);
     return Nr;
}
 
void String::Write()
{
     printf("%s\n",D);
}

// Append a string, until you get \0
void String::Append(const char *s)
{
     if (!s) return;
     long n=strlen(s);
     if (!D) D=(char*)malloc(n+1);
     else D=(char*)realloc(D,N+n+1);
     strncpy(D+N,s,n);
     N+=n;
     D[N]=0;
}

// Append n chars from a string
void String::Append(const char *s, long n)
{
      if (!s) return;
      D=(char*)realloc(D,N+n+1);
      strncpy(D+N,s,n);
      N+=n;
}

// Append another string
void String::Append(const String &S)
{
     Append(S.D);
}

// Append the string resulting from printing "i" with a given format
void String::Append_F(const char *formato, long i)
{
     char *s=(char*)malloc(100);
     sprintf(s,formato,i);
     Append(s);
     free(s);
}

// Append the string resulting from printing "i" with a given format
void String::Append_F(const char *formato, double x)
{
     char *s=(char*)malloc(100);
     sprintf(s,formato,x);
     Append(s);
     free(s);
}

// Return S from chars i0 to i1, both included (first char is 0!!!)
String String::Part(long i0, long i1) const
{
     String S2;
     if (i0>i1) return S2; 
     long n=i1-i0+1;
     S2.Append(D+i0,n);
     return S2;
}

String& String::operator=(const String &S)
{
     //if (this==&S) return *this;
     Copy(*this,S);
     return *this;
}

String& String::operator=(const char *z)
{
     Create(z);
     return *this;
}

ostream& operator<< (ostream &out, String &S)
{
     out << S.D;
     return out;
}


void Copy(String &S, const String &S1)
{
     S.Destroy();
     S.N=S1.N;
     if (S1.D!=NULL)
     {
	  S.D=(char*)malloc(S1.N+1);
	  strncpy(S.D,S1.D,S1.N);
	  S.D[S.N]=0;
     }
}

void Copy(String &S, const char *z)
{
     S.Destroy();
     if (!z) return;
     long n=strlen(z);
     S.D=(char*)malloc(n+1);
     strcpy(S.D,z);
     S.N=n;
}

void Copy(String &S, const char *z, long n)
{
     S.Destroy();
     if (!z) return;
     S.D=(char*)malloc(n+1);
     strncpy(S.D,z,n);
     S.D[n]=0;
     S.N=n;
}

// remove blanks from beginning and end of char, free the pointer
void String::Strip_Blanks()
{
     long i0=0;
     while(D[i0]==' ')
	  i0++;  // now D[i0]!=' '
     long i1=strlen(D)-1;
     while(D[i1]==' ')
	  i1--;
     long n=i1-i0+1;
     char *z=(char*)malloc(n+1);
     strncpy(z,D+i0,n);
     z[n]=0;
     Destroy();
     D=z;
     N=n;
}

void String::Append(const char q)
{
     D=(char*)realloc(D,N+2);
     D[N]=q;
     D[N+1]=0;
     N++;
}

bool operator==(const String &S1, const String &S2)
{
     if (S1.N!=S2.N) return false;
     return !strcmp(S1.D,S2.D);
}

bool operator!=(const String &S1, const String &S2)
{
     return !(S1==S2);
}

bool String::Is_Here(const char *z)
{
     if (!N) return false;
     return(memcmp(D,z,strlen(z))==0);
}

bool String::Is_Here(const String &S)
{
     return Is_Here(S.D);
}

bool String::Is_There(const char *z)
{
     if (!N) return false;
     if (strstr(D,z)!=NULL) return true;
     else return false;
}

bool String::Is_There(const String &S)
{
     return Is_There(S.D);
}

// return the substring from the n-th appearance of q0 
// until the next appearance of q1
// if q0==0, it means "from the beginning"; if q1==0, "to the end"
// if no appropriate q0 is found, returns empty String
// if no appropriate q1 is found, returns from q0 to end of the String
String String::Token(char q0, char q1, long n) const
{
     long i0, i1; // limiting indices
     if (!q0) i0=-1;
     else
     {
	  i0=Find_Nth(q0,n);
	  if (i0==-1) { String S; return S; } // return empty string  
     }
     String S=Part(i0+1,N);
     if (!q1) i1=S.N;
     else
	  i1=S.Find_Nth(q1,1);
     if (i1==-1) i1=S.N;
     return S.Part(0,i1-1);
}

long String::Get_Long(char q0, char q1, long n) const
{

     String S=Token(q0,q1,n);
     return S.To_Long();
}

double String::Get_Double(char q0, char q1, long n) const
{
     String S=Token(q0,q1,n);
     return S.To_Double();
}

// count how many appearances of the given character
long String::Count(char q) const
{
     long n=0;
     for (long i=0;i<=N;i++) 
	  if (D[i]==q) n++;
     return n;
}

// find the n-th appearance of the character q in the string
// -1 means not found
long String::Find_Nth(char q, long n) const
{
     long vez=0, i;
     for (i=0;i<N;i++) 
     {
	  if (D[i]==q) vez++;
	  if (vez==n) break;
     }
     if (i==N) return -1;
     return i;
}

long String::To_Long() const
{
     return (long)atoi(D);
}

double String::To_Double() const
{
     return strtod(D,NULL); 
}

void String::To_LowerCase()
{
     for (long i=0;i<N;i++)
	  D[i]=tolower(D[i]);
}

void String::To_UpperCase()
{
     for (long i=0;i<N;i++)
	  D[i]=toupper(D[i]);
}
