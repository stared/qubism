// Routines for parsing 
// Last modification: 061017
#ifndef PARSING
#define PARSING

#include"parsing.h"

char* get_until(const char* s, const char c[])
// take a chain and a set of delimiters "c". Returns a new chain which
// contains the first characters of "s" until the first occurrence of "c"
{
     int i=-1;
     do
     {
	  i++;
     }while(memcmp(s+i,c,strlen(c))!=0 && s[i]!='\0');
     char *copia=(char*)malloc(i+1);
     memcpy(copia,s,i);
     copia[i]='\0';
     return copia;
}

char* get_token(const char* s, const char c[], int num=1)
// get_token takes a chain, a pair of delimiters and a number.
// Finds the "num"-th occurrence of c[0] and returns a new string
// containing the characters between that c[0] and the first appearance
// of c[1] afterwards.
// '^' is a symbol for a wildcard delimiter (terminator)
{
     int i, vez=0, l=(s ? strlen(s) : 0);
     char *cad;
     if (c[0]=='^') // from the beginning until the num-th time of c[1]
     {
	  if (c[1]=='^') // return a copy of the whole chain!
	  {
	       cad=(char*)malloc(l+1);
	       strcpy(cad,s);
	       return cad;
	  }
	  for (i=0;i<l;i++) 
	  {
	       if (s[i]==c[1]) vez++;
	       if (vez==num) break;
	  }
	  if (i==l) return NULL;
	  cad=(char*)malloc(i+1);
	  memcpy(cad,s,i);
	  cad[i]='\0';
	  return cad;
     }
     for (i=0;i<l;i++)
     {
	  if (s[i]==c[0]) vez++;
	  if (vez==num) break;
     }
     if (i==l) return NULL;
     int i0=i;
     if (c[1]=='^') // from c[0] to the end
     {
	  cad=(char*)malloc(l-i0+1);
	  memcpy(cad,s+i0+1,l-i0);
	  return cad;
     }
     do // typical case
     {
	  i++;
     }while(i<l && s[i]!=c[1]);
     if (i==l) return NULL; // c[1] was not found!
     cad=(char*)malloc(i-i0+1);
     memcpy(cad,s+i0+1,i-i0-1);
     cad[i-i0-1]='\0';
     return cad;
}

bool is_here(const char *s, const char *c)
// "s" is the big string and "c" is the smaller
// returns TRUE if "c" is the beginning of "s".
{
     if (!s) return false;
     return(memcmp(s,c,strlen(c))==0);
}

bool is_there(const char *s, const char *c)
// returns TRUE if "c" is somewhere inside "s"
{
     if (!s) return false;
     if (strstr(s,c)!=NULL) return true;
     else return false;
}

int get_int(const char *s, const char c[], int num=1)
// returns an integer within "s" delimited by c[0]..c[1]
{
     char *res;
     res=get_token(s,c,num);
     int result=atoi(res);
     free(res);
     return result;
}

double get_double(const char *s, const char c[], int num=1)
{
     char *res;
     res=get_token(s,c,num);
     double result=atof(res);
     free(res);
     return result;
}
     
int find_first(const char *s, char c)
// Returns the index to the first appearance of "c" into "s"
// if not found, returns -1.
{
     if (strlen(s)==0) return -1;
     int i=-1, length=strlen(s);
     do
     {
	  i++;
     }while(s[i]!=c && i<length);
     if (i<length) return i; 
     else return -1; 
}

char* part(const char *s, int i0, int i1)
// returns a new string object [i0,i1]
{
     if (i0>i1) return (char*)NULL;
     char *d=(char*)malloc(i1-i0+2);
     for (int i=0;i<=i1-i0;i++)
	  d[i]=s[i0+i];
     d[i1-i0+1]='\0';
     return d;
}

char* get_argument(const char* z, const char c[], int n)
// gets the n-th argument of chain "z", where each argument is 
// something included into the delimiters of "c". For example, if
// c="[]", then get_argument("[dasd[asd]asd][xxx]","[]",2) nos
// devuelve "xxx".
{
     int i=0;
     int parenthesis=0;
     int argument=0;
     int l=strlen(z);
     int j0=0, j1=0;
     do
     {
	  if (z[i]==c[0])
	  {
	       if (parenthesis==0) { j0=i+1; }
	       parenthesis++;
	  }
	  if (z[i]==c[1]) 
	  {
	       parenthesis--;
	       if (parenthesis==0)
	       {
		    argument++;
		    j1=i-1;
	       }
	  }
	  i++;
     }while(i<=l && argument<n);
     if (argument!=n) return (char*)NULL;
     return(part(z,j0,j1));
}
	  

int get_order(char*& command, char*& arg, char*& program)
// read a command and its set of arguments
// program "advances" an order!
{
     int i,command_end, total_end=strlen(program);
        
     command_end=find_first(program,'(');
     command=part(program,0,command_end-1);
     
     i=command_end;
     int parentesis=1;
     do
     {
	  i++;
	  if (program[i]=='(') parentesis++;
	  if (program[i]==')') parentesis--;
     }while(parentesis>=1 && i<total_end);
     arg=part(program,command_end+1,i-1);
     
     if (i>=total_end-2) {program=NULL; return 1;}
     program+=i+1;
     return 0;     
}
     
	  
char* compose(char*c1, char*c2)
{
     int l1=strlen(c1), l2=strlen(c2);
     char *nuevo=(char*)malloc(l1+l2+1);
     memcpy(nuevo,c1,l1);
     memcpy(nuevo+l1,c2,l2);
     nuevo[l1+l2]='\0';
     return nuevo;
}

void left_add(char*& principal, char* piece)
{
     char *nuevo=compose(piece,principal);
     free(principal);
     principal=nuevo;
}

void right_add(char*& principal, char* piece)
{
     char *nuevo=compose(principal,piece);
     free(principal);
     principal=nuevo;
}

//////////////////////////////////////////////////////////////////////
// Reading the file
//////////////////////////////////////////////////////////////////////

// Adds a character to a string at position "k",
// *N contains the mallocked size. Reallocates if needed.
void add_char(char **puntero, size_t *N, int k, char c)
{
     if (k>=(int)*N)
     {
          puntero[0]=(char*)realloc(puntero[0],k+1);
          (*N)=k;
     }
     puntero[0][k]=c;
}


// get_line, reads a full line from a file and returns the full size
// of the line, possibly zero.
// N must contain the mallocked size of *puntero, 
// Return value is the new size of puntero (N if it was enough)
// or -1 if smth failed or the end of the file was reached bfr reading
ssize_t get_line(char **puntero, size_t *N, FILE *archivo)
{
     if (feof(archivo)) return -1;
     char c;
     if (!fscanf(archivo,"%c",&c))
	  merror("Error in get_line");
     if (c=='\n') return 0;
     int i=0;
     do
     {
          add_char(puntero,N,i,c);
          if (!fscanf(archivo,"%c",&c))
	       merror("Error in get_line");
          i++;
     }while(!feof(archivo) && c!='\n');
     add_char(puntero,N,i,0);
     return i;
}


int get_next_line(char **z, FILE *fich)
{     
     size_t st_size=(size_t)1024; 
     z[0]=(char*)realloc(z[0],st_size);
     int error;
     error=(int)get_line(z,&st_size,fich);
     return (error);
}

#endif




