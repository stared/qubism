// Parsing routines
// 111130
#ifndef TEXT_H
#define TEXT_H
#include"common.h"

class String
{
public:
     long N;
     char *D;
     String();
     String(const char *);
     String(const String &);
     ~String();
     void Start();
     void Create(const char *);
     void Destroy();

     long Get_Line(FILE *); // read a line from a file, ret -1 is failure
     void Write();

     void Append(const char *);
     void Append(const char *, long);
     void Append(const String &);
     void Append(const char);

     void Append_F(const char *, long i);
     void Append_F(const char *, double f);

     String& operator=(const String&);
     String& operator=(const char *);

     void Strip_Blanks();

     bool Is_Here(const char *);
     bool Is_Here(const String &);
     bool Is_There(const char *);
     bool Is_There(const String &);

     String Token(char q0, char q1, long n) const;
     String Part(long i0, long i1) const;

     long Get_Long(char q0, char q1, long n) const;
     double Get_Double(char q0, char q1, long n) const;

     long Count(char q) const;
     long Find_Nth(char q, long n) const;
     long To_Long() const;
     double To_Double() const;
     
     void To_LowerCase(); 
     void To_UpperCase();
};

void Copy(String &S, const String &S1);
void Copy(String &S, const char *z);
void Copy(String &S, const char *z, long n);

ostream& operator<< (ostream &out, String &S);
bool operator==(const String &S1, const String &S2);


#endif
