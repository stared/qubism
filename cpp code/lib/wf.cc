// Routines to deal with wavefunctions in a standard way
// 110914

#include"wf.h"

void Write_WF(const CVector &Psi, long L, long s, const char *filename)
{
     FILE *fich=fopen(filename,"wt");
     if (!fich) merror("Couldn't create file for writing\n");
     fprintf(fich,"#L %ld\n",L);
     fprintf(fich,"#s %ld\n",s);
     for (long i=1;i<=Psi.N;i++)
	  fprintf(fich,"%ld %16.12g %16.12g\n",i-1,real(Psi(i)),imag(Psi(i)));
     fclose(fich);
}

// Read a wf: the CVector, L and s (number of states per site)
void Read_WF(CVector &V, long &L, long &s, const char *name)
{
     // First, find out the size that we'll need
     FILE *fich=fopen(name,"rt");
     if (!fich) merror("File does not exist\n");

     char *z=(char*)malloc(100);
     s=2; L=0; // default values
     bool createdvector=false;
     int err=0;
     do
     {
          err=get_next_line(&z,fich);
          if (err==-1) break;
          if (z[0]=='#')
	  {
	       if (z[1]=='L') 
	       {
		    L=get_int(z," ^",1);
	       }
	       if (z[1]=='s')
	       {
		    s=get_int(z," ^",1);
	       }
	       continue;
	  }
	  if (!createdvector) // the first line without a '#' contains data
	  {
	       if (s==2) V.Create(1<<L);
	       else V.Create((long)pow((double)s,(double)L));
	       V.Zero();
	       createdvector=true;
	  }
	  long i; double x, y;
	  sscanf(z,"%ld %lg %lg\n",&i,&x,&y);
	  if (i>=V.N) merror("Error reading wf file\n");
	  V(i+1)=x+M_I*y;
     }while(err!=-1);
     fclose(fich);
}
