////////////////////////////////////////////////////////
// hvb++ 1.0
// Copyleft: Javier Rodríguez Laguna
// 080725

// Graph structure, without geometry! Only topology!
// Last modification: 100418

#ifndef GRAPH
#define GRAPH

#include"graph.h"

/////////////////////////////////////////////////////////////////////////
///////      IMPLEMENTATIONS
/////////////////////////////////////////////////////////////////////////

Graph::Graph()
{
     N=0; Nl=0;
     V=(long**)NULL; I=(long*)NULL;
     updated=false;
}

Graph::Graph(long n)
{
     N=0; Nl=0; V=(long**)NULL; I=(long*)NULL; updated=false;
     for (long i=1;i<=n;i++)
	  New_Site();
}

Graph::Graph(const Graph &G)
{
     N=0;
     Copy(*this,G);
}

void Copy(Graph &G, const Graph &G2)
{
     if (G.N!=0) G.Destroy();
     G.N=G2.N;
     G.Nl=G2.Nl;
     G.updated=G2.updated;
     G.V=(long**)malloc((G2.N+1)*sizeof(long*));
     for (long i=1;i<=G.N;i++)
     {
	  long n=G2.V[i][0];
	  G.V[i]=(long*)malloc( (n+1)*sizeof(long));
	  memcpy(G.V[i],G2.V[i],(n+1)*sizeof(long));
     }
    if (G2.I!=NULL)
     {
	  G.I=(long*)malloc((G2.N+1)*sizeof(long));
	  memcpy(G.I,G2.I,(G2.N+1)*sizeof(long));
     }
}

Graph::~Graph()
{
     Destroy();
}

void Graph::Destroy()
{
     if (N==0) return;
     for (long i=1;i<=N;i++)
	  free(V[i]);
     free(V);
     if (updated) free(I);
     V=NULL;
     I=NULL;
     updated=false;
     N=0; Nl=0;
}

void Graph::Clean()
{
     if (!N) return;
     if (updated) { free(I); I=(long*)NULL; }
     for (long i=1;i<=N;i++)
	  free(V[i]);
     for (long i=1;i<=N;i++)
     {
	  V[i]=(long*)malloc(sizeof(long));
	  V[i][0]=0;
     }
     updated=false;
     Nl=0;
}
	  
bool Graph::New_Link(long s1, long s2)
{
     if (s1<1 || s1>N || s2<1 || s2>N) return 0;
     if (Is_Link(s1,s2)) return 0;
     long n1=V[s1][0], n2=V[s2][0];
     n1++; n2++;
     V[s1]=(long*)realloc(V[s1],(n1+1)*sizeof(long));
     V[s2]=(long*)realloc(V[s2],(n2+1)*sizeof(long));
     V[s1][0]++; V[s2][0]++;
     V[s1][n1]=s2;
     V[s2][n2]=s1;
     Nl++; 
     updated=false;
     return 1; 
}

void Graph::Remove_Link(long s1, long s2)
{
     if (!Is_Link(s1,s2)) return;
     List L1(V[s1]);
     List L2(V[s2]);
     L1.Remove_All(s2);
     L2.Remove_All(s1);
     free(V[s1]);
     free(V[s2]);
     V[s1]=(long*)malloc((L1.N+1)*sizeof(long));
     V[s2]=(long*)malloc((L2.N+1)*sizeof(long));
     memcpy(V[s1],L1.D,(L1.N+1)*sizeof(long));
     memcpy(V[s2],L2.D,(L2.N+1)*sizeof(long));
     updated=false;
     Nl--;
}

// stem is 0 by default
long Graph::New_Site(long stem)
{
     N++;
     if (N==1)
	  V=(long**)malloc(2*sizeof(long*));
     else
	  V=(long**)realloc(V,(N+1)*sizeof(long*));
   
     V[N]=(long*)malloc(sizeof(long));
     V[N][0]=0;
     if (stem!=0 && stem<=N-1)
	  New_Link(N,stem);
     return N;
}

void Graph::Add_Sites(long n)
{
     V=(long**)realloc(V,(N+n+1)*sizeof(long*));
     for (long i=1;i<=n;i++)
     {
	  V[N+i]=(long*)malloc(sizeof(long));
	  V[N+i][0]=0;
     }
     N+=n;
}
     
bool Graph::Is_Link(long s1, long s2) const
// returns 1 if there is a link from s1 to s2.
{
     long i=0, n1=V[s1][0];
     if (n1==0) return 0; 
     do
	  i++;
     while(V[s1][i]!=s2 && i<n1);
     return (V[s1][i]==s2);
}
	  
void Graph::Write() const 
{
     printf("Graph with %ld sites and %ld links\n",N,Nl);
     for (long i=1;i<=N;i++)
     {
	 printf("%ld (%ld): ",i,V[i][0]);
	 for (long j=1;j<=V[i][0];j++)
	     printf("%ld ",V[i][j]);
	 printf("\n");
     }
     printf("\n");
     if (!updated) return;
     printf("Links: \n");
     long k1, k2;
     for (long i=1;i<=Nl;i++)
     {
	  Get_Link_Sites(k1,k2,i);
	  printf("Link %ld: (%ld-%ld)\n",i,k1,k2);
     }
     printf("\n");
}
     
void Graph::Save(const char *name) const
{
     FILE *fich=fopen(name,"wt");
     fprintf(fich,"Graph with %ld sites and %ld links\n",N,Nl);
     for (long i=1;i<=N;i++)
     {
	 fprintf(fich,"%ld (%ld): ",i,V[i][0]);
	 for (long j=1;j<=V[i][0];j++)
	     fprintf(fich,"%ld ",V[i][j]);
	 fprintf(fich,"\n");
     }
     fprintf(fich,"\n");
     if (!updated) {  fclose(fich); return; }
     fprintf(fich,"Links: \n");
     long k1, k2;
     for (long i=1;i<=Nl;i++)
     {
	  Get_Link_Sites(k1,k2,i);
	  fprintf(fich,"Link %ld: (%ld-%ld)\n",i,k1,k2);
     }
     fprintf(fich,"\n");
     fclose(fich);
}


Graph& Graph::operator=(const Graph &G)
{
     if (&G==this) return *this;
     Copy(*this,G);
     return *this;
}
     
long Graph::operator()(long n) const
{
     return V[n][0];
}

long Graph::operator()(long n, long m) const
{
     if (m>V[n][0]) return 0;
     return V[n][m];
}

long& Graph::operator()(long n, long m)
{
     if (m>V[n][0]) 
     {
	  printf("Error accessing Graph data!\n");
	  exit(1);
     }
     return V[n][m];
}

long Graph::Degree(long p) const// number of neighbours of site p
{
     return V[p][0];
}

long Graph::Neighbour(long p, long k) const// k-th neighbour of site p
{
     return V[p][k];
}


void Graph::Update_Link_Index() 
// Puts to use the I vector of longs.
// I[n] means the number of links of all sites with index less than "n"
{
     updated=true;
     I=(long*)realloc(I,(N+1)*sizeof(long));
     I[1]=0;
     for (long i=2;i<=N;i++)
     {
	  I[i]=I[i-1];
	  for (long j=1;j<=V[i-1][0];j++)
	       if (V[i-1][j]>i-1) I[i]++; 
               // only links (s1-s2) with s2>s1 count
     }
}

long Graph::Get_Link_Index(long s1,long s2) const
// Get the index of the link corresponding to sites s1-s2
// Returns 0 if sites are not linked
{

     if (s1<=0 || s1>N || s2<=0 || s2>N) return 0;
     if (!updated) 
     {
	  merror("Link Index was not updated for Get_Link_Index!\n");
//	  Update_Link_Index();
//	  merror("Please, use first Update_Link_Index Function!!\n");
	  // exit(1);
     }

     long i1, i2;
     i1=(s1<s2 ? s1 : s2);
     i2=(s1<s2 ? s2 : s1);
     long index=I[i1]; 
     long j=0;
     do 
     {
	  j++; 
	  if (V[i1][j]>i1) index++;
     } while(V[i1][j]!=i2 && j<V[i1][0]);
     if (V[i1][j]!=i2) 
	  return 0;
     return index;
}

void Graph::Get_Link_Sites(long &s1, long &s2, long index) const
// returns s1 and s2 as the sites corresponding to link number "index"
{

     if (!updated)
     {
	  merror("Link Index was not updated for Get_Link_Sites!\n");
//	  Update_Link_Index();
	  
     }

     // a binary search
     long sa=1, sb=N;

#ifdef DEBUG
     if (!(I[sb]>=index && I[sa]<index)) 
	  merror("Check in Get_Link_Sites failed!\n");
#endif

     do
     {
	  long sc=(sb+sa)/2;
	  if (I[sc]>=index) sb=sc;
	  else sa=sc;
	  
     }while(sb-sa>1);
     s1=sa;

     long index_0=I[s1];
     long i=0, j=0; // s2 shall always be higher than s1
     do
     {
	  j++;
	  s2=V[s1][j];
	  if (s2>s1) i++;
     }while(index_0+i< index);
}

// Not quite efficient, but at least it works
// Should be optimized!
Graph Graph::Remove_Sites(const List &L) const
{
     Graph G(*this);
     for (long i=1;i<=L.N;i++)
	  G=G.Remove_Site(L(i));
     return G;
}

// Returning a graph perhaps is not the best idea!
// Let's see how it develops...
Graph Graph::Remove_Site(long n) const
{
     Graph G(*this);
     for (long k=1;k<=V[n][0];k++)
	  G.Remove_Link(n,V[n][k]);
     return G;
}
	  
long Graph::Find_Longest_Path_From(List &P, long is) const
{
     if (V[is][0]==0) 
     {
	  P.Append(is);
	  return P.N;
     }
     
     Graph Gis=Remove_Site(is);
     List Q, Q2;
     
     long lmax=0;
     for (long j=1;j<=V[is][0];j++)
     {
	  Q.Append(is);
	  long l=Gis.Find_Longest_Path_From(Q,V[is][j]);
	  if (l>lmax)
	  {
	       lmax=l;
	       Copy(Q2,Q);
	  }
	  Q.Destroy();
     }
     P.Append(Q2);
     return P.N;
}

List Graph::Find_Longest_Path() const
{
     long lmax=0;
     List P, Q;
     
     for (long i=1;i<=N;i++)
     {
	  long l=Find_Longest_Path_From(Q,i);
	  if (l>lmax) 
	  { 
	       lmax=l; 
	       Copy(P,Q); 
	  }
	  Q.Destroy();
     }
     return P;
}

Matrix Graph::Adjacency_Matrix() const
{
     Matrix A(N);
     for (long i=1;i<=N;i++)
 	  for (long j=i+1;j<=N;j++)
	       if (Is_Link(i,j)) A(i,j)=A(j,i)=1.0;
     return A;
}

// void Booleanize(Matrix &M)
// {
//      for (long i=1;i<=M.N1;i++)
// 	  for (long j=1;j<=M.N2;j++)
// 	       if (M(i,j)!=0.0) M(i,j)=1.0;
// }

// List Graph::Connected_Component(long p) const
// {
//      List C;
//      Matrix A=Adjacency_Matrix();
//      Matrix A2=A;
//      Matrix M=Unit(N)+A;
//      for (long k=1;k<=N;k++)
//      {
// 	  A*=A2;
// 	  M+=A;
// 	  Booleanize(A);
// 	  Booleanize(M);
//      }
//      // now M=A+A^2+A^3+...
     
//      for (long i=1;i<=N;i++)
//      if (M(i,p)>0.0) C.Append(i);
//      return C;
// }

// Connected component not in list
// Auxiliary function to Connected_Component
// Add to L site p and all its connected sites which are not already in L
void CC_Not_In_List(List &L, const Graph &G, long p)
{
     L.Append(p);
     for (long k=1;k<=G(p);k++)
	  if (!L.Find(G(p,k))) CC_Not_In_List(L,G,G(p,k));
}

// Better way to find the connected component, we don't need the
// adjacency matrix!!!!
List Graph::Connected_Component(long p) const
{
     List L;
     CC_Not_In_List(L,*this,p);
     return L;
}

// Prim's (aka DJP) algorithm seems suitable
Graph Graph::Minimum_Spanning_Tree() const
{
     Graph G2(*this);
     G2.Clean(); // remove all links
   
     List Tree;
     List Links=List_Interval(1,Nl);
     
     long p=0;
     do
     {
	  p++;
     }while(!V[p][0]);
     Tree.Append(p);
     long s1, s2;
     long ausgang=0;
     long treen_old=0;
     do
     {
	  long i=1;
	  do
	  {
	       Get_Link_Sites(s1,s2,Links(i));
//	       printf("Sites %d and %d\n",s1,s2);
	       // If the link connects the graph and the rest?
	       if ((Tree.Find(s1) && !Tree.Find(s2)) ||
		   (!Tree.Find(s1) && Tree.Find(s2)))
	       {
		    //	    printf("Fulfill condition\n");
		    G2.New_Link(s1,s2);
		    if (Tree.Find(s1)) 
		    {
			 Tree.Append(s2);
			 // printf("Appending %d to Tree\n",s2);
		    }
		    else
		    {
			 Tree.Append(s1);
			 //printf("Appending %d to Tree\n",s1); 
		    }
		    Links.Remove_All(Links(i));
	       }
	       i++;
	  }while(i<=Links.N);
	  // Check if all remaining links are longernal to the tree
	  long suma=0;
	  for (long i=1;i<=Links.N;i++)
	  {
	       Get_Link_Sites(s1,s2,Links(i));
	       if (Tree.Find(s1) && Tree.Find(s2)) suma++;
	  }
	  if (suma==Links.N) ausgang=1;
	  else if (Tree.N==treen_old)
	  {
	       printf("Graph was not connected, "
		      "Spanning tree applies to one component.\n");
	       ausgang=1;
	  }
	  treen_old=Tree.N;
     }while(!ausgang);
     return G2;
}


void Graph::Load(const char *name)
{
     FILE *fich=fopen(name,"rt");
     long n, nl;
     if (!fscanf(fich,"Graph with %ld sites and %ld links\n",&n,&nl))
	  merror("Error in Graph::Load");
     Destroy();
     for (long i=1;i<=n;i++)
     	  New_Site();
     for (long i=1;i<=n;i++)
     {
	  long i2, d;
	  if (!fscanf(fich,"%ld (%ld):",&i2,&d))
	       merror("Error in Graph::Load");
	  if (i2!=i) printf("Error reading graph\n");
	  V[i]=(long*)realloc(V[i],(d+1)*sizeof(long));
	  V[i][0]=d;
	  Nl+=d;
	  for (long j=1;j<=d;j++)
	       if (!fscanf(fich,"%ld",&V[i][j]))
		    merror("Error in Graph::Load");
     }
     Nl/=2;
     fclose(fich);
}

// Implementation of Dijkstra's algorithm
void Graph::Find_Path(List &P, long s1, long s2) const
{
     List D(N), Prev(N); // distances from s1
     for (long i=1;i<=N;i++)
     {
	  D(i)=LONG_MAX; // infinity
	  Prev(i)=0; // undefined
     }
     D(s1)=0; // s1 to s1 is zero...
     List S; // empty list, sites whose D(i) is OK
     List Q=List_Interval(1,N); // all sites, opposite of S
     while(Q.N) // while Q is not empty
     {
	  // Extract the site with the minimum D within Q
	  long imin=0, dmin=LONG_MAX;
	  for (long i=1;i<=Q.N;i++)
	       if (D(Q(i))<dmin) { imin=i; dmin=D(Q(i)); }
	  if (imin==0) printf("Help!\n");
	  long u=Q(imin);
	  S.Append(u);
	  Q.Remove_All(u);
	  for (long k=1;k<=Degree(u);k++)
	  {
	       long v=Neighbour(u,k);
	       if (D(v)>D(u)+1) // if links are weighted, this should be modif
	       {
		    D(v)=D(u)+1;
		    Prev(v)=u;
	       }
	  }
	  if (u==s2) break;
     }
     long p=s2;
     while(Prev(p))
     {
	  P.Append(p);
	  p=Prev(p);
     }
     P.Append(p);
}

long Graph::Distance(long s1, long s2) const
{
     List P;
     Find_Path(P,s1,s2);
     return P.N;
}

//////////////////////////////////////////////////////////////////////
// A few concrete graphs
//////////////////////////////////////////////////////////////////////

void Linear_Graph(Graph &G, long N)
{
     G.Destroy();
     G.New_Site();
     for (long i=2;i<=N;i++)
	 G.New_Site(i-1);
     G.Update_Link_Index();
}

void Linear_Graph_PBC(Graph &G, long N)
{
     G.Destroy();
     Linear_Graph(G,N);
     G.New_Link(1,N);
     G.Update_Link_Index();
}

void Square_Graph(Graph &G, long lx, long ly)
{
     G.Destroy();
     long N=lx*ly;
     G.Add_Sites(N);
    
     // graph index: (x-1)*Ly+y+1
     for (long x=1;x<=lx;x++)
         for (long y=1;y<=ly;y++)
	 {
	      if (x>1) G.New_Link((x-1)*ly+y,(x-2)*ly+y);
	      if (y>1) G.New_Link((x-1)*ly+y,(x-1)*ly+y-1);
	 }
     G.Update_Link_Index();

}

void Square_Graph_PBC(Graph &G, long lx, long ly)
{
     G.Destroy();
     long N=lx*ly;
     G.Add_Sites(N);

     // graph index: (x-1)*Ly+y+1
     for (long x=1;x<=lx;x++)
         for (long y=1;y<=ly;y++)
	 {
	      if (x>1) G.New_Link((x-1)*ly+y,(x-2)*ly+y);
	      else G.New_Link((x-1)*ly+y,(ly-1)*ly+y);
	      
	      if (y>1) G.New_Link((x-1)*ly+y,(x-1)*ly+y-1);
	      else G.New_Link((x-1)*ly+y,(x-1)*ly+ly);
	     
	 }
     G.Update_Link_Index();
}

void Complete_Graph(Graph &G, long N)
{
     G.Destroy();
     G.Add_Sites(N);
     for (long i=1;i<N;i++)
	  for (long j=i+1;j<=N;j++)
	       G.New_Link(i,j);
     G.Update_Link_Index();
}

// Implementation of Dijkstra's algorithm
// V is the vector of distances for each link
List Find_Path(const Graph &G, const Vector &V, long s1, long s2)
{
     List P;
     long N=G.N;
     Vector D(N); 
     List Prev(N); // distances from s1 
     for (long i=1;i<=N;i++)
     {
	  D(i)=1e20; // infinity
	  Prev(i)=0; // undefined
     }
     D(s1)=0.0; // s1 to s1 is zero...
     List S; // empty list, sites whose D(i) is OK
     List Q=List_Interval(1,N); // all sites, opposite of S
     while(Q.N) // while Q is not empty
     {
	  // Extract the site with the minimum D within Q
	  long imin=0; double dmin=1e20;
	  for (long i=1;i<=Q.N;i++)
	       if (D(Q(i))<dmin) { imin=i; dmin=D(Q(i)); }
	  if (imin==0) printf("Help!\n");
	  long u=Q(imin);
	  S.Append(u);
	  Q.Remove_All(u);
	  for (long k=1;k<=G.Degree(u);k++)
	  {
	       long v=G(u,k);
	       long jx=G.Get_Link_Index(u,v);
	       if (D(v)>D(u)+V(jx)) 
	       {
		    D(v)=D(u)+V(jx);
		    Prev(v)=u;
	       }
	  }
	  if (u==s2) break;
     }
     long p=s2;
     while(Prev(p))
     {
	  P.Append(p);
	  p=Prev(p);
     }
     P.Append(p);
     return P;
}

/////////////////////////////////////////////////////////////////
// Get a Graph from a coordinate Matrix and a cutoff distance
// Using usual euclidean distance, can be modified...
Graph Matrix_To_Graph(const Matrix &M)
{
     if (M.N1!=M.N2) merror("Need square Matrix for Matrix_To_Graph\n");
     long N=M.N1;
     Graph G(N);
     for (long i=1;i<=N;i++)
	  for (long j=1;j<=N;j++)
	       if (M(i,j)!=0.0) G.New_Link(i,j);
     G.Update_Link_Index();
     return G;
}

Graph Network_To_Graph(const Matrix &M, double dcutoff)
{
     long N=M.N2;
     Graph G(N);
     Matrix D(N);
     for (long i=1;i<=N;i++)
     {
	  Vector Vi=M.Col(i);
	  for (long j=i+1;j<=N;j++)
	  {
	       Vector Vj=M.Col(j);
	       double dist=(Vi-Vj).Norm();
	       if (dist<=dcutoff) G.New_Link(i,j);
	  }
     }
     G.Update_Link_Index();
     return G;
}


#endif
