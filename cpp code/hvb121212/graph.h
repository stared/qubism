////////////////////////////////////////////////////////
// hvb++ 1.0
// Copyleft: Javier Rodríguez Laguna
// 080725-101112

// Graph class

#ifndef GRAPH_HEADER
#define GRAPH_HEADER

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<limits.h>
#include"common.h"
#include"matrix.h"

////////////////////////////////////////////////////////////////////////////
////  HEADER
////////////////////////////////////////////////////////////////////////////

class Graph
{
 public:
     long N; // sites number
     long Nl; // number of links
     long **V; // neighborhood structure
     long *I; // indices for the links
     bool updated; // if the link indices were created
     
     Graph();
     Graph(long N);
     Graph(const Graph &);
     ~Graph();
     long Degree(long p) const; // number of neighbours of site p
     long Neighbour(long p, long k) const; // k-th neighbour of site p
     void Destroy();
     void Clean(); // delete all links, but retain basic structure
     bool New_Link(long, long);
     void Remove_Link(long, long);
     void Add_Sites(long n); // create n new sites
     long New_Site(long n=0); // n is the stem of the site
     bool Is_Link(long, long) const;
     void Write() const;
     Graph& operator=(const Graph &);

     long operator()(long) const;
     long& operator()(long,long);
     long operator()(long,long) const;

     Matrix Adjacency_Matrix() const;
     List Connected_Component(long) const;
     long Find_Longest_Path_From(List &P, long) const;
     List Find_Longest_Path() const;
     Graph Remove_Site(long) const;
     Graph Remove_Sites(const List &) const;
     Graph Minimum_Spanning_Tree() const;
     void Find_Path(List &P, long s1, long s2) const; // Dijkstra's algorithm
     long Distance(long s1, long s2) const;
     void Load(const char*);
     void Save(const char*) const;

     void Update_Link_Index();
     long Get_Link_Index(long,long) const;
     void Get_Link_Sites(long&,long&,long) const;
};

	 
void Copy(Graph&, const Graph&);

//////////////////////////////////////////////////////////////////////
// A few concrete graphs
void Linear_Graph(Graph &G, long N);
void Linear_Graph_PBC(Graph &G, long N);
void Square_Graph(Graph &G, long lx, long ly);
void Square_Graph_PBC(Graph &G, long lx, long ly);
void Complete_Graph(Graph &G, long N);
// void Hex_Graph(Graph &G, long lx, long ly);
// void Hex_Graph_PBC(Graph &G, long lx, long ly);

/////////////////////////////////////////////////////////////////
// A few handy routines
Graph Effective_Graph(const Graph &G, const List &L);
Graph Remove_Site(const Graph &G, long n);

/////////////////////////////////////////////////////////////////
// Useful algorithms
// Dijkstra's algorithm to find the shortest path between two
// sites of a graph with distances given by V.
List Find_Path(const Graph &G, const Vector &V, long s1, long s2);

/////////////////////////////////////////////////////////////////
// Get a Graph from a coordinate Matrix and a cutoff distance
// Using usual euclidean distance, can be modified...
Graph Matrix_To_Graph(const Matrix &M);
Graph Network_To_Graph(const Matrix &M, double dcutoff);

#endif





