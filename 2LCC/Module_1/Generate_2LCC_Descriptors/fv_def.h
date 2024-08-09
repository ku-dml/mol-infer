/************************************************************
  fv_def.h
************************************************************/

/********** include **********/
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <queue>
using namespace std;

/********** macros **********/
#define DEBUG
#undef DEBUG

#define RHO 2
#define Undef -1
#define INTERIOR 0
#define EXTERIOR 1

/********** typedefs **********/
typedef struct _Param Param;
typedef struct _Vertex VertexStruct;
typedef struct _Edge Edge;
typedef struct _ChemGraph ChemGraph;
typedef map<string, int> Atom2Int;

/********** structures **********/
struct _Param{
  string sdf_file;
  string csv_file;
  string norm_csv_file;
  string fringe_file;
  bool   test_use;
  string test_sdf;
  string test_csv;
  string test_norm_csv;
  bool norm=false;
  bool std=false;
};

struct _Vertex{
  int index;      // index
  string alpha;   // symbol of chemical element = atom+valence
  string atom;    // atom
  int valence;    // valence (defined as bondsum-eledeg)
  int bondsum;    // sum of bonds over the edges incident to the vertex
  int eledeg;     // eledeg (read from SDF)
  int ht;         // height
  int in_ex;      // whether interior or exterior; undefined for hydrogen
  int deg;        // degree
  vector < Edge > edge; // list of connecting edges
};

struct _Edge{
  int v;      // index of connected vertex
  int bond;   // bond multiplicity
  int in_ex;  // whether interior or exterior; undefined for edge incident to H
};

struct _ChemGraph{
  string CID;          // CID
  int n;               // number of vertices
  int m;               // number of edges
  vector < VertexStruct > V; // vector of vertices
};


/********** global variables **********/
Atom2Int Mass;
