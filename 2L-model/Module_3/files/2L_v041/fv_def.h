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
typedef struct _Vertex Vertex;
typedef struct _Edge Edge;
typedef struct _ChemGraph ChemGraph;
typedef map<string, int> Atom2Int;

/********** structures **********/
struct _Param{
  string sdf_file;
  string csv_file;
  string fringe_file;
  bool norm=false;
  bool std=false;
};

struct _Vertex{
  int index;            // index
  string alpha;         // atom
  int ht;               // height
  int in_ex;            // whether interior or exterior
  int deg;              // degree
  vector < Edge > edge; // list of connecting edges
};

struct _Edge{
  int v;      // index of connected vertex
  int bond;   // bond multiplicity
  int in_ex;  // whether interior or exterior
};

struct _ChemGraph{
  string CID;          // CID
  int n;               // number of vertices
  int m;               // number of edges
  vector < Vertex > V; // vector of vertices
};
