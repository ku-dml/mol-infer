/************************************************************
  Feature vector generator for 2L (two-layered) model 

  Copyright 2021
  Discrete Mathematics Laboratory at Kyoto University (KU-DML)
  Released Under the MIT License
  https://opensource.org/licenses/mit-license.php

 ************************************************************/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <queue>

#include "./fringe/TopologyGraph.hpp"
#include "./fringe/RootedGraph.hpp"

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
typedef struct _Vertex VertexStr;
typedef struct _Edge Edge;
typedef struct _ChemGraph ChemGraph;
typedef map<string, int> Atom2Int;

/********** structures **********/
struct _Param{
  string sdf_file;
  string output_stem;
  string desc_csv;
  string desc_norm_csv;
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
  vector < VertexStr > V; // vector of vertices
};


/***** functions for preprocess *****/
Atom2Int init_mass(){
  Atom2Int Mass;
  Mass["e*"] = 0;
  Mass["B"]  = 108;
  Mass["C"]  = 120;
  Mass["O"]  = 160;
  Mass["N"]  = 140;
  Mass["F"]  = 190;
  Mass["Si"] = 280;
  Mass["P"]  = 310;
  Mass["S"]  = 320;
  Mass["Cl"] = 355;
  Mass["V"]  = 510;
  Mass["Br"] = 800;
  Mass["Cd"] = 1124;
  Mass["I"]  = 1270;
  Mass["Hg"] = 2006;
  Mass["Pb"] = 2072;
  Mass["Al"] = 269;
  return Mass;
}

Atom2Int init_valence(){
  Atom2Int Val;
  Val["e*"] = 2;
  Val["B"]  = 3;
  Val["C"]  = 4;
  Val["O"]  = 2;
  Val["N"]  = 3;
  Val["F"]  = 1;
  Val["Si"] = 4;
  Val["P"]  = 3;
  Val["S"]  = 2;
  Val["Cl"] = 1;
  Val["Val"]  = 3;
  Val["Br"] = 1;
  Val["Cd"] = 2;
  Val["I"]  = 1;
  Val["Hg"] = 2;
  Val["Pb"] = 2;
  Val["Al"] = 3;
  return Val;
}

vector < ChemGraph > read_sdf(const string& fname){
  vector < ChemGraph > GraphSet;

  //open the file
  ifstream infile;
  try {
    infile.open(fname);
  } catch (const exception& e) {
    cerr << "Couldn't open '" << fname << "' for reading!" << endl;
    throw(e);
  }

  if (! infile.is_open()) {
    cerr << "Error, infile not initialized!" << endl;
    throw(-1);
  }

  //read the file
  int count = -1;
  int flag = 0;
  string line;

  while (getline(infile, line)) {

    if (flag == 0) {
      ChemGraph g;
      GraphSet.push_back(g);
      count++;
      flag = 1;
    }

    while (flag == 1) {
      GraphSet[count].CID = line;
      // Read and discard the next two lines of the file
      getline(infile, line);
      getline(infile, line);

      int n, m;
      stringstream st;
      char *threeChars = new char[4];
      threeChars[3] = '\0';
      getline(infile, line);
      st << line;
      st.read(threeChars, 3);
      n = atoi(threeChars);
      st.read(threeChars, 3);
      m = atoi(threeChars);
      delete[] threeChars;

      GraphSet[count].n = n;
      GraphSet[count].m = m;
      GraphSet[count].V.resize(n);
      
      // a local variable to store the atom symbol
      string alpha;
      // read atom data
      for (int i = 0; i < n; ++i) {
	getline(infile, line);
	stringstream st1;
	// temporary storage for the 3D data from a line in the sdf file
	char *tdData = new char[31];
	st1 << line;
	// read 3D coordinate data in a temporary string storage
	st1.read(tdData, 31); // read the 3D information from the line
	st1 >> alpha;
	GraphSet[count].V[i].index = i;
	GraphSet[count].V[i].alpha = alpha;
	GraphSet[count].V[i].ht = Undef;
	GraphSet[count].V[i].in_ex = Undef;
	GraphSet[count].V[i].edge.clear();
	GraphSet[count].V[i].deg = 0;
	delete[] tdData;
      }

      // read bond data
      for (int j = 0; j < m; ++j) {
	// local variable for vertex ids
	int v1, v2, mul;
	getline(infile, line);
	stringstream st2;
	char *vChars = new char[4];
  vChars[3] = '\0';
	st2 << line;
	st2.read(vChars, 3);
	v1 = atoi(vChars); // id, 1-n of vertex
	st2.read(vChars, 3);
	v2 = atoi(vChars); // id, 1-n of vertex
	st2.read(vChars, 3);
	mul = atoi(vChars); // v1-v2 multiplicity
	delete[] vChars;
	
	// account for off-by-one in the indexing
	Edge e1,e2;
	e1.v = v2-1;
	e1.bond = mul;
	e1.in_ex = Undef;
	e2.v = v1-1;
	e2.bond = mul;
	e2.in_ex = Undef;
	GraphSet[count].V[v1-1].edge.push_back(e1);
	GraphSet[count].V[v2-1].edge.push_back(e2);
      }

      // compute degree (but will be modified for H-suppression)
      for (int i = 0; i < n; ++i)
	GraphSet[count].V[i].deg = GraphSet[count].V[i].edge.size();
	
      flag = 2;
    }

    // The end of the graph information in the file is marked with a "$$$$"
    if (line == "$$$$" && flag == 2) {
      flag = 0;
    }
  }
  return GraphSet;
}


ChemGraph suppress_H(const ChemGraph& G){
  ChemGraph H;
  vector < int > iso;
  int *GtoH = new int[G.n];
  int n=0;

  // CID of a generated graph H is identical to that of G
  H.CID = G.CID;

  // construct the mapping of vertex index from G to H
  for(int k=0;k<G.n;k++)
    if(G.V[k].alpha=="H")
      GtoH[k] = Undef;
    else
      GtoH[k] = true;

  // examine whether isolated points appear by suppression
  for(int k=0;k<G.n;k++){
    bool is_isolated = true;
    if(GtoH[k]==Undef)
      continue;
    for(auto e:G.V[k].edge)
      if(GtoH[e.v]==true){
	is_isolated = false;
	break;
      }
    if(is_isolated)
      iso.push_back(k);
  }

  // organize GtoH 
  for(auto k:iso)
    GtoH[k] = Undef;
  for(int k=0;k<G.n;k++)
    if(GtoH[k]!=Undef)
      GtoH[k] = n++;
    
  H.n = n;
  H.m = 0;
  H.V.resize(n);

  for(int k=0;k<G.n;k++){
    if(GtoH[k]==Undef)
      continue;
    int i=GtoH[k];
    H.V[i].index = i;
    H.V[i].alpha = G.V[k].alpha;
    H.V[i].ht = Undef;
    H.V[i].in_ex = Undef;
    H.V[i].edge.clear();
    for(int j=0;j<G.V[k].deg;j++){
      if(GtoH[G.V[k].edge[j].v]==Undef)
	continue;
      Edge e;
      e.v = GtoH[G.V[k].edge[j].v];
      e.bond = G.V[k].edge[j].bond;
      e.in_ex = Undef;
      H.V[i].edge.push_back(e);

      if(e.bond < 1 || e.bond > 4)
	exit(1);
    }
    H.V[i].deg = H.V[i].edge.size();
    H.m += H.V[i].deg;
  }

  H.m = H.m/2;

#ifdef DEBUG
  for(int i=0;i<H.n;i++){
    cout << i << "\t" << H.V[i].alpha << " -->";
    for(int j=0;j<H.V[i].deg;j++)
      cout << " " << H.V[i].edge[j].v << "(" << H.V[i].edge[j].bond << ")";
    cout <<"\n";
  }
  cout << H.n << " " << H.m << "\n";
#endif
  
  delete[] GtoH;
  return H;
}


ChemGraph compute_ht_in_ex(ChemGraph& G){
  int *tmpdeg = new int[G.n];
  int h = 0;

  // initialize ht and in_ex
  for(int i=0;i<G.n;i++){
    G.V[i].ht = Undef;
    G.V[i].in_ex = INTERIOR;
    tmpdeg[i] = G.V[i].deg;
  }

  // compute ht and in_ex of vertices
  while(true){
    bool is_all_assigned = true;
    vector < int > deg_one_vertices;
    deg_one_vertices.clear();
    for(int i=0;i<G.n;i++){
      if(G.V[i].ht!=Undef || tmpdeg[i]>1)
	continue;
      is_all_assigned = false;
      deg_one_vertices.push_back(i);
    }
    if(is_all_assigned)
      break;
    for(auto i:deg_one_vertices){
      G.V[i].ht = h;
      if(h<RHO)
	G.V[i].in_ex = EXTERIOR;
      for(auto e:G.V[i].edge)
	tmpdeg[e.v]--;
    }
    h++;
  }

  // compute in_ex of edges
  for(auto& v: G.V)
    for(auto& e: v.edge)
      if(v.in_ex==INTERIOR && G.V[e.v].in_ex==INTERIOR)
	e.in_ex = INTERIOR;
      else
	e.in_ex = EXTERIOR;
      
#ifdef DEBUG
  for(int i=0;i<G.n;i++){
    cout << i << " (" << G.V[i].ht << ":" << G.V[i].in_ex << ")";
    for(auto e:G.V[i].edge)
      cout << " --> " << e.v << "(" << G.V[e.v].ht << ":" << e.in_ex << ":" << e.bond << ")";
    cout << endl;
  }
#endif
  
  delete[] tmpdeg;
  return G;
}

/***** functions for postprocess *****/

// write data to csv file
void output_csv(const string csv_file,
		const vector < ChemGraph >& GS,
		const vector < string >& ftr_name,
		const vector < vector < double > >& ftr_val){

  int num_ftrs = ftr_name.size();
  ofstream fout(csv_file, ios::out);

  // header
  fout << "CID";
  for(int j=0;j<num_ftrs;j++)
    fout << "," << ftr_name[j];
  fout << endl;

  // output feature values
  for(int k=0;k<GS.size();k++){
    int num_vals = ftr_val[k].size();
    if(num_ftrs != num_vals){
      fprintf(stderr, "error: number of features (%d) != length of feature vectors (%d)\n", num_ftrs, num_vals);
      exit(1);
    }
    fout << GS[k].CID;
    for(int j=0;j<num_vals;j++)
      fout << "," << ftr_val[k][j];
    fout << endl;
  }
  fout.close();
}


/***** functions for computing feature values *****/
void get_num_vertices(const vector < ChemGraph >& GS,
		      vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
  for(int k=0;k<num_graphs;k++)
    ftr_val[k].push_back(GS[k].n);
}


void get_num_interior_vertices(const vector < ChemGraph >& GS,
			       vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
  for(int k=0;k<num_graphs;k++){
    int num_interiors = 0;
    for(int i=0;i<GS[k].n;i++)
      if(GS[k].V[i].in_ex==INTERIOR)
	num_interiors++;
    ftr_val[k].push_back(num_interiors);
  }
}


void get_avg_mass(Atom2Int& Mass,
		  const vector < ChemGraph >& GS,
		  vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
  for(int k=0;k<num_graphs;k++){
    double avg = 0.;
    for(int i=0;i<GS[k].n;i++)
      avg += (double)Mass[GS[k].V[i].alpha];
    avg = avg/(double)GS[k].n;
    ftr_val[k].push_back(avg);
  }
}


void get_interior_deg_G(const vector < ChemGraph >& GS,
			vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
  for(int k=0;k<num_graphs;k++){
    int num[4] = {0,0,0,0};
    for(auto v:GS[k].V)
      if(v.in_ex==INTERIOR)
	num[v.deg-1]++;
    for(auto d:num)
      ftr_val[k].push_back(d);
  }
}


void get_interior_deg_G_in(const vector < ChemGraph >& GS,
			   vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
  for(int k=0;k<num_graphs;k++){
    int num[4] = {0,0,0,0};
    for(auto v:GS[k].V){
      int deg = 0;
      if(v.in_ex==INTERIOR)
	for(auto e:v.edge)
	  if(GS[k].V[e.v].in_ex==INTERIOR)
	    deg++;
      if(deg>0)
	num[deg-1]++;
    }
    for(auto d:num)
      ftr_val[k].push_back(d);
  }
}


void get_interior_bond(const vector < ChemGraph >& GS,
		       vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
  for(int k=0;k<num_graphs;k++){
    int num[2] = {0,0};
    for(auto v:GS[k].V){
      if(v.in_ex==INTERIOR)
	for(auto e:v.edge){
	  if(v.index < e.v && e.in_ex==INTERIOR && e.bond==2)
	    num[0]++;
	  else if(v.index < e.v && e.in_ex==INTERIOR && e.bond==3)
	    num[1]++;
	}
    }
    for(auto d:num)
      ftr_val[k].push_back(d);
  }
}


void get_num_hydrogens(Atom2Int& Val,
		       const vector < ChemGraph >& GS,
		       vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
  for(int k=0;k<num_graphs;k++){
    int val_sum = 0;
    int bond_sum = 0;
    for(auto v:GS[k].V){
      val_sum += Val[v.alpha];
      for(auto e:v.edge)
	bond_sum += e.bond;
    }
    // CAUTION: bond of {u,v} is added to bond_sum from both u and v
    //          thus you should not double bond_sum 
    ftr_val[k].push_back(val_sum-bond_sum);
  }
}


void get_num_hydrodegs(Atom2Int& Val,
		       const vector < ChemGraph >& GS,
		       vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
	  
  for(int k=0;k<num_graphs;k++){
    int dist[4] = {0,0,0,0};
    for(auto v:GS[k].V){
      int val = Val[v.alpha];
      int bond_sum = 0;
      for(auto e:v.edge)
	bond_sum += e.bond;
      dist[val-bond_sum]++;
      //cout << v.alpha << " " << val << " " << bond_sum << endl;
    }
    for(int d=0;d<=3;d++)
      ftr_val[k].push_back(dist[d]);
  }
}


Atom2Int get_exist_elements(const vector < ChemGraph >& GS, int in_ex){
  Atom2Int F;
  F.clear();
  for(auto G:GS)
    for(auto v:G.V){
      if(v.in_ex != in_ex)
	continue;
      if(F.find(v.alpha) == F.end())
	F[v.alpha] = 0;
      F[v.alpha]++;      
    }
  return F;
}


void get_frequency(Atom2Int &Lambda, const vector < ChemGraph >& GS,
		   int in_ex, vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
  int num_elems = Lambda.size();
  int *freq = new int[num_elems];

  for(int k=0;k<num_graphs;k++){
    for(int a=0;a<num_elems;a++)
      freq[a] = 0;
    for(auto v:GS[k].V)
      if(v.in_ex==in_ex){
	int a=0;
	for(auto itr=Lambda.begin();itr!=Lambda.end();itr++){
	  if(v.alpha == itr->first)
	    break;
	  a++;
	}
	freq[a]++;
      }
    for(int a=0;a<num_elems;a++)
      ftr_val[k].push_back(freq[a]);	      
  }

  delete[] freq;
}

string get_EC_string(Atom2Int &Mass, const string a_0, const int deg_0,
		     const string a_1, const int deg_1, const int bond){
  string former_atom, latter_atom;
  int former_deg, latter_deg;
  if(a_0 == a_1){
    former_atom = a_0;
    latter_atom = a_1;
    if(deg_0 <= deg_1){
      former_deg = deg_0;
      latter_deg = deg_1;
    }
    else{
      former_deg = deg_1;
      latter_deg = deg_0;
    }
  }
  else if(Mass[a_0] < Mass[a_1]){
    former_atom = a_0;
    latter_atom = a_1;
    former_deg = deg_0;
    latter_deg = deg_1;
  }
  else{
    former_atom = a_1;
    latter_atom = a_0;
    former_deg = deg_1;
    latter_deg = deg_0;
  }
  return former_atom + to_string(former_deg) + "_" + latter_atom + to_string(latter_deg) + "_" + to_string(bond);
}


Atom2Int get_exist_EC(Atom2Int& Mass, const vector < ChemGraph >& GS,
		      const int in_ex){
  Atom2Int EC;
  EC.clear();
  for(auto G:GS)
    for(auto v:G.V)
      for(auto e:v.edge)
	if(v.index < e.v && e.in_ex == in_ex){
	  VertexStr u = G.V[e.v];
	  string key = get_EC_string(Mass, v.alpha, v.deg, u.alpha, u.deg, e.bond);
	  if(EC.find(key)==EC.end())
	    EC[key] = 0;
	  EC[key]++;
	}
  return EC;
}


void get_EC(Atom2Int& Mass, Atom2Int& EC_in, const vector < ChemGraph >& GS,
	    const int in_ex, vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
  int num_ECs = EC_in.size();
  int *freq = new int[num_ECs];

  for(int k=0;k<num_graphs;k++){
    for(int a=0;a<num_ECs;a++)
      freq[a] = 0;
    for(auto v:GS[k].V)
      for(auto e:v.edge)
	if(v.index < e.v && e.in_ex == in_ex){
	  VertexStr u = GS[k].V[e.v];
	  string key = get_EC_string(Mass, v.alpha, v.deg, u.alpha, u.deg, e.bond);
	  int a=0;
	  for(auto itr=EC_in.begin();itr!=EC_in.end();itr++){
	    if(key == itr->first)
	      break;
	    a++;
	  }
	  freq[a]++;
	}
    for(int a=0;a<num_ECs;a++)
      ftr_val[k].push_back(freq[a]);	      
  }

  delete[] freq;  
  
}

/***** the structure TS_with_freq is introduced to sort TS by freq *****/
class TS_with_freq{
public:
  TreeSignature TS;
  int freq;
  TS_with_freq(TreeSignature _TS, int _freq){
    this->TS = _TS;
    this->freq = _freq;
  };
};
  
bool operator<(const TS_with_freq& T1, const TS_with_freq& T2){
  return T1.freq < T2.freq;
};
/***********************************************************************/


void write_result_csv(
		      const string outputfilename, 
		      const unordered_map <TreeSignature, int>& TS_map, 
		      const int& total_fringe_tree_num
		      ){
    
  ofstream outputfile(outputfilename, ios::out);

  outputfile << "# non-isomorphic fringe trees = ," << TS_map.size() << "\n";
  outputfile << "# of all Fringe trees = ," << total_fringe_tree_num << "\n";
  outputfile << "\n";
  outputfile << "(atom depth)-seq, Multiplicity seq, Frequency, Percentage" << "\n";

  for (auto& tmp : TS_map){
    auto& TS = tmp.first;
    auto& freq = tmp.second;

    string str_tmp = "";
    int i = 0;
    for (auto& a : TS.delta){
      if (i % 2 == 0){
	str_tmp += map_int2atomic.at(a) + " ";
      } else {
	str_tmp += to_string(a) + " ";
      }
      ++i;   
    }
    outputfile << str_tmp << ",";
        
    str_tmp = "";
    for (auto& a : TS.mu){
      str_tmp += to_string(a) + " ";
    }
    outputfile << str_tmp << ",";
     
    outputfile << freq << "," << (double)freq / total_fringe_tree_num <<  "\n";      
  }

  outputfile.close();
}

void write_result(
		  const string outputfilename, 
		  const unordered_map <TreeSignature, int>& TS_map, 
		  const int& total_fringe_tree_num
		  ){
    
  ofstream outputfile(outputfilename, ios::out);

  outputfile << TS_map.size() << "\n";
  outputfile << total_fringe_tree_num << "\n";

  for (auto& tmp : TS_map){
    auto& TS = tmp.first;
    auto& freq = tmp.second;

    int i = 0;
    for (auto& a : TS.delta){
      if (i % 2 == 0){
	outputfile << map_int2atomic.at(a) << " ";
      } else {
	outputfile << a << " ";
      }
      ++i;   
    }
    outputfile << "\n";
    for (auto& a : TS.mu){
      outputfile << a << " ";
    }
    outputfile << "\n";
     
    outputfile << freq << " " << (double)freq / total_fringe_tree_num <<  "\n";      
  }

  outputfile.close();
}

void read_data(const string outputfilename, unordered_map <TreeSignature, int>& TS_map, int& total_fringe_tree_num){
  ifstream infile;
  string line;

  try {
    infile.open(outputfilename);
  } catch (const exception& e){
    return;
  }

  if (!infile.is_open()){
    return;
  }

  if (infile.eof()){
    return;
  }

  stringstream st1;
  getline(infile, line);
  st1 << line;
  int n;
  st1 >> n;
  st1.clear();
  getline(infile, line);
  st1 << line;
  st1 >> total_fringe_tree_num;
  st1.clear();

  for (int i = 0; i < n; ++i){
    TreeSignature TS;

    TS.delta.clear();
    TS.mu.clear();

    stringstream st2;
    getline(infile, line);
    st2 << line;

    string str_tmp;
    int tmp2;
    while (st2 >> str_tmp >> tmp2){
      TS.delta.push_back(map_atomic_number.at(str_tmp));
      TS.delta.push_back(tmp2);
    }

    stringstream st3;
    getline(infile, line);
    st3 << line;

    int tmp3;
    while (st3 >> tmp3){
      TS.mu.push_back(tmp3);
    }

    stringstream st4;
    getline(infile, line);
    int freq;
    st4 << line;
    st4 >> freq;

    TS_map.emplace(TS, freq);
  }
    
  infile.close();
}


/********** main of computing fc **********/
void fc_main_v02(string sdf_file,
		 string fringe_file, int rho,
		 vector < string >& fc_CID,
		 vector < string >& fc_ftr_name,
		 vector < vector < double > >& fc_ftr_val){
  int maxd_check = 1;
  ifstream infile;

  // read the sdf file
  try {
    infile.open(sdf_file);
  } catch (const exception& e){
    cerr << "Cannot open '" << sdf_file << "' for reading!" << endl;
    throw(e);
  }
  if (!infile.is_open()){
    cerr << "Error, infile not initialized!" << endl;
    throw(-1);
  }

  // preparation
  unordered_map <string, int> atom_map = map_atomic_number;
  unordered_map <TreeSignature, int> TS_map;
  TS_map.clear();
  map <string, unordered_map <TreeSignature, int>> CID_TS_map;
  CID_TS_map.clear();
  
  int total_fringe_tree_num = 0;

  // read chemical compounds one by one
  while (! infile.eof()){
    ChemicalGraph _g;
    int cond = _g.readFromFile(infile);
    if (cond == 0){
      break;
    }
    if (cond == -1){
      continue;
    }
    _g.HSuppress();
    int n = _g.EffectiveAtomNum();
    if (n == 0){
      continue;
    }
    if (maxd_check != 0 && _g.maxd() > 4){
      continue;
    }
    auto g = RemoveNullVertex(_g);
        
    vector <RootedTree> FTs = getFringeTrees(g, rho);

    unordered_map <TreeSignature, int> TS_map_tmp;
    TS_map_tmp.clear();

    for (auto& RT : FTs){
      ++total_fringe_tree_num;
      auto RT_tmp = RemoveNullVertex(RT);
      RootedGraph RG(RT_tmp);
      RootedMultiTree RMT(RG, RT_tmp);

      for (Vertex u = 0; u < RMT.n; ++u){
	RMT.label[u] = atom_map.at(RT_tmp.graph.alpha[u]);
      }
      TreeSignature TS = RMT.getSignature();
      ++TS_map[TS];
      ++TS_map_tmp[TS];
    }
    CID_TS_map[_g.CID] = TS_map_tmp;
    fc_CID.push_back(_g.CID); // 0. Store CID in order to check consistency
  }

  // 1. Construct queue in order to sort TreeSignatures
  priority_queue<TS_with_freq, vector<TS_with_freq>> Q;
  for(auto itr=TS_map.begin();itr!=TS_map.end();itr++)
    Q.push(TS_with_freq(itr->first,itr->second));

  // 2. Store frequency and determine feature names
  //    ... and output fringe-trees to fringe_file if it is specified (v02)
  int num_f_trees=0;
  unordered_map <TreeSignature, int> Ordered_TS_map;
  vector < TreeSignature > Ordered_TS;
  ofstream fout; 
  if(fringe_file!=""){
    fout.open(fringe_file);
    fout << "FID,(atom depth)-seq,Multiplicity seq" << endl;
  }

  
  while (!Q.empty()) {
    TS_with_freq T = Q.top();
    Ordered_TS_map[T.TS] = num_f_trees;
    Ordered_TS.push_back(T.TS);
    string fname = "FC_" + to_string(num_f_trees+1) + "_" + to_string(T.freq);
    fc_ftr_name.push_back(fname);
    if(fringe_file!=""){
      fout << num_f_trees+1 << ",";
      for (int i = 0; i < T.TS.delta.size(); i++)
	if(i%2==0)
	  fout << " " << map_int2atomic[T.TS.delta[i]];
	else
	  fout << " " << T.TS.delta[i];
      fout << ",";
      for (int i = 0; i < T.TS.mu.size(); i++)
	fout << " " << T.TS.mu[i];
      fout << endl;
    }
    num_f_trees++;
    Q.pop();
  }
  if(fringe_file!="")
    fout.close();
    
  // 3. Compute frequncy for each chemical compound
  for(auto CID:fc_CID){
    vector < double > v;
    v.resize(num_f_trees);
    for(int j=0;j<num_f_trees;j++)
      v[j] = 0.;
    for(int j=0;j<num_f_trees;j++){
      TreeSignature TS = Ordered_TS[j];
      v[j] = CID_TS_map[CID][TS];
    }
    fc_ftr_val.push_back(v);
  }

  
}


void normalize( vector < vector < double > >& ftr_val ){
  int n,m;
  n = ftr_val.size();
  m = ftr_val[0].size();
  for(int j=0;j<m;j++){
    double min=ftr_val[0][j];
    double max=ftr_val[0][j];
    for(int i=1;i<n;i++)
      if(ftr_val[i][j] < min)
        min = ftr_val[i][j];
      else if(ftr_val[i][j] > max)
        max = ftr_val[i][j];
    for(int i=0;i<n;i++)
      if(min==max)
        ftr_val[i][j] = 0.;
      else
        ftr_val[i][j] = (ftr_val[i][j]-min) / (max-min);
  }
}

/********** read arguments **********/
void read_arguments(int argc, char *argv[], Param& Par){
  if(argc!=3){
    fprintf(stderr, "usage: %s (input.sdf)(OUTPUT)\n\n", argv[0]);
    fprintf(stderr, "The program generates the following three files:\n");
    fprintf(stderr, "\t(1) OUTPUT_fringe.txt\n");
    fprintf(stderr, "\t(2) OUTPUT_desc.csv\n");
    fprintf(stderr, "\t(3) OUTPUT_desc_norm.csv\n\n");
    fprintf(stderr, "Module 2 requires only (3)\n");
    fprintf(stderr, "Module 3 requires all of (1), (2) and (3)\n");
    fprintf(stderr, "\n");
    exit(1);
  }  
  Par.sdf_file = argv[1];
  Par.output_stem = argv[2];
  Par.desc_csv = Par.output_stem + "_desc.csv";
  Par.desc_norm_csv = Par.output_stem + "_desc_norm.csv";
  Par.fringe_file = Par.output_stem + "_fringe.txt";
}


/********** main **********/
int main(int argc, char *argv[]){
  Param Par;
  vector < ChemGraph > Original,GS;
  Atom2Int Mass,Val;
  vector < string > ftr_name;
  vector < vector < double > > ftr_val;

  
  /***** read arguments *****/
  read_arguments(argc, argv, Par);
  
  /***** initialize *****/
  // initialize hashes for mass and valence
  Mass = init_mass();
  Val = init_valence();
  // read input.sdf
  Original = read_sdf(Par.sdf_file);
  // suppress hydrogens for our formulation
  for(auto G:Original)
    GS.push_back(suppress_H(G));
  // compute vertex height and decide if it is interior/exterior
  for(int k=0;k<GS.size();k++)
    GS[k] = compute_ht_in_ex(GS[k]);
  // resize ftr_val
  ftr_val.resize(GS.size());

  
  /***** obtain feature names & compute feature values *****/  
  // dcp_1(G): number n(G) of vertices in G
  ftr_name.push_back("n");
  get_num_vertices(GS, ftr_val);

  // dcp_2(G): number |V^int(G)| of interior vertices in G
  ftr_name.push_back("n_in");
  get_num_interior_vertices(GS, ftr_val);

  // dcp_3(G): average of mass* over all non-H atoms in G
  ftr_name.push_back("ms");
  get_avg_mass(Mass, GS, ftr_val);
  
  // dcp_i(G) for i=4 to 7: number dg_d(G) of interior vertices of deg d in G
  for(int d=1;d<=4;d++){
    ostringstream desc_name;
    desc_name << "dg_" << d;
    ftr_name.push_back(desc_name.str());
  }
  get_interior_deg_G(GS, ftr_val);
    
  // dcp_i(G) for i=8 to 11: number dg^int_d(G) of
  //                         interior vertices of deg d in G^int
  for(int d=1;d<=4;d++){
    ostringstream desc_name;
    desc_name << "dg_in_" << d;
    ftr_name.push_back(desc_name.str());
  }
  get_interior_deg_G_in(GS, ftr_val);
  
  // dcp_i(G) for i=12 to 15: number hydg_d(G) of vertices v in G of hydro-degree deg_hyd(v)=d
  for(int d=0;d<=3;d++){
    ostringstream desc_name;
    desc_name << "dg_hyd_" << d;
    ftr_name.push_back(desc_name.str());
  }
  get_num_hydrodegs(Val, GS, ftr_val);
  
  // dcp_i(G) for i=16 and 17: number bd^int_m of interior-edges with bond m
  ftr_name.push_back("bd_in_2");
  ftr_name.push_back("bd_in_3");
  get_interior_bond(GS, ftr_val);

  // dcp_i(G): number of chemical element in the interior vertices
  Atom2Int Lambda_in;
  Lambda_in = get_exist_elements(GS, INTERIOR);
  for(auto itr=Lambda_in.begin();itr!=Lambda_in.end();itr++)
    ftr_name.push_back("na_in_"+itr->first);
  get_frequency(Lambda_in, GS, INTERIOR, ftr_val);
  
  // dcp_i(G): number of chemical element in the exterior vertices
  Atom2Int Lambda_ex;
  Lambda_ex = get_exist_elements(GS, EXTERIOR);
  for(auto itr=Lambda_ex.begin();itr!=Lambda_ex.end();itr++)
    ftr_name.push_back("na_ex_"+itr->first);
  get_frequency(Lambda_ex, GS, EXTERIOR, ftr_val);

  // dcp_i(G): frequency of EC (edge-configurations) in interior edges
  Atom2Int EC_in;
  EC_in = get_exist_EC(Mass, GS, INTERIOR);
  for(auto itr=EC_in.begin();itr!=EC_in.end();itr++)
    ftr_name.push_back(itr->first);
  get_EC(Mass, EC_in, GS, INTERIOR, ftr_val);

  // dcp_i(G): frequency of FC (fringe-configurations) in exterior edges
  vector < string > fc_CID;
  vector < string > fc_ftr_name;
  vector < vector < double > > fc_ftr_val;
  fc_main_v02(Par.sdf_file, Par.fringe_file, RHO, fc_CID, fc_ftr_name, fc_ftr_val);
  for(int j=0;j<fc_ftr_name.size();j++)
    ftr_name.push_back(fc_ftr_name[j]);  
  for(int k=0;k<GS.size();k++){
    if(GS[k].CID != fc_CID[k]){
      cerr << "error: CID is not consistent (" << GS[k].CID << "!=" << fc_CID[k] << ")\n";
      exit(1);
    }
    for(auto x:fc_ftr_val[k])
      ftr_val[k].push_back(x);
  }
  cout << Par.fringe_file << " is generated successfully." << endl;

  
  /***** postprocess *****/
  output_csv(Par.desc_csv, GS, ftr_name, ftr_val);  
  cout << Par.desc_csv << " is generated successfully." << endl;
  normalize(ftr_val);
  output_csv(Par.desc_norm_csv, GS, ftr_name, ftr_val);  
  cout << Par.desc_norm_csv << " is generated successfully." << endl;
  return 0;
}
