/************************************************************
  fv_common.cpp
************************************************************/

#include "fv_def.h"
#include "fv_common.h"

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
      char *threeChars = new char[3];
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
	char *vChars = new char[3];
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
	  Vertex u = G.V[e.v];
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
	  Vertex u = GS[k].V[e.v];
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

