/************************************************************
  fv_common.hpp
  polymer
************************************************************/

/***** functions for preprocess *****/
Atom2Int init_mass(){
  // Mass[a] is defined to be the floor of 10 times real mass of a
  // see http://www.inv.co.jp/~yoshi/kigou/gensi.html
  Atom2Int Mass;
  Mass["e*"] = 0;
  Mass["H"] = 10;
  Mass["B"]  = 108;
  Mass["C"]  = 120;
  Mass["O"]  = 160;
  Mass["N"]  = 140;
  Mass["F"]  = 190;
  Mass["Si"] = 280;
  Mass["P"]  = 309;
  Mass["S"]  = 320;
  Mass["Cl"] = 354;
  Mass["V"]  = 509;
  Mass["Br"] = 799;
  Mass["Cd"] = 1124;
  Mass["I"]  = 1269;
  Mass["Hg"] = 2006;
  Mass["Pb"] = 2072;
  Mass["Al"] = 269;
  return Mass;
}

/*****
// standard valence should not be assumed in the 2LMH model. 
Atom2Int init_valence(){
  Atom2Int Val;
  Val["e*"] = 2;
  Val["H"]  = 1;
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
*****/

/***** read SDF *****/
vector < ChemGraph > read_sdf(const string& fname){
  vector < ChemGraph > GraphSet;

  /* open the SDF file to be read */
  ifstream infile;
  try {
    infile.open(fname);
  } catch (const exception& e) {
    cerr << "error: failed to open '" << fname << "' for reading." << endl;
    throw(e);
  }
  if (! infile.is_open()) {
    cerr << "error: at read_sdf(): '" << fname << "' is not initialized." << endl;
    throw(-1);
  }

  /* read the SDF file named "fname" */
  string line;

  while (getline(infile, line)) {
    // prepare a new graph
    ChemGraph g;
    int n, m;
    stringstream st_stream;

    // read CID
    g.CID = line;
    getline(infile, line); // dummy
    getline(infile, line); // dummy
    
    // read n&m
    char *threeChars = new char[3];
    getline(infile, line);
    st_stream << line;
    st_stream.read(threeChars, 3);
    n = atoi(threeChars);
    st_stream.read(threeChars, 3);
    m = atoi(threeChars);
    delete[] threeChars;
    g.n = n;
    g.m = m;
    g.V.resize(n);
    
    g.is_polymer = false;

    if (g.CID[g.CID.size() - 1] == 'p'){
      g.is_polymer = true;
    }

    // read atom block (vertex)
    for (int i = 0; i < n; ++i) {
      stringstream str;
      string st;
      getline(infile, line);
      g.V[i].index = i; // index
      str << line;
      char *tdData = new char[31];
      str.read(tdData, 31); // read the 3D information from the line; ignored
      delete[] tdData;
      str >> g.V[i].atom; // atom
      str >> st; // dummy
      str >> st;

      switch(stoi(st)){ // eledeg
      case 0: g.V[i].eledeg = 0; break; 
      case 1: g.V[i].eledeg = 3; break; 
      case 2: g.V[i].eledeg = 2; break; 
      case 3: g.V[i].eledeg = 1; break; 
      case 5: g.V[i].eledeg = -1; break; 
      case 6: g.V[i].eledeg = -2; break; 
      case 7: g.V[i].eledeg = -3; break; 
      case 4:
	cerr << "error: " << g.CID << " has a radical atom. Please eliminate it. " << endl;
      default:
	cerr << "\t" << g.CID << " is strange. Halt." << endl;
	exit(1);
      }
      g.V[i].alpha = ""; 
      g.V[i].valence = 0;   // computed later
      g.V[i].bondsum = 0;   // computed later
      g.V[i].ht = Undef;    // computed in other function
      g.V[i].in_ex = Undef; // computed in other function
      g.V[i].edge.clear();  // computed later
      g.V[i].deg = 0;       // computed later
    }
    
    // read bond data
    for (int j = 0; j < m; ++j) {
      stringstream str;
      int v1, v2, mul;
      char *vChars = new char[3];
      getline(infile, line);
      str << line;
      str.read(vChars, 3);
      v1 = atoi(vChars); // id, 1-n of vertex
      str.read(vChars, 3);
      v2 = atoi(vChars); // id, 1-n of vertex
      str.read(vChars, 3);
      mul = atoi(vChars); // v1-v2 multiplicity
      delete[] vChars;

      // check feasibility
      if(mul==4){
	cerr << "error: aromatic edges cannot be treated. " << endl;
	cerr << "       please convert them into single/double-bond edges." << endl;
	exit(1);
      }
      
      // account for off-by-one in the indexing
      Edge e1,e2;
      e1.v = v2-1;
      e1.bond = mul;
      e1.in_ex = Undef; // computed in other function
      e1.link_edge = false;
      e1.connecting_edge = false;
      e2.v = v1-1;
      e2.bond = mul;
      e2.in_ex = Undef; // computed in other function
      e2.link_edge = false;
      e2.connecting_edge = false;

      if (j == m-1 && g.is_polymer){ // assume that the last edge in SDF file is a link edge and connectting edge
        e1.link_edge = true;
        e2.link_edge = true;

        e1.connecting_edge = true;
        e2.connecting_edge = true;
      }

      g.V[v1-1].edge.push_back(e1);
      g.V[v2-1].edge.push_back(e2);
    }
    
  
    // The end of the graph information in the file is marked with a "$$$$"
    while(getline(infile, line)){
      stringstream str;
      string st;
      if (line == "$$$$") {
	// compute deg, bondsum, valence and alpha of each vertex
	for (int i = 0; i < n; ++i){
	  g.V[i].deg = g.V[i].edge.size();
	  for(auto e:g.V[i].edge)
	    g.V[i].bondsum += e.bond;
	  g.V[i].valence = g.V[i].bondsum - g.V[i].eledeg;
	  g.V[i].alpha = g.V[i].atom + to_string(g.V[i].valence);
	  if(g.V[i].atom == "H" && g.V[i].valence>1){
	    cerr << "error: hydrogen with valence>1 is found." << endl;
	    cerr << "       probably it is out of scope." << endl;
	    exit(1);
	  }
	}
	GraphSet.push_back(g);
	break;
      }
      str << line;
      str >> st;
      if(st=="M"){
	str >> st;
	if(st=="CHG"){
	  int charges,i,chg;
	  str >> st;
	  charges = stoi(st);
	  for(int c=0;c<charges;c++){
	    str >> st;
	    i = stoi(st)-1;
	    str >> st;
	    chg = stoi(st);
	  }
	}
	else if(st=="RAD"){
	  cerr << "error: radical molecule is contained." << endl;
	  cerr << "       Remove it by preprocessor." << endl;
	  exit(1);
	}	
      }
    }
  }  
  return GraphSet;
}


ChemGraph suppress_H(const ChemGraph& G) {
	ChemGraph H;
	vector < int > iso;
	int *GtoH = new int[G.n];
	int n = 0;

	// CID of a generated graph H is identical to that of G
	H.CID = G.CID;
  H.is_polymer = G.is_polymer;

	// construct the mapping of vertex index from G to H
	for (int k = 0; k<G.n; k++)
		if (G.V[k].atom == "H")
			GtoH[k] = Undef;
		else
			GtoH[k] = true;

	for (int k = 0; k<G.n; k++)
		if (GtoH[k] != Undef)
			GtoH[k] = n++;

	H.n = n;
	H.m = 0;
	H.V.resize(n);

	for (int k = 0; k<G.n; k++) {
		if (GtoH[k] == Undef)
			continue;
		int i = GtoH[k];
		H.V[i].index = i;
		H.V[i].alpha = G.V[k].alpha;
		H.V[i].atom = G.V[k].atom;
		//H.V[i].bondsum = 0;
		H.V[i].bondsum = G.V[k].bondsum;
		H.V[i].eledeg = G.V[k].eledeg;
		H.V[i].ht = Undef;
		H.V[i].in_ex = Undef;
		H.V[i].edge.clear();
		for (int j = 0; j<G.V[k].deg; j++) {
			if (GtoH[G.V[k].edge[j].v] == Undef)
				continue;
			Edge e;
			e.v = GtoH[G.V[k].edge[j].v];
			e.bond = G.V[k].edge[j].bond;
			e.in_ex = Undef;
      e.link_edge = G.V[k].edge[j].link_edge;
      e.connecting_edge = G.V[k].edge[j].connecting_edge;
			H.V[i].edge.push_back(e);
			//H.V[i].bondsum += e.bond;      
			if (e.bond < 1 || e.bond > 4) {
				cerr << "error at suppress_H" << endl;
				exit(1);
			}
		}
		H.V[i].valence = H.V[i].bondsum - H.V[i].eledeg;
		H.V[i].deg = H.V[i].edge.size();
		H.m += H.V[i].deg;
}

	H.m = H.m / 2;

#ifdef DEBUG
	cout << H.CID << endl;
	for (int i = 0; i<H.n; i++) {
		cout << i << "\t" << H.V[i].alpha << " -->";
		for (int j = 0; j<H.V[i].deg; j++)
			cout << " " << H.V[i].edge[j].v << "(" << H.V[i].edge[j].bond << ")";
		cout << "\n";
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
    if(G.V[i].atom=="H"){
      cerr << "error: compute_ht_in_ex should be used for H-suppressed graph." << endl;
      exit(1);
    }
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

void get_rank(const vector < ChemGraph >& GS,
        vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
  for(int k=0;k<num_graphs;k++){
    int rank;
    rank = GS[k].m - (GS[k].n-1);
    if(rank<0){
      cerr << "error: CID "<< GS[k].CID << "is strange; n=" << GS[k].n << ", m=" << GS[k].m << ", rank="<< rank << endl;
      exit(EXIT_FAILURE);
    }
    ftr_val[k].push_back(rank);
  }
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

void get_num_link_edge(const vector < ChemGraph >& GS,
             vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
  for(int k=0;k<num_graphs;k++){
    int num_link_edge = 0;
    for(int i=0;i<GS[k].n;i++){
      for (auto& e : GS[k].V[i].edge){
        if (e.link_edge) ++num_link_edge;
      }
    }
    num_link_edge /= 2;
    ftr_val[k].push_back(num_link_edge);
  }
}


void get_avg_mass(Atom2Int& Mass,
		  const vector < ChemGraph >& GS,
		  vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
  for(int k=0;k<num_graphs;k++){
    double avg = 0.;
    for(int i=0;i<GS[k].n;i++)
      avg += (double)Mass[GS[k].V[i].atom];
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
      int d;
      for(auto e:v.edge)
	bond_sum += e.bond;
      d = val-bond_sum;
      if(d<0 || d>3){
	cerr << "error: val=" << val << ", bond_sum=" << bond_sum << " at atom "<< v.alpha << "(" << v.index << ") of CID=" << GS[k].CID << endl;
	exit(1);
      }
      dist[d]++;
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
	if(a>=0 && a<num_elems) // this is needed since this element may not be contained in Lambda when GS is a test set
	  freq[a]++;
      }
    for(int a=0;a<num_elems;a++)
      ftr_val[k].push_back(freq[a]);	      
  }

  delete[] freq;
}


string get_EC_string(Atom2Int &Mass,
	const string atom_0, const int deg_0, const int val_0,
	const string atom_1, const int deg_1, const int val_1,
	const int bond) {
	int prior = 0;
	if (Mass[atom_0] > Mass[atom_1])
		prior = 1;
	else if (Mass[atom_0] == Mass[atom_1]) {
		if (val_0 > val_1)
			prior = 1;
		else if (val_0 == val_1) {
			if (deg_0 > deg_1)
				prior = 1;
		}
	}

#ifdef DEBUG
	if (prior == 0)
		cout << atom_0 + to_string(val_0) + "_" + to_string(deg_0) + "_"
		+ atom_1 + to_string(val_1) + "_" + to_string(deg_1) + "_"
		+ to_string(bond) << endl;

	// prior==1 case:
	cout << atom_1 + to_string(val_1) + "_" + to_string(deg_1) + "_"
		+ atom_0 + to_string(val_0) + "_" + to_string(deg_0) + "_"
		+ to_string(bond) << endl;

#endif


	if (prior == 0)
		return atom_0 + to_string(val_0) + "_" + to_string(deg_0) + "_"
		+ atom_1 + to_string(val_1) + "_" + to_string(deg_1) + "_"
		+ to_string(bond);

	// prior==1 case:
	return atom_1 + to_string(val_1) + "_" + to_string(deg_1) + "_"
		+ atom_0 + to_string(val_0) + "_" + to_string(deg_0) + "_"
		+ to_string(bond);
}


string get_CS_string(Atom2Int &Mass,
	const string atom_0, const int deg_0, const int val_0) {

	return atom_0 + to_string(val_0) + "_" + to_string(deg_0);
}



Atom2Int get_exist_EC(Atom2Int& Mass, const vector < ChemGraph >& GS,
		      const int in_ex){
  Atom2Int EC;
  EC.clear();
  for(auto G:GS)
    for(auto v:G.V)
      for(auto e:v.edge)
	if(v.index < e.v && e.in_ex == in_ex){
	  VertexStruct u = G.V[e.v];
	  string key = get_EC_string(Mass, v.atom, v.deg, v.valence, u.atom, u.deg, u.valence, e.bond);
	  if(EC.find(key)==EC.end())
	    EC[key] = 0;
	  EC[key]++;
	}
  return EC;
}

Atom2Int get_exist_EC_link_edge(Atom2Int& Mass, const vector < ChemGraph >& GS,
          const int in_ex){
  Atom2Int EC;
  EC.clear();
  for(auto G:GS)
    for(auto v:G.V)
      for(auto e:v.edge)
  if(e.link_edge && v.index < e.v && e.in_ex == in_ex){
    VertexStruct u = G.V[e.v];
    string key = get_EC_string(Mass, v.atom, v.deg, v.valence, u.atom, u.deg, u.valence, e.bond);
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
	  VertexStruct u = GS[k].V[e.v];
	  string key = get_EC_string(Mass, v.atom, v.deg, v.valence, u.atom, u.deg, u.valence, e.bond);
	  int a=0;
	  for(auto itr=EC_in.begin();itr!=EC_in.end();itr++){
	    if(key == itr->first)
	      break;
	    a++;
	  }
	  if(a>=0 && a<num_ECs) // this is needed since this EC may not be contained in EC_in when GS is a test set
	    freq[a]++;
	}
    for(int a=0;a<num_ECs;a++)
      ftr_val[k].push_back(freq[a]);	      
  }

  delete[] freq;  
  
}

void get_EC_link_edge(Atom2Int& Mass, Atom2Int& EC_in, const vector < ChemGraph >& GS,
      const int in_ex, vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
  int num_ECs = EC_in.size();
  int *freq = new int[num_ECs];

  for(int k=0;k<num_graphs;k++){
    for(int a=0;a<num_ECs;a++)
      freq[a] = 0;
    for(auto v:GS[k].V)
      for(auto e:v.edge)
  if(e.link_edge && v.index < e.v && e.in_ex == in_ex){
    VertexStruct u = GS[k].V[e.v];
    string key = get_EC_string(Mass, v.atom, v.deg, v.valence, u.atom, u.deg, u.valence, e.bond);
    int a=0;
    for(auto itr=EC_in.begin();itr!=EC_in.end();itr++){
      if(key == itr->first)
        break;
      a++;
    }
    if(a>=0 && a<num_ECs) // this is needed since this EC may not be contained in EC_in when GS is a test set
      freq[a]++;
  }
    for(int a=0;a<num_ECs;a++)
      ftr_val[k].push_back(freq[a]);        
  }

  delete[] freq;  
  
}

Atom2Int get_exist_LeafAC(const vector < ChemGraph >& GS){
  Atom2Int LeafAC;
  LeafAC.clear();
  for(auto G:GS)
    for(auto v:G.V){
      if(v.deg != 1)
  continue;
      Edge e = v.edge[0];
      string key = v.alpha + "_" + G.V[e.v].alpha + "_" + to_string(e.bond); 
      if(LeafAC.find(key)==LeafAC.end())
  LeafAC[key] = 0;
      LeafAC[key]++;
    }
  return LeafAC;
}


void get_LeafAC(Atom2Int& AC_leaf,
    const vector < ChemGraph >& GS,
    vector < vector < double  > >& ftr_val){
  int num_graphs = GS.size();
  int num_ACs = AC_leaf.size();
  int *freq = new int[num_ACs];
  for(int k=0;k<num_graphs;k++){
    for(int a=0;a<num_ACs;a++)
      freq[a] = 0;
    for(auto v:GS[k].V){
      if(v.deg != 1)
  continue;
      Edge e = v.edge[0];
      string key = v.alpha + "_" + GS[k].V[e.v].alpha + "_" + to_string(e.bond); 
      int a=0;
      for(auto itr=AC_leaf.begin();itr!=AC_leaf.end();itr++){
  if(key == itr->first)
    break;
  a++;
      }
      if(a>=0 && a<num_ACs) // this is needed since this AC may not be contained in AC_leaf when GS is a test set
  freq[a]++;
    }
    for(int a=0;a<num_ACs;a++)
      ftr_val[k].push_back(freq[a]);        
  }
  delete[] freq;  
}


Atom2Int get_exist_CS_connecting_vertices(Atom2Int& Mass, const vector < ChemGraph >& GS,
          const int in_ex){
  Atom2Int CS;
  CS.clear();
  for(auto G:GS)
    for(auto v:G.V)
      for(auto e:v.edge)
  if(e.connecting_edge && v.index < e.v && e.in_ex == in_ex){
    VertexStruct u = G.V[e.v];
    string key_1 = get_CS_string(Mass, v.atom, v.deg, v.valence);
    if(CS.find(key_1)==CS.end())
      CS[key_1] = 0;
    CS[key_1]++;
    string key_2 = get_CS_string(Mass, u.atom, u.deg, u.valence);
    if(CS.find(key_2)==CS.end())
      CS[key_2] = 0;
    CS[key_2]++;
  }
  return CS;
}


void get_CS_connecting_vertices(Atom2Int& Mass, Atom2Int& CS, const vector < ChemGraph >& GS,
      const int in_ex, vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
  int num_CSs = CS.size();
  int *freq = new int[num_CSs];

  for(int k=0;k<num_graphs;k++){
    for(int a=0;a<num_CSs;a++)
      freq[a] = 0;
    for(auto v:GS[k].V)
      for(auto e:v.edge)
  if(e.connecting_edge && v.index < e.v && e.in_ex == in_ex){
    VertexStruct u = GS[k].V[e.v];
    string key_1 = get_CS_string(Mass, v.atom, v.deg, v.valence);
    int a=0;
    for(auto itr=CS.begin();itr!=CS.end();itr++){
      if(key_1 == itr->first)
        break;
      a++;
    }
    if(a>=0 && a<num_CSs) // this is needed since this CS may not be contained in CS when GS is a test set
      freq[a]++;
    
    string key_2 = get_CS_string(Mass, u.atom, u.deg, u.valence);
    int b=0;
    for(auto itr=CS.begin();itr!=CS.end();itr++){
      if(key_2 == itr->first)
        break;
      b++;
    }
    if(b>=0 && a<num_CSs) // this is needed since this CS may not be contained in CS when GS is a test set
      freq[b]++;
  }
    for(int a=0;a<num_CSs;a++)
      ftr_val[k].push_back(freq[a]);    
  }

  delete[] freq;  
  
}

void add_electric_density(const vector < ChemGraph >& GS, vector < vector < double > >& ftr_val){
  int num_graphs = GS.size();
  for(int k=0;k<num_graphs;k++){
    string electric_density = GS[k].CID;

    vector<string> v;
    stringstream ss{electric_density};
    string buf;
    while (getline(ss, buf, '_')) {
      v.push_back(buf);
    }

    int _electric_density = stoi(v[v.size()-1]);

    ftr_val[k].push_back(_electric_density);
  }
}




/***** functions for transformation *****/
void normalize( vector < vector < double > >& ftr_val,
		double *min_val, double *max_val ){
  int n,m;
  n = ftr_val.size();
  m = ftr_val[0].size();
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      if(min_val[j] == max_val[j])
	ftr_val[i][j] = 0.;
 //      else if(ftr_val[i][j] > max_val[j])
	// ftr_val[i][j] = 1.;
 //      else if(ftr_val[i][j] < min_val[j])
	// ftr_val[i][j] = 0.;
      else
	ftr_val[i][j] = (ftr_val[i][j]-min_val[j]) /
	  (max_val[j]-min_val[j]);
  }
}

/***** functions for link edge (2LMM_polymer) *****/
void find_C(vector <int>& C_already, const ChemGraph& G, bool C_found, vector <int>& C){
  if (C_found){
    return;
  }

  int u = C_already[C_already.size() - 1];

  if (C_already[0] == u){
    C = C_already;
    C_found = true;
    return;
  }

  for (auto& E : G.V[u].edge){
    if ((C_already.size() > 2 && E.v == C_already[0]) || 
      find(C_already.begin(), C_already.end(), E.v) == C_already.end()
    ){
      auto C_new = C_already;
      C_new.push_back(E.v);
      find_C(C_new, G, C_found, C);
      if (C_found) break;
    }
  }
}

void find_maximal_component(
  set <int>& comps, const int& v1, const int& v2, const int& _v1, const int& _v2, const ChemGraph& G
){
  // find the maximal componenet containing vertex v1 when ignoring the edge (v1, v2) and (_v1, _v2)
  comps.clear();
  map <int, bool> reached;
  for (auto& v : G.V){
    reached.emplace(v.index, false);
  }

  deque <int> qq;
  qq.push_back(v1);
  // comps.push_back(v1);

  while (!qq.empty()){
    int s = qq.back();
    qq.pop_back();
    if (reached[s]) continue;
    reached[s] = true;
    comps.insert(s);

    for (auto& e : G.V[s].edge){
      if (!reached[e.v] && 
        (s != v1 || e.v != v2) && (s != v2 || e.v != v1) &&
        (s != _v1 || e.v != _v2) && (s != _v2 || e.v != _v1)
      ){
        qq.push_back(e.v);
      }
    }
  }
}

void find_link_edge(vector <ChemGraph>& GS){

  for (auto& G : GS){
    if (!G.is_polymer){
      continue;
    }

    // 1. find the link edge e1 = (v1, v2) that already known
    int v1, v2;
    bool flag = false;
    for (auto& V : G.V){
      for (auto& E : V.edge){
        if (E.link_edge){
          v1 = V.index;
          v2 = E.v;
          flag = true;
          break;
        }
      }
      if (flag) break;
    }

    // 2. find a circle C that contains e1
    vector <int> C;
    vector <int> C_already({v1, v2});

    int num_v2_appear = 0;
    for (auto& e : G.V[v1].edge){
      if (e.v == v2){
        ++num_v2_appear;
      }
    }
   // cout << num_v2_appear << endl;
    if (num_v2_appear == 1) {
      find_C(C_already, G, false, C);
    } else {
      // multi edge case
      C.push_back(v1);
      C.push_back(v2);
      C.push_back(v1);
    }

    // 3. check every edge in C whether it is a link edge (check connectivity after ignoring e1 and itself)
    for (int i = 1; i < C.size() - 1; ++i){
      int _v1 = C[i];
      int _v2 = C[i + 1];

      set <int> comps;
      find_maximal_component(comps, v1, v2, _v1, _v2, G);

      if (comps.size() < G.n){
        // edge (_v1, _v2) is a link edge
        for (auto& E : G.V[_v1].edge){
          if (E.v == _v2){
            E.link_edge = true;
            break;
          }
        }
        for (auto& E : G.V[_v2].edge){
          if (E.v == _v1){
            E.link_edge = true;
            break;
          }
        }
      }
    }

  }
}
