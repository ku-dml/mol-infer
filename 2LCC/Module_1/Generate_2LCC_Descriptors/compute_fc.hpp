#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>

#include "./fringe/TopologyGraph.hpp"
#include "./fringe/RootedGraph.hpp"

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
int fc_main_v03(string sdf_file,
		string fringe_file,
		string test_sdf, 
		int rho,
		vector < string >& fc_CID,
		vector < string >& fc_ftr_name,
		vector < vector < double > >& fc_ftr_val,
		vector < vector < double > >& fc_test_val){
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
    ChemicalGraph _g; // must with H
    int cond = _g.readFromFile(infile, false, true);
    if (cond == 0){
      break;
    }
    if (cond == -2){
      continue;
    }
    _g.atom_name();

    auto g = _g;     // g is H-fully-attached 
    _g.HSuppress();  // _g is H-suppressed
    
    int n = _g.EffectiveAtomNum();
    if (n == 0){
      continue;
    }
    if (maxd_check != 0 && _g.maxd() > 4){
      continue;
    }
    // auto g = RemoveNullVertex(_g);
        
    vector <RootedTree> FTs = getFringeTrees_H(g, _g, rho);

    unordered_map <TreeSignature, int> TS_map_tmp;
    TS_map_tmp.clear();

    for (auto& RT : FTs){
      ++total_fringe_tree_num;
      auto RT_tmp = RemoveNullVertex(RT);
      RootedGraph RG(RT_tmp);
      RootedMultiTree RMT(RG, RT_tmp);
      for (Vertex u = 0; u < RMT.n; ++u){
	//RMT.label[u] = atom_map.at(RT_tmp.graph.alpha[u]);
		  int len = RT_tmp.graph.alpha[u].size();
		  int ind = 0;
		  for (int i = 0; i < len; ++i) {
			  if (RT_tmp.graph.alpha[u][i] >= '0' && RT_tmp.graph.alpha[u][i] <= '9') {
				  ind = i;
				  break;
			  }
		  }
		  string str_tmp = RT_tmp.graph.alpha[u].substr(0, ind);
		  int val_tmp = stoi(RT_tmp.graph.alpha[u].substr(ind));

		  RMT.label[u] = val_tmp * 200 + atom_map.at(str_tmp); 
		  // we are assigning colors as val(a)*200+atomic_number
      }
      TreeSignature TS = RMT.getSignature();
      ++TS_map[TS];     // TS_map stores whole f-trees in the data set
      ++TS_map_tmp[TS]; // TS_map_tmp stores f-trees for the current chem-graph
    }
    CID_TS_map[_g.CID] = TS_map_tmp;
    fc_CID.push_back(_g.CID); // 0. Store CID in order to check consistency
  }

  // 1. Construct queue in order to sort TreeSignatures
  priority_queue<TS_with_freq, vector<TS_with_freq>> Q;
  for(auto itr=TS_map.begin();itr!=TS_map.end();itr++)
    Q.push(TS_with_freq(itr->first,itr->second));

  // 2. Store frequency and determine feature names
  int num_f_trees=0;
  unordered_map <TreeSignature, int> Ordered_TS_map;
  vector < TreeSignature > Ordered_TS;
  ofstream fout; 
  if(fringe_file!=""){
    fout.open(fringe_file);
//    fout << "FID,(atom depth)-seq,Multiplicity seq, Charge seq" << endl;
  }

  while (!Q.empty()) {
    TS_with_freq T = Q.top();
    Ordered_TS_map[T.TS] = num_f_trees;
    Ordered_TS.push_back(T.TS);
    string fname = "FC_" + to_string(num_f_trees + 1);
    for (int i = 0; i < T.TS.delta.size(); i++) {
		if (i % 2 == 0) {
			int atom_num = T.TS.delta[i] % 200;
			int atom_val = (T.TS.delta[i] - atom_num) / 200;
			fname = fname + "_" + map_int2atomic[atom_num] + to_string(atom_val);
		}
      else
	fname = fname + "_" + to_string(T.TS.delta[i]);
    }
    fname = fname + "_/";
    for (int i = 0; i < T.TS.mu.size(); i++)
      fname = fname + "_" + to_string(T.TS.mu[i]);
	fname = fname + "_/";
	for (int i = 0; i < T.TS.chg.size(); i++)
		fname = fname + "_" + to_string(T.TS.chg[i]);    
	fc_ftr_name.push_back(fname);
    if(fringe_file!=""){
      fout << num_f_trees+1 << ",";
      for (int i = 0; i < T.TS.delta.size(); i++)
		  if (i % 2 == 0) {
			  int atom_num = T.TS.delta[i] % 200;
			  int atom_val = (T.TS.delta[i] - atom_num) / 200;
			  fout << " " << map_int2atomic[atom_num];
			  fout << atom_val;
		  }
	else
	  fout << " " << T.TS.delta[i];
      fout << ",";
      for (int i = 0; i < T.TS.mu.size(); i++)
	fout << " " << T.TS.mu[i];
	  fout << ",";
	  for (int i = 0; i < T.TS.chg.size(); i++)
		  fout << " " << T.TS.chg[i];
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
  infile.close();
  
  /********** processing of test set **********/
  if(test_sdf=="")
    return num_f_trees;

  vector < string > Test_CID;
  // read the sdf file
  try{
    infile.open(test_sdf);
  }catch (const exception& e){
    cerr << "Cannot open '" << sdf_file << "' for reading!" << endl;
    throw(e);
  }
  if(!infile.is_open()){
    cerr << "Error, infile not initialized!" << endl;
    throw(-1);
  }
  
  // preparation
  map <string, unordered_map <TreeSignature, int>> Test_CID_TS_map;
  Test_CID_TS_map.clear();
  
  // read chemical compounds one by one
  while (! infile.eof()){
    ChemicalGraph _g;
    int cond = _g.readFromFile(infile, false, true);
    if (cond == 0){
      break;
    }
    if (cond == -2){// ignoring RAD
      continue;
    }
    _g.atom_name();
    auto g = _g;     // g is H-fully-attached 
    _g.HSuppress();  // _g is H-suppressed
    
    int n = _g.EffectiveAtomNum();
    if (n == 0){
      continue;
    }
    if (maxd_check != 0 && _g.maxd() > 4){ // check part
      continue;
    }
        
    vector <RootedTree> FTs = getFringeTrees_H(g, _g, rho);
    unordered_map <TreeSignature, int> TS_map_tmp;
    TS_map_tmp.clear();

    for (auto& RT : FTs){
      auto RT_tmp = RemoveNullVertex(RT);
      RootedGraph RG(RT_tmp);
      RootedMultiTree RMT(RG, RT_tmp);
      for (Vertex u = 0; u < RMT.n; ++u){
	//RMT.label[u] = atom_map.at(RT_tmp.graph.alpha[u]);
	int len = RT_tmp.graph.alpha[u].size();
	int ind = 0;
	for (int i = 0; i < len; ++i) {
		if (RT_tmp.graph.alpha[u][i] >= '0' && RT_tmp.graph.alpha[u][i] <= '9') {
			ind = i;
			break;
		}
	}
	string str_tmp = RT_tmp.graph.alpha[u].substr(0, ind);
	int val_tmp = stoi(RT_tmp.graph.alpha[u].substr(ind));

	RMT.label[u] = val_tmp * 200 + atom_map.at(str_tmp);
	  }
      TreeSignature TS = RMT.getSignature();
      ++TS_map_tmp[TS];
    }
    Test_CID_TS_map[_g.CID] = TS_map_tmp;
    Test_CID.push_back(_g.CID);
  }
   
  // 3. Compute frequncy for each chemical compound
  for(auto CID:Test_CID){
    vector < double > v;
    v.resize(num_f_trees);
    for(int j=0;j<num_f_trees;j++)
      v[j] = 0.;
    for(int j=0;j<num_f_trees;j++){
      TreeSignature TS = Ordered_TS[j];
      v[j] = Test_CID_TS_map[CID][TS];
    }
    fc_test_val.push_back(v);
  }

  return num_f_trees; 
}
