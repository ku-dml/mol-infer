/************************************************************
  Feature vector generator for 2LMM model

  2LMM: two-layered model with 
    - hydrogen;
    - multi-valence chemical elements; 
    - cations; and 
    - anions

  Ver: 1.0
  - first version that starts from 2L_v074

  (v011, v012) some bugs are fixed.
  (v013) fringe tree part is fixed.
  (v014) output is changed a bit. 
  (v015) significant bug in v014 is removed. 
  (v016) leaf_ac is implemented
  (v017) the order of descriptors is corrected. 
  (v018) rank is implemented

 ************************************************************/

/********** include header files **********/
#include "fv_def.h"
#include "fv_common.hpp"
#include "compute_fc.hpp"

#define N_INDEX 0

// main routine of computing fc 
int fc_main_v03(string sdf_file,
		string fringe_file,
		string test_sdf, 
		int rho, 
		vector < string >& fc_CID,
		vector < string >& fc_ftr_name,
		vector < vector < double > >& fc_ftr_val,
		vector < vector < double > >& fc_test_val);


/********** read arguments **********/
void read_arguments(int argc, char *argv[], Param& Par){
  if(argc!=3 && argc!=5){
    fprintf(stderr, "usage: %s (input.sdf)(OUTPUT)[(test.sdf)(TEST)]\n\n", argv[0]);
    fprintf(stderr, "Caution: in this model, input.sdf MUST BE completely H-defined.\n\n");
    fprintf(stderr, "The program outputs:\n");
    fprintf(stderr, "\tOUTPUT_desc.csv: unnormalized feature vectors\n");
    fprintf(stderr, "\tOUTPUT_desc_norm.csv: normalized feature vectors\n");
    fprintf(stderr, "\tOUTPUT_fringe.txt: text file that contains fringe trees\n");
    fprintf(stderr, "\nIf you specify (test.sdf) and (TEST),\n");
    fprintf(stderr, "TEST_desc.csv and TEST_desc_norm.csv are output accordingly.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  ostringstream ost;
  
  Par.sdf_file = argv[1];
  ost << argv[2] << "_desc.csv";
  Par.csv_file = ost.str();
  ost.str("");
  ost.clear();
  
  ost << argv[2] << "_desc_norm.csv";
  Par.norm_csv_file = ost.str();
  ost.str("");
  ost.clear();
  
  ost << argv[2] << "_fringe.txt";
  Par.fringe_file = ost.str();
  ost.str("");
  ost.clear();

  if(argc==3){
    Par.test_use = false;
    Par.test_sdf = "";
    Par.test_csv = "";
  }
  else{
    Par.test_use = true;
    Par.test_sdf = argv[3];
    ost << argv[4] << "_desc.csv";
    Par.test_csv = ost.str();
    ost.str("");
    ost.clear();
    ost << argv[4] << "_desc_norm.csv";
    Par.test_norm_csv = ost.str();
    ost.str("");
    ost.clear();
  }
    
}


/********** main **********/
int main(int argc, char *argv[]){
  Param Par;
  vector < ChemGraph > Train_with_H, Train_without_H;
  vector < ChemGraph > Test_with_H, Test_without_H;
  vector < string > ftr_name;
  vector < vector < double > > ftr_val, test_val;
  
  /***** read arguments *****/
  read_arguments(argc, argv, Par);
  
  /***** initialize *****/
  // initialize hashes for mass
  Mass = init_mass();

  /**** read input.sdf *****/
  Train_with_H = read_sdf(Par.sdf_file);
  ftr_val.resize(Train_with_H.size());
  // prepare H-suppressed version of training set
  for(auto G:Train_with_H)
    Train_without_H.push_back(suppress_H(G));
  for(int k = 0;k<Train_without_H.size(); k++)
    Train_without_H[k] = compute_ht_in_ex(Train_without_H[k]);  

  /***** read test set if necessary *****/
  if(Par.test_use){
    Test_with_H = read_sdf(Par.test_sdf);
    test_val.resize(Test_with_H.size());
    // prepare H-suppressed version of test set
    for(auto G:Test_with_H)
      Test_without_H.push_back(suppress_H(G));
    for(int k=0;k<Test_without_H.size();k++)
      Test_without_H[k] = compute_ht_in_ex(Test_without_H[k]);  
  }
  
  /***** obtain feature names & compute feature values *****/  
  // 1. dcp_i(G): number n(G) of non-hydrogen vertices in G
  ftr_name.push_back("n");
  get_num_vertices(Train_without_H, ftr_val);
  if(Par.test_use)
    get_num_vertices(Test_without_H, test_val);

  // 2. dcp_i(G): rank r(G) of G
  ftr_name.push_back("rank");
  get_rank(Train_with_H, ftr_val);
  if(Par.test_use)
    get_rank(Test_with_H, test_val);
  
  
  // 3. dcp_i(G): number |V^int(G)| of interior vertices in G
  ftr_name.push_back("n_in");
  get_num_interior_vertices(Train_without_H, ftr_val);
  if(Par.test_use)
    get_num_interior_vertices(Test_without_H, test_val);

  // 4. dcp_i(G): average of mass* over all atoms, including H, in G
  ftr_name.push_back("ms");
  get_avg_mass(Mass, Train_with_H, ftr_val);
  if(Par.test_use)
    get_avg_mass(Mass, Test_with_H, test_val);

  // 5. dcp_i(G): number dg_d(G) of interior vertices of deg d in G
  for(int d=1;d<=4;d++){
    ostringstream desc_name;
    desc_name << "dg_" << d;
    ftr_name.push_back(desc_name.str());
  }
  get_interior_deg_G(Train_without_H, ftr_val);
  if(Par.test_use)
    get_interior_deg_G(Test_without_H, test_val);

  // 6. dcp_i(G): number dg^int_d(G) of interior vertices of deg d in G^int
  for(int d=1;d<=4;d++){
    ostringstream desc_name;
    desc_name << "dg_in_" << d;
    ftr_name.push_back(desc_name.str());
  }
  get_interior_deg_G_in(Train_without_H, ftr_val);
  if(Par.test_use)
    get_interior_deg_G_in(Test_without_H, test_val);

  // 7. dcp_i(G): number bd^int_m of interior-edges with bond m
  ftr_name.push_back("bd_in_2");
  ftr_name.push_back("bd_in_3");
  get_interior_bond(Train_without_H, ftr_val);
  if(Par.test_use)
    get_interior_bond(Test_without_H, test_val);
  
  // 8. dcp_i(G): number of chemical element in the interior vertices
  Atom2Int Lambda_in;
  Lambda_in = get_exist_elements(Train_without_H, INTERIOR);
  for(auto itr=Lambda_in.begin();itr!=Lambda_in.end();itr++)
    ftr_name.push_back("na_in_"+itr->first);
  get_frequency(Lambda_in, Train_without_H, INTERIOR, ftr_val);
  if(Par.test_use)
    get_frequency(Lambda_in, Test_without_H, INTERIOR, test_val);
  
  // 9. dcp_i(G): number of chemical element in the exterior vertices
  Atom2Int Lambda_ex;
  Lambda_ex = get_exist_elements(Train_without_H, EXTERIOR);
  for(auto itr=Lambda_ex.begin();itr!=Lambda_ex.end();itr++)
    ftr_name.push_back("na_ex_"+itr->first);
  get_frequency(Lambda_ex, Train_without_H, EXTERIOR, ftr_val);
  if(Par.test_use)
    get_frequency(Lambda_ex, Test_without_H, EXTERIOR, test_val);

 
  // 10. dcp_i(G): frequency of EC (edge-configurations) in interior edges
  Atom2Int EC_in;
  EC_in = get_exist_EC(Mass, Train_without_H, INTERIOR);
  for(auto itr=EC_in.begin();itr!=EC_in.end();itr++)
    ftr_name.push_back(itr->first);
  get_EC(Mass, EC_in, Train_without_H, INTERIOR, ftr_val);
  if(Par.test_use)
    get_EC(Mass, EC_in, Test_without_H, INTERIOR, test_val);
  
  // 11. dcp_i(G): frequency of FC (fringe-configurations)
  vector < string > fc_CID;
  vector < string > fc_ftr_name;
  vector < vector < double > > fc_ftr_val, fc_test_val;
  int num_FC;
  num_FC = fc_main_v03(Par.sdf_file, Par.fringe_file, Par.test_sdf,
		       RHO, fc_CID, fc_ftr_name, fc_ftr_val, fc_test_val);
  for(int j=0;j<fc_ftr_name.size();j++)
    ftr_name.push_back(fc_ftr_name[j]);  
  for(int k=0;k<Train_with_H.size();k++){
    if(Train_with_H[k].CID != fc_CID[k]){
      cerr << "error: CID is not consistent (" << Train_with_H[k].CID << "!=" << fc_CID[k] << ")\n";
      exit(1);
    }
    for(auto x:fc_ftr_val[k])
      ftr_val[k].push_back(x);
  }
  if(Par.test_use)
    for(int k=0;k<Test_with_H.size();k++)
      for(auto x:fc_test_val[k])
	test_val[k].push_back(x);

  // 12. dcp_i(G): frequency of AC (adjacency-configurations) in leaf edges
  Atom2Int AC_leaf;
  AC_leaf = get_exist_LeafAC(Train_without_H);
  for(auto itr=AC_leaf.begin();itr!=AC_leaf.end();itr++)
    ftr_name.push_back("LeafAC_"+itr->first);
  get_LeafAC(AC_leaf, Train_without_H, ftr_val);
  if(Par.test_use)
    get_LeafAC(AC_leaf, Test_without_H, test_val);

  
  /***** output important stats *****/
  cout << "INPUT_SDF:\t" << Par.sdf_file << endl;
  cout << "NUM_VECTORS:\t" << Train_without_H.size() << endl;
  cout << "NUM_EC:\t" << EC_in.size() << endl;
  cout << "NUM_LEAF_AC:\t" << AC_leaf.size() << endl;
  cout << "NUM_FC:\t" << num_FC << endl;

  int dim = ftr_name.size();
  double *min_val = new double[dim];
  double *max_val = new double[dim];
  for(int j=0;j<dim;j++){
    min_val[j] = ftr_val[0][j];
    max_val[j] = ftr_val[0][j];
    for(int k=1;k<Train_without_H.size();k++)
    if(ftr_val[k][j] < min_val[j])
      min_val[j] = ftr_val[k][j];
    else if(ftr_val[k][j] > max_val[j])
      max_val[j] = ftr_val[k][j];
  }
  cout << "N_RANGE:\t" << min_val[N_INDEX] << "\t" << max_val[N_INDEX] << endl;
  
  /***** postprocess *****/
  cout << Par.fringe_file << " is generated." << endl;

  // write data to csv file
  output_csv(Par.csv_file, Train_with_H, ftr_name, ftr_val);
  cout << Par.csv_file << " is generated." << endl;

  // normalization and output 
  normalize(ftr_val, min_val, max_val);
  output_csv(Par.norm_csv_file, Train_with_H, ftr_name, ftr_val);
  cout << Par.norm_csv_file << " is generated." << endl;

  /***** processing the test set *****/
  if(Par.test_use){
    cout << "TEST_NUM_VECTORS:\t" << Test_without_H.size() << endl;
    // write data to csv file
    output_csv(Par.test_csv, Test_with_H, ftr_name, test_val);
    cout << Par.test_csv << " is generated." << endl;
    
    // normalization and output 
    normalize(test_val, min_val, max_val);
    output_csv(Par.test_norm_csv, Test_with_H, ftr_name, test_val);
    cout << Par.test_norm_csv << " is generated." << endl;
  }
  
  delete[] min_val;
  delete[] max_val;

  return 0;
}
