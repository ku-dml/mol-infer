/************************************************************
  Feature vector generator for 2L (two-layered) model 

  Ver: 4.1
  - bugs are fixed

  Ver: 4.0
  - You can deal with polymer data

  Ver: 3.0
  - the correspondence between index and fringe-tree structure is output

  Ver: 2.0
  - correspond to 2L-model_v4
  - distribution of hydrogen-degrees is appended

 ************************************************************/

/********** include header files **********/
#include "fv_def.h"
#include "fv_common.h"

// main routine of computing fc 
void fc_main_v02(string sdf_file,
		 string fringe_file,
		 int rho,
		 vector < string >& fc_CID,
		 vector < string >& fc_ftr_name,
		 vector < vector < double > >& fc_ftr_val);


/********** read arguments **********/
void read_arguments(int argc, char *argv[], Param& Par){
  if(argc<3 || argc%2==0){
    fprintf(stderr, "usage: %s (input.sdf)(output.csv)[options]\n\n", argv[0]);
    fprintf(stderr, "options:\n");
    fprintf(stderr, "\t-F:\tfile name for fringe-trees (none)\n");
    fprintf(stderr, "\t-N:\tnormalization (false)\n");
    fprintf(stderr, "\t-S:\tstandardization (false)\n");
    fprintf(stderr, "\n");
    exit(1);
  }
  
  Par.sdf_file = argv[1];
  Par.csv_file = argv[2];
  Par.fringe_file = "";
  for(size_t c=3;c<argc;c++){
    if(argv[c][0]!='-'){
      fprintf(stderr, "error: the option <%s> is illegal.\n", argv[c]);
      exit(1);
    }
    if(argv[c][1]=='F'){
      Par.fringe_file = argv[c+1];
      c++;
    }
    else if(argv[c][1]=='N'){
      fprintf(stderr, "warning: -N is under construction.\n");
    }
    else if(argv[c][1]=='S'){
      fprintf(stderr, "warning: -S is under construction.\n");
    }
    else{
      fprintf(stderr, "error: option <%s> is not available.\n", argv[c]);
      exit(1);
    }
  }
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

  /***** postprocess *****/
  // write data to csv file
  output_csv(Par.csv_file, GS, ftr_name, ftr_val);  
  
  return 0;
}
