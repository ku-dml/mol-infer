/************************************************************
  fv_common.h
************************************************************/

/***** functions for preprocess *****/
Atom2Int init_mass();
AtomBond2Int init_convalent_radius();
vector < ChemGraph > read_sdf(const string& fname);
ChemGraph suppress_H(const ChemGraph& G);
ChemGraph compute_ht_in_ex(ChemGraph& G, bool with_H);


/***** functions for postprocess *****/
void output_csv(const string csv_file,
		const vector < ChemGraph >& GS,
		const vector < string >& ftr_name,
		const vector < vector < double > >& ftr_val);


/***** functions for computing feature values *****/
void get_num_vertices(const vector < ChemGraph >& GS,
		      vector < vector < double > >& ftr_val);
void get_rank(const vector < ChemGraph >& GS,
	      vector < vector < double > >& ftr_val);
void get_num_interior_vertices(const vector < ChemGraph >& GS,
			       vector < vector < double > >& ftr_val);
void get_avg_mass(Atom2Int& Mass,
		  const vector < ChemGraph >& GS,
		  vector < vector < double > >& ftr_val);
void get_interior_deg_G(const vector < ChemGraph >& GS,
			vector < vector < double > >& ftr_val);
void get_interior_deg_G_in(const vector < ChemGraph >& GS,
			   vector < vector < double > >& ftr_val);
void get_interior_bond(const vector < ChemGraph >& GS,
		       vector < vector < double > >& ftr_val);
void get_num_hydrogens(Atom2Int& Val,
		       const vector < ChemGraph >& GS,
		       vector < vector < double > >& ftr_val);
void get_num_hydrodegs(Atom2Int& Val,
		       const vector < ChemGraph >& GS,
		       vector < vector < double > >& ftr_val);
Atom2Int get_exist_elements(const vector < ChemGraph >& GS, int in_ex);
void get_frequency(Atom2Int &Lambda, 
		   const vector < ChemGraph >& GS, int in_ex,
		   vector < vector < double > >& ftr_val);
Atom2Int get_exist_LeafAC(const vector < ChemGraph >& GS);
void get_LeafAC(Atom2Int& AC_leaf,
		const vector < ChemGraph >& GS,
		vector < vector < double  > >& ftr_val);
string get_EC_string(Atom2Int &Mass,
		     const string atom_0, const int deg_0, const int val_0,
		     const string atom_1, const int deg_1, const int val_1,
		     const int bond);
Atom2Int get_exist_EC(Atom2Int& Mass, const vector < ChemGraph >& GS,
		      const int in_ex);
void get_EC(Atom2Int& Mass, Atom2Int& EC_in, const vector < ChemGraph >& GS,
	    const int in_ex, vector < vector < double > >& ftr_val);

Atom2Int get_exist_opt_R(const vector < ChemGraph >& GS);
void get_opt_R(const AtomBond2Int& CR, const Atom2Int M, 
	       const vector < ChemGraph >& GS,
	       vector < vector < double > >& ftr_val);

/***** functions for transformation *****/
void normalize( vector < vector < double > >& ftr_val,
		double *min_val, double *max_val );
