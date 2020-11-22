// The main header file of the code.
// This files contains functions used to generate fringe trees, the sets W1^(h), W2^(h), W3^(h) and W4^(h),
// finding feasible pairs and outputting in SDF format.

#ifndef INCLUDES_FRINGE_TREE_HPP_
#define INCLUDES_FRINGE_TREE_HPP_

#include <iostream>
#include <cstdlib>
#include <vector>
#include <iomanip>

#include "data_structures.hpp"
#include "tools.hpp"

using namespace std;

double _start_time;
size_t time_limit = 3600;    // time limit used for one step of merging

bool _stop = false;

bool _stop_gen = false;
size_t num_limit = 100;   //  limit of the number of graph to output
unsigned long long UB_limit = 10000000;

size_t total_num = 0;

unordered_set <string> T_W1_F1;
unordered_set <string> T_W2_F1;
unordered_set <string> T_W3_F1;

// A bound using operaton (4)
bool simple_bound(
	const input_info& IN_INFO,
	vector_3D <unsigned short>& CVE_in,
	vector_3D <unsigned short>& CVE_out
){
	size_t M = IN_INFO.num_kind_atoms;
	vector <size_t> nb(M + 1, 0);

	for (size_t i = 1; i <= M; ++i){
		for (size_t j = 1; j <= M; ++j){
        	for (size_t k = 1; k <= 3; ++k){
            	if (i == j){
            		nb[i] += CVE_out[i][j][k];
            	} else {
            		nb[i] += CVE_out[i][j][k];
            	}
            }    
        }
    }

    for (size_t i = 1; i <= M; ++i){
    	if (IN_INFO.rv.CVE_out[i][i][0] - 
    		CVE_out[i][i][0] > IN_INFO.nb[i] - nb[i]){
    		return false;
    	}
    }
    return true;
}

/**
 * A recursive function to generate all possible assignments
 * of resources to a BMT in lexicographically descending order
 */
void extendTree(
		const input_info& IN_INFO,
		FringeTree &T, 		// Fringe tree being constructed
		const size_t &next_index,	// Index in T.seq to assign
									// multiplicity and color
		const size_t &residual_valence, // Residual valence at v_2 of BMT
		vector <FringeTree>& Tau, // To store the collection of enumerated trees
		const size_t &dmax,
		const size_t root_degree,      // 1 when end fringe tree, 2 when internal fringe tree
		const bool bound_check = false
){
	if (next_index >= dmax + 1) return; // We've assigned colors to each vertex in BMT
	// else:
	size_t parent_index = 1;
	MultCol &last_assignment = T.seq[next_index - 1]; // no need to make a copy
	Mult last_mult = last_assignment.first;
	Color last_color = last_assignment.second;
	Color parent_color = T.seq[parent_index].second;
	// Number of different atoms, for iteration
	Color max_color = IN_INFO.atoms.size()-1;
	if (next_index == 2){
		last_mult = 3;
		last_color = max_color;
	}

	// In order to get left-heavy trees, for the next assignment
	// we only consider MultColor >= last assignment
	for (Mult next_mult = last_mult; next_mult >= 0; --next_mult) {
		if (next_mult <= residual_valence) {
			for (Color color_it = max_color; color_it >= 0; --color_it) {
				Color next_color = color_it;
				if (next_mult == last_mult && next_color > last_color) continue;
				if (next_mult > 0 && next_color == 0) continue;
				if (next_mult == 0 && next_color > 0) continue;
				// Check if the candidate atom has valence
				// at least the next multiplicity

				if (next_mult <= IN_INFO.atoms[next_color].VALENCE) {
					// Check if there are resources
					if (IN_INFO.rv.CVE_out[parent_color][next_color][next_mult] > T.rv.CVE_out[parent_color][next_color][next_mult]
						&&
						IN_INFO.rv.CVE_out[next_color][next_color][0] > T.rv.CVE_out[next_color][next_color][0]) {
						// All conditions are satisfied
						MultCol
							next_assignment = MultCol(next_mult, next_color);
						T.seq[next_index] = next_assignment;
						++T.rv.CVE_out[parent_color][next_color][next_mult];
						++T.rv.CVE_out[next_color][parent_color][next_mult];
						++T.rv.CVE_out[next_color][next_color][0];
						if (next_color > 0) ++T.num_verts;
						// store the newly generated fringe tree
						
						// recurse

						if (next_index <= dmax) {
							extendTree( IN_INFO, T, next_index+1,
										residual_valence - next_mult,
										Tau, dmax, root_degree );
						} // if
						T.calc_bc(IN_INFO, dmax, root_degree, 2);
						if (dc_compare_out(IN_INFO, T.rv.bc_num_out)
						 // && simple_dc_bound(IN_INFO, T)
							&& ((!bound_check) || (bound_check && simple_bound(IN_INFO, T.rv.CVE_in, T.rv.CVE_out)))
						){
							Tau.push_back(T);	
						}

						// restore variables after recursion
						MultCol empty_assignment(0, 0);
						T.seq[next_index] = empty_assignment;
						--T.rv.CVE_out[parent_color][next_color][next_mult];
						--T.rv.CVE_out[next_color][parent_color][next_mult];
						--T.rv.CVE_out[next_color][next_color][0];
						if (next_color > 0) --T.num_verts;
					}
				} // if
			} // for
		} // if
	} // for loop over next_mult
} // function extendTree

/**
 * An initializer function that assigns color a to the root of
 * a BMT, then iterates
 */
void
rootedFringeTree(
		const input_info& IN_INFO,
		Color a, 					// Color for the root of the fringe tree
		vector<FringeTree> &Tau,		// store all generated trees
		const size_t &dmax,
		const size_t root_degree,
		const bool bound_check = false
		){
	size_t M = IN_INFO.num_kind_atoms;
	FringeTree T(M, dmax);
	T.seq[0] = MultCol(0, a);
	if (a != 0) ++T.num_verts;
	T.rv.CVE_in[a][a][0] = 1;
	if (IN_INFO.rv.CVE_in[a][a][0] == 0) return;
	// Tau.push_back(T);

	Color max_color = IN_INFO.atoms.size() - 1;
	Mult max_val;
	if (IN_INFO.atoms[a].VALENCE < 3 + root_degree){
		max_val = IN_INFO.atoms[a].VALENCE - root_degree;
	} else {
		max_val = 3;
	}
	for (Mult m = max_val; m >= 1; --m) {
		for (Color b = max_color; b >= 1; --b) {
			// Check resources and valence of b
			if (IN_INFO.rv.CVE_out[b][b][0] > T.rv.CVE_out[b][b][0]
				&&
				IN_INFO.rv.CVE_out[a][b][m] > 0
				&&
				m <= IN_INFO.atoms[b].VALENCE) {
				MultCol first_assignment(m, b);
				size_t residual_valence = IN_INFO.atoms[b].VALENCE - m;
				size_t next_index = 2;
				++T.rv.CVE_out[b][b][0];
				++T.rv.CVE_out[a][b][m];
				++T.rv.CVE_out[b][a][m];
				if (b > 0) ++T.num_verts;
				T.seq[1] = MultCol(m, b);
				T.root_val = m;
				if (b != 0)
					extendTree( IN_INFO, T, next_index,
								residual_valence, Tau, dmax, root_degree, bound_check);
				Tau.push_back(T);
				T.root_val = 0;
				--T.rv.CVE_out[b][b][0];
				--T.rv.CVE_out[a][b][m];
				--T.rv.CVE_out[b][a][m];
				if (b > 0) --T.num_verts;
			} // if there are available resources
		} // for color b
	} // for multiplicity m
	T.seq[1] = MultCol(0, 0);
	// T.rv.CVE[a][a][0] = 0;
	Tau.push_back(T);
} // function rootedFringeTree

/**
* The above is the old code to get output of step 1 of Fringe-tree
* This part should be modified by new trie and degree based bounds in phase 1
*/

//Generating Fringe-trees for a given color a: Step 2 //
vector <FringeTree> GenAllFringeTrees(
	const input_info& IN_INFO, //information from input file
	const size_t &index_a, // index of col a is given in Lambda; don't use 0 index
	vector<FringeTree> &AllFringeTrees, // all fringe trees of size at most 5
	const size_t &dmax, // max degree
	const size_t root_degree = 1,
	const bool height_check = true,
	const bool bound_check = true
	)
{
	size_t m = AllFringeTrees.size();
	vector<FringeTree> AllReqFringeTrees;// output of step 2
	
	for (size_t i = 0; i < m; ++i){
		AllFringeTrees[i].calc_bc(IN_INFO, dmax, root_degree, 1);
		AllFringeTrees[i].calc_deg(dmax, root_degree);
	}

	for (size_t i = 0; i < m; ++i) {
		size_t num_vert_Ti = AllFringeTrees[i].num_verts;
		resource_vector& ri = AllFringeTrees[i].rv;

		if (((height_check && num_vert_Ti >= 3 &&   // added by Zhu: the height of the firnge tree should be 2, which means
								   // that at least one of the used trees has at least 3 vertices
			num_vert_Ti <= 4) || (!height_check))
			&&
			TestLessEqAB(ri, IN_INFO.rv) &&
			AllFringeTrees[i].seq[1].first + root_degree <= IN_INFO.atoms[index_a].VALENCE) {
			FringeTree Tip(AllFringeTrees[i]);
			if (Tip.num_verts > 1) {
				Tip.root_degree = root_degree + 1;
			} else {
				Tip.root_degree = root_degree;
			}
			if ((!bound_check) || (bound_check && simple_bound(IN_INFO, ri.CVE_in, ri.CVE_out))){
				AllReqFringeTrees.push_back(Tip);
			}
		}
		for (size_t j = i; j < m; ++j) {
			size_t num_vert_Tj = AllFringeTrees[j].num_verts;

			//getting sum of resourses of Ti and Tj
			vector_3D <unsigned short> w_ijPlus_in = SUMinus3D(AllFringeTrees[i].rv.CVE_in,
				AllFringeTrees[j].rv.CVE_in, true);
			vector_3D <unsigned short> w_ijPlus_out = SUMinus3D(AllFringeTrees[i].rv.CVE_out,
				AllFringeTrees[j].rv.CVE_out, true);
			--w_ijPlus_in[index_a][index_a][0];
			// getting multiplicity between root and its children
			size_t k_i = AllFringeTrees[i].seq[1].first;
			size_t k_j = AllFringeTrees[j].seq[1].first;

			if ((k_i + k_j + root_degree <= IN_INFO.atoms[index_a].VALENCE) &&
				TestLessEqAB(w_ijPlus_in, IN_INFO.rv.CVE_in) &&
				TestLessEqAB(w_ijPlus_out, IN_INFO.rv.CVE_out) 
			){
				for (size_t h = j; h < m; ++h) {
					int num_vert_Th = AllFringeTrees[h].num_verts;
					//getting sum of resourses of Ti, Tj, T_h
					vector_3D <unsigned short> w_ijhPlus_in = SUMinus3D(w_ijPlus_in,
						AllFringeTrees[h].rv.CVE_in, true);
					vector_3D <unsigned short> w_ijhPlus_out = SUMinus3D(w_ijPlus_out,
						AllFringeTrees[h].rv.CVE_out, true);
					--w_ijhPlus_in[index_a][index_a][0];

					size_t k_h = AllFringeTrees[h].seq[1].first;

					if ((k_i + k_j + k_h + root_degree <= IN_INFO.atoms[index_a].VALENCE) &&
						TestLessEqAB(w_ijhPlus_in, IN_INFO.rv.CVE_in) &&
						TestLessEqAB(w_ijhPlus_out, IN_INFO.rv.CVE_out) &&
						((height_check &&
						(num_vert_Ti + num_vert_Tj + num_vert_Th - 2 <= 8) &&
						(num_vert_Ti >= 3 || num_vert_Tj >= 3 || num_vert_Th >= 3)) || (!height_check)) && // added by Zhu: height
						(num_vert_Ti > 1 && num_vert_Tj > 1 && num_vert_Th > 1)) { // added by Zhu: none of the trees can not be a tree with only one vertex (only root)
							
						vector <MultCol> TiTjTh;
						// constructing fringe tree sequence with Ti, Tj and Tk,
						//excluding "a" from Tj and Th
						for (const auto &s : AllFringeTrees[i].seq) {
							TiTjTh.push_back(s);
						}
						vector <size_t> P = { j, h };
						for (const auto &p : P) {
							for ( size_t q = 1;
									q < AllFringeTrees[p].seq.size();
									++q ) {
								//Note:AllFringeTrees[x].seq.size() is fixed to 5
								TiTjTh.push_back(AllFringeTrees[p].seq[q]);
							}
						}
						FringeTree Tijkp;
						Tijkp.seq = TiTjTh;
						Tijkp.rv.CVE_in = w_ijhPlus_in;
						Tijkp.rv.CVE_out = w_ijhPlus_out;
						Tijkp.num_verts =
								num_vert_Ti + num_vert_Tj + num_vert_Th - 2;
						Tijkp.calc_bc(IN_INFO, dmax, root_degree, 1);
						Tijkp.calc_deg(dmax, root_degree);
						Tijkp.root_degree = root_degree + 3;
						Tijkp.root_val = k_i + k_j + k_h;
						if (
							((!bound_check) || (bound_check && simple_bound(IN_INFO, w_ijhPlus_in, w_ijhPlus_out))) && 
							TestLessEqAB(Tijkp.rv, IN_INFO.rv)){
							AllReqFringeTrees.push_back(Tijkp);
						}
					}
				}
				
				if (((height_check && num_vert_Ti + num_vert_Tj - 1 <= 6 &&
					(num_vert_Ti >= 3 || num_vert_Tj >= 3)) || (!height_check))
					&& // added by Zhu: height
					(num_vert_Ti > 1 && num_vert_Tj > 1)) {  // added by Zhu: both of the trees can not be a tree with only one vertex (only root), or will be duplicated
					vector <MultCol> Tij;
					for (const auto &s : AllFringeTrees[i].seq) {
						Tij.push_back(s);
					}
					for (size_t jp = 1; jp < AllFringeTrees[j].seq.size(); ++jp) {
						Tij.push_back(AllFringeTrees[j].seq[jp]);
					}
					FringeTree Tijp;
					Tijp.seq = Tij;
					Tijp.rv.CVE_in = w_ijPlus_in;
					Tijp.rv.CVE_out = w_ijPlus_out;
					Tijp.num_verts = num_vert_Ti + num_vert_Tj - 1;

					Tijp.calc_bc(IN_INFO, dmax, root_degree, 1);
					Tijp.calc_deg(dmax, root_degree);

					Tijp.root_degree = root_degree + 2;
					Tijp.root_val = k_i + k_j;
					if (
						((!bound_check) || (bound_check && simple_bound(IN_INFO, w_ijPlus_in, w_ijPlus_out))) &&  
						TestLessEqAB(Tijp.rv, IN_INFO.rv)){
						AllReqFringeTrees.push_back(Tijp);
					}
				}
			}//if ((k_i + k_j <= IN_INFO.atoms[index_a].VALENCE - root_degree) &&
		}//for j
	}//for i
	return AllReqFringeTrees;
}

// A function to generate W1^(0)
void GenW1(const input_info& IN_INFO,
	vector <unordered_set <root_status>>& all_rs_W1,
	unordered_set <root_status>& all_rs_W1_F1,
	unordered_map <pair <unsigned short, size_t>, vector <size_t>>& map_FT_W1,
	vector <unordered_map <root_status, map <size_t, pair <map_rv, size_t>>>>& map_W1,
	unordered_map <root_status, unordered_set <size_t>>& map_W1_F1,
	vector_2D <FringeTree>& All_FT_W1,
	vector_2D <PathTrie>& all_trie_node
	// root status : the information of the root (color, deg)
){
	size_t M = IN_INFO.num_kind_atoms;
	map_W1[0].clear();
	
	for (unsigned short a = 1; a <= M; ++a){
		vector <FringeTree> Tau = {};

		rootedFringeTree(IN_INFO, a, Tau, IN_INFO.dmax, 1, false);
		All_FT_W1[a] = GenAllFringeTrees(IN_INFO, a, Tau, IN_INFO.dmax, 1, true, false);

		for (size_t i = 0; i < All_FT_W1[a].size(); ++i){
			auto& T = All_FT_W1[a][i];
			root_status rs(a, T.root_val, T.root_degree);
			
			resource_vector_1D rv1D = rv2rv1D(T.rv);
			size_t ind = find_trie_index(rv1D, all_trie_node[0]);

			map_FT_W1[make_pair(a, ind)].push_back(i);

			if (map_W1[0][rs].find(ind) == map_W1[0][rs].end()){
				map_W1[0][rs][ind] = make_pair(map_rv(a, ind), 1);
			} else {
				++map_W1[0][rs][ind].second;
			}

			all_rs_W1[0].emplace(a, T.root_val, T.root_degree);

			// check if T in F1 or not
			if (T_W1_F1.find(get_fringe_tree_string(T)) != T_W1_F1.end()){
				all_rs_W1_F1.emplace(a, T.root_val, T.root_degree);
				map_W1_F1[rs].insert(ind);
			}
		}// for i
	}// for a
}

// A function to generate W2^(0)
void GenW2_0(const input_info& IN_INFO,
	vector <unordered_set <half_path_status>>& all_rs_W2,
	unordered_set <half_path_status>& all_hps_F1,
	unordered_map <pair <unsigned short, size_t>, vector <size_t>>& map_FT_W2,
	vector <unordered_map <half_path_status, map <size_t, pair <map_rv, size_t>>>>& map_W2,
	unordered_map <half_path_status, unordered_set <size_t>>& map_W2_F1,
    unordered_map <half_path_status, vector <size_t>>& map_W2_ind,
	vector_2D <FringeTree>& All_FT_W2,
	vector_2D <PathTrie>& all_trie_node
	// half path status : (len, left root, right root)
){
	size_t M = IN_INFO.num_kind_atoms;
	map_W2[0].clear();

	for (unsigned short a = 1; a <= M; ++a){
		vector <FringeTree> Tau = {};
		rootedFringeTree(IN_INFO, a, Tau, IN_INFO.dmax, 2, true);
		All_FT_W2[a] = GenAllFringeTrees(IN_INFO, a, Tau, IN_INFO.dmax, 2, false, true);

		for (size_t i = 0; i < All_FT_W2[a].size(); ++i){
			auto& T = All_FT_W2[a][i];
			root_status rs(a, T.root_val, T.root_degree);
			half_path_status hps(rs, rs);
			all_rs_W2[0].emplace(rs, rs);

			resource_vector_1D rv1D = rv2rv1D(T.rv);
			size_t ind = find_trie_index(rv1D, all_trie_node[0]);

			map_FT_W2[make_pair(a, ind)].push_back(i);

			if (map_W2[0][hps].find(ind) == map_W2[0][hps].end()){
				map_W2[0][hps][ind] = make_pair(map_rv(a, ind), 1);
			} else {
				++map_W2[0][hps][ind].second;
			}

			if (T_W2_F1.find(get_fringe_tree_string(T)) != T_W2_F1.end()){
				all_hps_F1.emplace(rs, rs);
				map_W2_F1[hps].insert(ind);
			}
		}//for i
	}// for a

	for (auto& hps : all_rs_W2[0]){
		map_W2_ind[hps].clear();
		for (auto& tmp : map_W2[0].at(hps)){
			map_W2_ind[hps].push_back(tmp.first);
		}
	}
}

// A function to generate W3^(0)
void GenW3_0(const input_info& IN_INFO,
	vector <unordered_set <root_status>>& all_rs_W3,
	unordered_set <root_status>& all_rs_W3_F1,
	unordered_map <pair <unsigned short, size_t>, vector <size_t>>& map_FT_W3,
	vector <unordered_map <root_status, map <size_t, pair <map_rv, size_t>>>>& map_W3,
    unordered_map <root_status, unordered_set <size_t>>& map_W3_F1,
	vector_2D <FringeTree>& All_FT_W3,
	vector_2D <PathTrie>& all_trie_node
){
	size_t M = IN_INFO.num_kind_atoms;
	map_W3[0].clear();

	for (unsigned short a = 1; a <= M; ++a){
		vector <FringeTree> Tau = {};
		rootedFringeTree(IN_INFO, a, Tau, IN_INFO.dmax, 3, false);
		All_FT_W3[a] = GenAllFringeTrees(IN_INFO, a, Tau, IN_INFO.dmax, 3, false, false);

		for (size_t i = 0; i < All_FT_W3[a].size(); ++i){
			auto& T = All_FT_W3[a][i];
			root_status rs(a, T.root_val, T.root_degree);
			all_rs_W3[0].emplace(rs);

			resource_vector_1D rv1D = rv2rv1D(T.rv);
			size_t ind = find_trie_index(rv1D, all_trie_node[0]);

			map_FT_W3[make_pair(a, ind)].push_back(i);

			if (map_W3[0][rs].find(ind) == map_W3[0][rs].end()){
				map_W3[0][rs][ind] = make_pair(map_rv(a, ind), 1);
			} else {
				++map_W3[0][rs][ind].second;
			}

			if (T_W3_F1.find(get_fringe_tree_string(T)) != T_W3_F1.end()){
				all_rs_W3_F1.emplace(a, T.root_val, T.root_degree);
				map_W3_F1[rs].insert(ind);
			}
		}
	}
}

// A function to generate W1^(h), h = 1,2,3,...
void merge_W1(const input_info& IN_INFO,
	vector <unordered_set <root_status>>& all_rs_W1,
	unordered_set <root_status>& all_rs_W1_F1,
	vector <unordered_set <half_path_status>>& all_rs_W2,   // a vector to store all  possible (len,  left root,  right root) 
	unordered_set <half_path_status>& all_hps_F1,
	vector <unordered_map <root_status, map <size_t, pair <map_rv, size_t>>>>& map_W1,
	unordered_map <root_status, unordered_set <size_t>>& map_W1_F1,
	vector <unordered_map <half_path_status, map <size_t, pair <map_rv, size_t>>>>& map_W2,
	unordered_map <half_path_status, unordered_set <size_t>>& map_W2_F1,
    unordered_map <half_path_status, vector <size_t>>& map_W2_ind,
	vector_2D <PathTrie>& all_trie_node
){
	size_t M = IN_INFO.num_kind_atoms;

	size_t len = IN_INFO.k1 - 1;

	for (size_t l = 1; l <= len; ++l){

		_start_time = globalTimeKeeper::tt.diff();
		_stop = false;

		size_t l1 = l - 1;
		size_t l2 = 0;

		map_W1[l].clear();

		size_t vector_size = 0;

		if (l == 1){
			for (auto& rhps : all_hps_F1){
				if(_stop) break;
				for (auto& rs : all_rs_W1_F1){
					if (_stop) break;
					for (size_t k = 1; k <= 3; ++k){
						if (_stop) break;
						// a possible a1,a2,d1,d2,k
						size_t a1 = rs.color;
						size_t a2 = rhps.lrs.color;
						size_t d1 = rs.deg;
						size_t d2 = rhps.lrs.deg;
						size_t v1 = rs.val;
						size_t v2 = rhps.lrs.val;

						size_t delta2 = 0;
						if (l2 == 0) delta2 = 1;

						if (
							v1 + k <= IN_INFO.atoms[a1].VALENCE 
							&& 
							v2 + k + delta2 <= IN_INFO.atoms[a2].VALENCE
						){
							
							auto& set_left = map_W1_F1.at(rs);
							auto& set_right = map_W2_F1.at(rhps);
							
							for (auto& _w2 : set_right){
								if (_stop) break;
								auto w2 = traverse(_w2, all_trie_node[l2], IN_INFO.rv1D_size);
								auto& _w2_second = map_W2[l2].at(rhps).at(_w2);
								for (auto& _w1 : set_left){
									if (globalTimeKeeper::tt.diff() - _start_time > time_limit){

										// time limit reached.
										// cout << "len = " << l << ". Time over, limit is " << time_limit << "s." << endl;
										_stop = true;
									}
									if (vector_size > UB_limit){
										_stop = true;
									}
									if (_stop) break;
									
									auto w1 = traverse(_w1, all_trie_node[l1], IN_INFO.rv1D_size);
									auto& _w1_second = map_W1[l1].at(rs).at(_w1);
									
									auto new_rv = SUMinus3D(w1, w2, true);

									++new_rv[rv1D_ind_CVE_in(a1, a2, k)];
									++new_rv[rv1D_ind_bc_in(IN_INFO, d1, d2, k)];

									if (
										TestLessEqAB(new_rv, IN_INFO.rv1D)    // w1 + w2 + 1_mu + 1_gamma <= x_star
									){
										root_status _rs(rhps.rrs);
										_rs.val += k;

										all_rs_W1[l].emplace(_rs);
										size_t ind = find_trie_index(new_rv, all_trie_node[l]);

										auto& seq1 = _w1_second.first;
										auto& seq2 = _w2_second.first;

										if (map_W1[l][_rs].find(ind) == map_W1[l][_rs].end()){
											++vector_size;
											map_W1[l][_rs][ind] = make_pair(map_rv(seq1, seq2, k, true), _w1_second.second * _w2_second.second);
										} else {
											map_W1[l][_rs][ind].second += _w1_second.second * _w2_second.second;
										}					
									} // if (TestLessEqAB)
								} // for _w1
							}// for _w2
						}
					} // for k
				}// for lhps
			}// for rhps


			for (auto& rhps : all_rs_W2[l2]){
				if(_stop) break;
				for (auto& rs : all_rs_W1_F1){
					if (_stop) break;
					for (size_t k = 1; k <= 3; ++k){
						if (_stop) break;
						// a possible a1,a2,d1,d2,k
						size_t a1 = rs.color;
						size_t a2 = rhps.lrs.color;
						size_t d1 = rs.deg;
						size_t d2 = rhps.lrs.deg;
						size_t v1 = rs.val;
						size_t v2 = rhps.lrs.val;

						size_t delta2 = 0;
						if (l2 == 0) delta2 = 1;

						if (
							v1 + k <= IN_INFO.atoms[a1].VALENCE 
							&& 
							v2 + k + delta2 <= IN_INFO.atoms[a2].VALENCE
						){
							
							auto& set_left = map_W1_F1.at(rs);
							auto& set_right = map_W2_ind.at(rhps);

							size_t index_size = map_W2_ind.at(rhps).size();
							size_t start_index = (index_size - 1) * (l - 1) / len;
							size_t ind = start_index;
							bool flag = true;
							
							while (ind != start_index || flag){
								flag = false;
								if (ind >= index_size){
									ind -= index_size;
								}
								if (_stop) break;
								auto& _w2 = map_W2_ind[rhps][ind];
								++ind;
								if (ind == index_size){
									ind = 0;
								}
							
								if (map_W2_F1[rhps].find(_w2) != map_W2_F1[rhps].end()) continue;
								auto w2 = traverse(_w2, all_trie_node[l2], IN_INFO.rv1D_size);
								auto& _w2_second = map_W2[l2].at(rhps).at(_w2);
								for (auto& _w1 : set_left){
									if (globalTimeKeeper::tt.diff() - _start_time > time_limit){

										// time limit reached.
										// cout << "len = " << l << ". Time over, limit is " << time_limit << "s." << endl;
										_stop = true;
									}
									if (vector_size > UB_limit){
										_stop = true;
									}
									if (_stop) break;
									
									auto w1 = traverse(_w1, all_trie_node[l1], IN_INFO.rv1D_size);
									auto& _w1_second = map_W1[l1].at(rs).at(_w1);
									
									auto new_rv = SUMinus3D(w1, w2, true);

									++new_rv[rv1D_ind_CVE_in(a1, a2, k)];
									++new_rv[rv1D_ind_bc_in(IN_INFO, d1, d2, k)];

									if (
										TestLessEqAB(new_rv, IN_INFO.rv1D)    // w1 + w2 + 1_mu + 1_gamma <= x_star
									){
										root_status _rs(rhps.rrs);
										_rs.val += k;

										all_rs_W1[l].emplace(_rs);
										size_t ind = find_trie_index(new_rv, all_trie_node[l]);

										auto& seq1 = _w1_second.first;
										auto& seq2 = _w2_second.first;

										if (map_W1[l][_rs].find(ind) == map_W1[l][_rs].end()){
											++vector_size;
											map_W1[l][_rs][ind] = make_pair(map_rv(seq1, seq2, k, true), _w1_second.second * _w2_second.second);
										} else {
											map_W1[l][_rs][ind].second += _w1_second.second * _w2_second.second;
										}					
									} // if (TestLessEqAB)
								} // for _w1
							}// for _w2
						}
					} // for k
				}// for lhps
			}// for rhps

			for (auto& rhps : all_hps_F1){
				if(_stop) break;
				for (auto& rs : all_rs_W1[l1]){
					if (_stop) break;
					for (size_t k = 1; k <= 3; ++k){
						if (_stop) break;
						// a possible a1,a2,d1,d2,k
						size_t a1 = rs.color;
						size_t a2 = rhps.lrs.color;
						size_t d1 = rs.deg;
						size_t d2 = rhps.lrs.deg;
						size_t v1 = rs.val;
						size_t v2 = rhps.lrs.val;

						size_t delta2 = 0;
						if (l2 == 0) delta2 = 1;

						if (
							v1 + k <= IN_INFO.atoms[a1].VALENCE 
							&& 
							v2 + k + delta2 <= IN_INFO.atoms[a2].VALENCE
						){
							
							auto& set_left = map_W1[l1].at(rs);
							auto& set_right = map_W2_F1.at(rhps);
							
							for (auto& _w2 : set_right){
								if (_stop) break;
								auto w2 = traverse(_w2, all_trie_node[l2], IN_INFO.rv1D_size);
								auto& _w2_second = map_W2[l2].at(rhps).at(_w2);
								for (auto& _w1 : set_left){
									if (globalTimeKeeper::tt.diff() - _start_time > time_limit){

										// time limit reached.
										// cout << "len = " << l << ". Time over, limit is " << time_limit << "s." << endl;
										_stop = true;
									}
									if (vector_size > UB_limit){
										_stop = true;
									}
									if (map_W1_F1[rs].find(_w1.first) != map_W1_F1[rs].end()) continue;
									if (_stop) break;
									
									auto w1 = traverse(_w1.first, all_trie_node[l1], IN_INFO.rv1D_size);
									
									auto new_rv = SUMinus3D(w1, w2, true);

									++new_rv[rv1D_ind_CVE_in(a1, a2, k)];
									++new_rv[rv1D_ind_bc_in(IN_INFO, d1, d2, k)];

									if (
										TestLessEqAB(new_rv, IN_INFO.rv1D)    // w1 + w2 + 1_mu + 1_gamma <= x_star
									){
										root_status _rs(rhps.rrs);
										_rs.val += k;

										all_rs_W1[l].emplace(_rs);
										size_t ind = find_trie_index(new_rv, all_trie_node[l]);

										auto& seq1 = _w1.second.first;
										auto& seq2 = _w2_second.first;

										if (map_W1[l][_rs].find(ind) == map_W1[l][_rs].end()){
											++vector_size;
											map_W1[l][_rs][ind] = make_pair(map_rv(seq1, seq2, k, true), _w1.second.second * _w2_second.second);
										} else {
											map_W1[l][_rs][ind].second += _w1.second.second * _w2_second.second;
										}					
									} // if (TestLessEqAB)
								} // for _w1
							}// for _w2
						}
					} // for k
				}// for lhps
			}// for rhps


			for (auto& rhps : all_rs_W2[l2]){
				if(_stop) break;
				for (auto& rs : all_rs_W1[l1]){
					if (_stop) break;
					for (size_t k = 1; k <= 3; ++k){
						if (_stop) break;
						// a possible a1,a2,d1,d2,k
						size_t a1 = rs.color;
						size_t a2 = rhps.lrs.color;
						size_t d1 = rs.deg;
						size_t d2 = rhps.lrs.deg;
						size_t v1 = rs.val;
						size_t v2 = rhps.lrs.val;

						size_t delta2 = 0;
						if (l2 == 0) delta2 = 1;

						if (
							v1 + k <= IN_INFO.atoms[a1].VALENCE 
							&& 
							v2 + k + delta2 <= IN_INFO.atoms[a2].VALENCE
						){
							
							auto& set_left = map_W1[l1].at(rs);
							auto& set_right = map_W2_ind.at(rhps);

							size_t index_size = map_W2_ind.at(rhps).size();
							size_t start_index = (index_size - 1) * (l - 1) / len;
							size_t ind = start_index;
							bool flag = true;
							
							while (ind != start_index || flag){
								flag = false;
								if (ind >= index_size){
									ind -= index_size;
								}
								if (_stop) break;
								auto& _w2 = map_W2_ind[rhps][ind];
								++ind;
								if (ind == index_size){
									ind = 0;
								}
							
								if (map_W2_F1[rhps].find(_w2) != map_W2_F1[rhps].end()) continue;
								auto w2 = traverse(_w2, all_trie_node[l2], IN_INFO.rv1D_size);
								auto& _w2_second = map_W2[l2].at(rhps).at(_w2);
								for (auto& _w1 : set_left){
									if (globalTimeKeeper::tt.diff() - _start_time > time_limit){

										// time limit reached.
										// cout << "len = " << l << ". Time over, limit is " << time_limit << "s." << endl;
										_stop = true;
									}
									if (vector_size > UB_limit){
										_stop = true;
									}
									if (map_W1_F1[rs].find(_w1.first) != map_W1_F1[rs].end()) continue;
									if (_stop) break;
									
									auto w1 = traverse(_w1.first, all_trie_node[l1], IN_INFO.rv1D_size);
									
									auto new_rv = SUMinus3D(w1, w2, true);

									++new_rv[rv1D_ind_CVE_in(a1, a2, k)];
									++new_rv[rv1D_ind_bc_in(IN_INFO, d1, d2, k)];

									if (
										TestLessEqAB(new_rv, IN_INFO.rv1D)    // w1 + w2 + 1_mu + 1_gamma <= x_star
									){
										root_status _rs(rhps.rrs);
										_rs.val += k;

										all_rs_W1[l].emplace(_rs);
										size_t ind = find_trie_index(new_rv, all_trie_node[l]);

										auto& seq1 = _w1.second.first;
										auto& seq2 = _w2_second.first;

										if (map_W1[l][_rs].find(ind) == map_W1[l][_rs].end()){
											++vector_size;
											map_W1[l][_rs][ind] = make_pair(map_rv(seq1, seq2, k, true), _w1.second.second * _w2_second.second);
										} else {
											map_W1[l][_rs][ind].second += _w1.second.second * _w2_second.second;
										}					
									} // if (TestLessEqAB)
								} // for _w1
							}// for _w2
						}
					} // for k
				}// for lhps
			}// for rhps
		} else {
			for (auto& rhps : all_hps_F1){
				if(_stop) break;
				for (auto& rs : all_rs_W1[l1]){
					if (_stop) break;
					for (size_t k = 1; k <= 3; ++k){
						if (_stop) break;
						// a possible a1,a2,d1,d2,k
						size_t a1 = rs.color;
						size_t a2 = rhps.lrs.color;
						size_t d1 = rs.deg;
						size_t d2 = rhps.lrs.deg;
						size_t v1 = rs.val;
						size_t v2 = rhps.lrs.val;

						size_t delta2 = 0;
						if (l2 == 0) delta2 = 1;

						if (
							v1 + k <= IN_INFO.atoms[a1].VALENCE 
							&& 
							v2 + k + delta2 <= IN_INFO.atoms[a2].VALENCE
						){
							
							auto& set_left = map_W1[l1].at(rs);
							auto& set_right = map_W2_F1.at(rhps);
							
							for (auto& _w2 : set_right){
								if (_stop) break;
								auto w2 = traverse(_w2, all_trie_node[l2], IN_INFO.rv1D_size);
								auto& _w2_second = map_W2[l2].at(rhps).at(_w2);
								for (auto& _w1 : set_left){
									if (globalTimeKeeper::tt.diff() - _start_time > time_limit){

										// time limit reached.
										// cout << "len = " << l << ". Time over, limit is " << time_limit << "s." << endl;
										_stop = true;
									}
									if (vector_size > UB_limit){
										_stop = true;
									}
									if (_stop) break;
									
									auto w1 = traverse(_w1.first, all_trie_node[l1], IN_INFO.rv1D_size);
									
									auto new_rv = SUMinus3D(w1, w2, true);

									++new_rv[rv1D_ind_CVE_in(a1, a2, k)];
									++new_rv[rv1D_ind_bc_in(IN_INFO, d1, d2, k)];

									if (
										TestLessEqAB(new_rv, IN_INFO.rv1D)    // w1 + w2 + 1_mu + 1_gamma <= x_star
									){
										root_status _rs(rhps.rrs);
										_rs.val += k;

										all_rs_W1[l].emplace(_rs);
										size_t ind = find_trie_index(new_rv, all_trie_node[l]);

										auto& seq1 = _w1.second.first;
										auto& seq2 = _w2_second.first;

										if (map_W1[l][_rs].find(ind) == map_W1[l][_rs].end()){
											++vector_size;
											map_W1[l][_rs][ind] = make_pair(map_rv(seq1, seq2, k, true), _w1.second.second * _w2_second.second);
										} else {
											map_W1[l][_rs][ind].second += _w1.second.second * _w2_second.second;
										}					
									} // if (TestLessEqAB)
								} // for _w1
							}// for _w2
						}
					} // for k
				}// for lhps
			}// for rhps


			for (auto& rhps : all_rs_W2[l2]){
				if(_stop) break;
				for (auto& rs : all_rs_W1[l1]){
					if (_stop) break;
					for (size_t k = 1; k <= 3; ++k){
						if (_stop) break;
						// a possible a1,a2,d1,d2,k
						size_t a1 = rs.color;
						size_t a2 = rhps.lrs.color;
						size_t d1 = rs.deg;
						size_t d2 = rhps.lrs.deg;
						size_t v1 = rs.val;
						size_t v2 = rhps.lrs.val;

						size_t delta2 = 0;
						if (l2 == 0) delta2 = 1;

						if (
							v1 + k <= IN_INFO.atoms[a1].VALENCE 
							&& 
							v2 + k + delta2 <= IN_INFO.atoms[a2].VALENCE
						){
							
							auto& set_left = map_W1[l1].at(rs);
							auto& set_right = map_W2_ind.at(rhps);

							size_t index_size = map_W2_ind.at(rhps).size();
							size_t start_index = (index_size - 1) * (l - 1) / len;
							size_t ind = start_index;
							bool flag = true;
							
							while (ind != start_index || flag){
								flag = false;
								if (ind >= index_size){
									ind -= index_size;
								}
								if (_stop) break;
								auto& _w2 = map_W2_ind[rhps][ind];
								++ind;
								if (ind == index_size){
									ind = 0;
								}
							
								if (map_W2_F1[rhps].find(_w2) != map_W2_F1[rhps].end()) continue;
								auto w2 = traverse(_w2, all_trie_node[l2], IN_INFO.rv1D_size);
								auto& _w2_second = map_W2[l2].at(rhps).at(_w2);
								for (auto& _w1 : set_left){
									if (globalTimeKeeper::tt.diff() - _start_time > time_limit){

										// time limit reached.
										// cout << "len = " << l << ". Time over, limit is " << time_limit << "s." << endl;
										_stop = true;
									}
									if (vector_size > UB_limit){
										_stop = true;
									}
									if (_stop) break;
									
									auto w1 = traverse(_w1.first, all_trie_node[l1], IN_INFO.rv1D_size);
									
									auto new_rv = SUMinus3D(w1, w2, true);

									++new_rv[rv1D_ind_CVE_in(a1, a2, k)];
									++new_rv[rv1D_ind_bc_in(IN_INFO, d1, d2, k)];

									if (
										TestLessEqAB(new_rv, IN_INFO.rv1D)    // w1 + w2 + 1_mu + 1_gamma <= x_star
									){
										root_status _rs(rhps.rrs);
										_rs.val += k;

										all_rs_W1[l].emplace(_rs);
										size_t ind = find_trie_index(new_rv, all_trie_node[l]);

										auto& seq1 = _w1.second.first;
										auto& seq2 = _w2_second.first;

										if (map_W1[l][_rs].find(ind) == map_W1[l][_rs].end()){
											++vector_size;
											map_W1[l][_rs][ind] = make_pair(map_rv(seq1, seq2, k, true), _w1.second.second * _w2_second.second);
										} else {
											map_W1[l][_rs][ind].second += _w1.second.second * _w2_second.second;
										}					
									} // if (TestLessEqAB)
								} // for _w1
							}// for _w2
						}
					} // for k
				}// for lhps
			}// for rhps
		}

		if (l1 < IN_INFO.k3){
			destroy(map_W1[l1]);
			if (l1 != 0) destroy(all_trie_node[l1]);
		}

	}// for l
}

// A function to generate W3^(h).
void merge_W3(const input_info& IN_INFO,
	vector <unordered_set <root_status>>& all_rs_W1,
	vector <unordered_set <root_status>>& all_rs_W3,  
	unordered_set <root_status>& all_rs_W3_F1,
	vector <unordered_map <root_status, map <size_t, pair <map_rv, size_t>>>>& map_W1,
	vector <unordered_map <root_status, map <size_t, pair <map_rv, size_t>>>>& map_W3,
	unordered_map <root_status, unordered_set <size_t>>& map_W3_F1,
	vector_2D <PathTrie>& all_trie_node
){
	size_t M = IN_INFO.num_kind_atoms;

	size_t l3 = IN_INFO.k3;
	size_t l2_max = (IN_INFO.diameter - 6) / 2;
	size_t l1_min = IN_INFO.diameter - 5 - l2_max;
	size_t l1_max = IN_INFO.diameter - 5 - l3;

	for (size_t l = l1_min; l <= l1_max; ++l){
	// size_t l = len;

		_start_time = globalTimeKeeper::tt.diff();
		_stop = false;

		size_t l1 = l - 1;
		size_t l2 = 0;

		map_W3[l].clear();

		size_t vector_size = 0;

		for (auto& rs3 : all_rs_W3_F1){
			if(_stop) break;
			for (auto& rs1 : all_rs_W1[l1]){
				if (_stop) break;
				for (size_t k = 1; k <= 3; ++k){
					if (_stop) break;
					// a possible a1,a2,d1,d2,k
					size_t a1 = rs1.color;
					size_t a2 = rs3.color;
					size_t d1 = rs1.deg;
					size_t d2 = rs3.deg;
					size_t v1 = rs1.val;
					size_t v2 = rs3.val;

					if (
						v1 + k <= IN_INFO.atoms[a1].VALENCE 
						&& 
						v2 + k + 2 <= IN_INFO.atoms[a2].VALENCE
					){
						
						auto& set_left = map_W1[l1].at(rs1);
						auto& set_right = map_W3_F1.at(rs3);

						for (auto& _w2 : set_right){
							if (_stop) break;
							auto w2 = traverse(_w2, all_trie_node[l2], IN_INFO.rv1D_size);
							auto& _w2_second = map_W3[l2].at(rs3).at(_w2);
							for (auto& _w1 : set_left){
								if (globalTimeKeeper::tt.diff() - _start_time > time_limit){

									// time limit reached.
									// cout << "len = " << l << ". Time over, limit is " << time_limit << "s." << endl;
									_stop = true;
								}
								if (vector_size > UB_limit){
									_stop = true;
								}
								if (_stop) break;
								
								auto w1 = traverse(_w1.first, all_trie_node[l1], IN_INFO.rv1D_size);
								
								auto new_rv = SUMinus3D(w1, w2, true);

								++new_rv[rv1D_ind_CVE_in(a1, a2, k)];
								++new_rv[rv1D_ind_bc_in(IN_INFO, d1, d2, k)];

								if (
									TestLessEqAB(new_rv, IN_INFO.rv1D)    // w1 + w2 + 1_mu + 1_gamma <= x_star
								){
									root_status _rs(rs3);
									_rs.val += k;

									all_rs_W3[l].emplace(_rs);
									size_t ind = find_trie_index(new_rv, all_trie_node[l]);

									auto& seq1 = _w1.second.first;
									auto& seq2 = _w2_second.first;

									if (map_W3[l][_rs].find(ind) == map_W3[l][_rs].end()){
										++vector_size;
										map_W3[l][_rs][ind] = make_pair(map_rv(seq1, seq2, k, true), _w1.second.second * _w2_second.second);
									} else {
										map_W3[l][_rs][ind].second += _w1.second.second * _w2_second.second;
									}					
								} // if (TestLessEqAB)
							} // for _w1
						}// for _w2
					}
				} // for k
			}// for lhps
		}// for rhps

		for (auto& rs3 : all_rs_W3[l2]){
			if(_stop) break;
			for (auto& rs1 : all_rs_W1[l1]){
				if (_stop) break;
				for (size_t k = 1; k <= 3; ++k){
					if (_stop) break;
					// a possible a1,a2,d1,d2,k
					size_t a1 = rs1.color;
					size_t a2 = rs3.color;
					size_t d1 = rs1.deg;
					size_t d2 = rs3.deg;
					size_t v1 = rs1.val;
					size_t v2 = rs3.val;

					if (
						v1 + k <= IN_INFO.atoms[a1].VALENCE 
						&& 
						v2 + k + 2 <= IN_INFO.atoms[a2].VALENCE
					){
						
						auto& set_left = map_W1[l1].at(rs1);
						auto& set_right = map_W3[l2].at(rs3);

						for (auto& _w2 : set_right){
							if (_stop) break;
							if (map_W3_F1[rs3].find(_w2.first) != map_W3_F1[rs3].end()) continue;
							auto w2 = traverse(_w2.first, all_trie_node[l2], IN_INFO.rv1D_size);
							for (auto& _w1 : set_left){
								if (globalTimeKeeper::tt.diff() - _start_time > time_limit){

									// time limit reached.
									// cout << "len = " << l << ". Time over, limit is " << time_limit << "s." << endl;
									_stop = true;
								}
								if (vector_size > UB_limit){
									_stop = true;
								}
								if (_stop) break;
								
								auto w1 = traverse(_w1.first, all_trie_node[l1], IN_INFO.rv1D_size);
								
								auto new_rv = SUMinus3D(w1, w2, true);

								++new_rv[rv1D_ind_CVE_in(a1, a2, k)];
								++new_rv[rv1D_ind_bc_in(IN_INFO, d1, d2, k)];

								if (
									TestLessEqAB(new_rv, IN_INFO.rv1D)    // w1 + w2 + 1_mu + 1_gamma <= x_star
								){
									root_status _rs(rs3);
									_rs.val += k;

									all_rs_W3[l].emplace(_rs);
									size_t ind = find_trie_index(new_rv, all_trie_node[l]);

									auto& seq1 = _w1.second.first;
									auto& seq2 = _w2.second.first;

									if (map_W3[l][_rs].find(ind) == map_W3[l][_rs].end()){
										++vector_size;
										map_W3[l][_rs][ind] = make_pair(map_rv(seq1, seq2, k, true), _w1.second.second * _w2.second.second);
									} else {
										map_W3[l][_rs][ind].second += _w1.second.second * _w2.second.second;
									}					
								} // if (TestLessEqAB)
							} // for _w1
						}// for _w2
					}
				} // for k
			}// for lhps
		}// for rhps
	}// for l
}

// A function to generate W4^(h)
void merge_W4(const input_info& IN_INFO,
	vector <unordered_set <root_status>>& all_rs_W1,
	vector <unordered_set <root_status>>& all_rs_W3, 
	vector <unordered_set <root_status>>& all_rs_W4,
	vector <unordered_map <root_status, map <size_t, pair <map_rv, size_t>>>>& map_W1,
	vector <unordered_map <root_status, map <size_t, pair <map_rv, size_t>>>>& map_W3,
	vector <unordered_map <root_status, map <size_t, pair <map_rv, size_t>>>>& map_W4,   
	vector_2D <PathTrie>& all_trie_node
){
	size_t M = IN_INFO.num_kind_atoms;

	size_t l2_min = IN_INFO.k3;
	size_t l2_max = (IN_INFO.diameter - 6) / 2;

	size_t l = IN_INFO.diameter - 4;
	size_t trie_l_ind = IN_INFO.k1 + 1;

	// cout << "all_trie_node[trie_l_ind].size() = " << all_trie_node[trie_l_ind].size() <<endl;

	for (size_t l2 = l2_max; l2 >= l2_min; --l2){

		_start_time = globalTimeKeeper::tt.diff();
		_stop = false;

		size_t l1 = l - 1 - l2;

		all_rs_W4[l1].clear();
		map_W4[l1].clear();    // means that the (l1+1)-th vertex is the +3-vertex.

		size_t vector_size = 0;

		for (auto& rs3 : all_rs_W3[l1]){
			if(_stop) break;
			for (auto& rs1 : all_rs_W1[l2]){
				if (_stop) break;
				for (size_t k = 1; k <= 3; ++k){
					if (_stop) break;
					// a possible a1,a2,d1,d2,k
					size_t a1 = rs3.color;
					size_t a2 = rs1.color;
					size_t d1 = rs3.deg;
					size_t d2 = rs1.deg;
					size_t v1 = rs3.val;
					size_t v2 = rs1.val;

					if (
						v1 + k + 1 <= IN_INFO.atoms[a1].VALENCE 
						&& 
						v2 + k <= IN_INFO.atoms[a2].VALENCE
					){
						
						auto& set_left = map_W3[l1].at(rs3);
						auto& set_right = map_W1[l2].at(rs1);

						for (auto& _w1 : set_left){
							if (_stop) break;

							auto w1 = traverse(_w1.first, all_trie_node[l1], IN_INFO.rv1D_size);
							for (auto& _w2 : set_right){
								if (globalTimeKeeper::tt.diff() - _start_time > time_limit){
									// time limit reached.
									// cout << "l1 = " << l1 << ". Time over, limit is " << time_limit << "s." << endl;
									_stop = true;
								}
								if (vector_size > UB_limit){
									_stop = true;
								}
								if (_stop) break;
								
								auto w2 = traverse(_w2.first, all_trie_node[l2], IN_INFO.rv1D_size);
								
								auto new_rv = SUMinus3D(w1, w2, true);
								++new_rv[rv1D_ind_CVE_in(a1, a2, k)];
								++new_rv[rv1D_ind_bc_in(IN_INFO, d1, d2, k)];

								if (
									TestLessEqAB(new_rv, IN_INFO.rv1D)    // w1 + w2 + 1_mu + 1_gamma <= x_star
								){
									root_status _rs(rs3);
									_rs.val += k;

									all_rs_W4[l1].emplace(_rs);
									size_t ind = find_trie_index(new_rv, all_trie_node[trie_l_ind]);

									auto& seq1 = _w1.second.first;
									auto& seq2 = _w2.second.first;

									if (map_W4[l1][_rs].find(ind) == map_W4[l1][_rs].end()){
										++vector_size;
										map_W4[l1][_rs][ind] = make_pair(map_rv(seq1, seq2, k, false), _w1.second.second * _w2.second.second);
									} else {
										map_W4[l1][_rs][ind].second += _w1.second.second * _w2.second.second;
									}					
								} // if (TestLessEqAB)
							} // for _w2
						}// for _w1
					}
				} // for k
			}// for rs1
		}// for rs3

		if (l2 != IN_INFO.k3){
			destroy(map_W1[l2]);
		}
		destroy(map_W3[l1]);

		if (l2 == 0) break;
	}// for l

	for (size_t i = 0; i <= IN_INFO.k1; ++i){
		if (i != IN_INFO.k3){
			destroy(all_trie_node[i]);
		}
	}
}

// A function used to generate SDF, especially adding a fringe tree to a vertex
void add_fringe_tree(const input_info& IN_INFO,
	FringeTree& T,
	size_t& graph_ind,
	vector_2D <size_t>& graph_adj,
	vector <size_t>& graph_col
){
	auto& prt = _prt[IN_INFO.dmax];

    vector <size_t> ind_map(T.seq.size());
    size_t root_ind = graph_ind;
    ind_map[0] = root_ind;
    graph_col[root_ind] = T.seq[0].second;

    for (size_t i = 1; i < T.seq.size(); ++i){
    	if (T.seq[i].first != 0){
    		++graph_ind;
    		graph_col[graph_ind] = T.seq[i].second;
    		ind_map[i] = graph_ind;
    		graph_adj[graph_ind][ind_map[prt[i]]] = T.seq[i].first;
    		graph_adj[ind_map[prt[i]]][graph_ind] = T.seq[i].first;
    	}
    }
}

void output_SDF(const input_info& IN_INFO,
	vector_2D <size_t>& graph_adj,
	vector <size_t>& graph_col,
	string outputfilename
){

	ofstream output(outputfilename, ios::app);
	size_t n = IN_INFO.n;
	output << total_num << "\n";
	output << "3 branches" << "\n";
	output << "3 branches" << "\n";
	output << std::setw(3) << n << std::setw(3) << n - 1 << "  0  0  0  0  0  0  0  0999 V2000 " << "\n";

	for (size_t i = 0; i < n; ++i){
		string atom_symbol = IN_INFO.atoms[graph_col[i]].NAME;
		output << "    0.0000    0.0000    0.0000" <<
				  std::setw(3) << atom_symbol << "  0  0  0  0  0  0  0  0  0  0  0  0" << "\n";
	}

	for (size_t i = 0; i < n; ++i){
		for (size_t j = i + 1; j < n; ++j){
			if (graph_adj[i][j] > 0){
				output << std::setw(3) << i + 1 << std::setw(3) << j + 1 << std::setw(3) << graph_adj[i][j] << "  0  0  0  0" << "\n";
			}
		}
	}
	output << "M  END" << "\n";
	output << "$$$$" << "\n";
	output.close();
}

void prepare_SDF(const input_info& IN_INFO,
	vector <int>& st,
	map_rv& seq,
	unordered_map <pair <unsigned short, size_t>, vector <size_t>>& map_FT_W1,
	unordered_map <pair <unsigned short, size_t>, vector <size_t>>& map_FT_W2,
	unordered_map <pair <unsigned short, size_t>, vector <size_t>>& map_FT_W3,
	vector_2D <FringeTree>& All_FT_W1,
	vector_2D <FringeTree>& All_FT_W2,
	vector_2D <FringeTree>& All_FT_W3,
	size_t l1,
	string outputfilename
){
	// A function used to prepare the graph for outputting into SDF format.

	vector_2D <size_t> graph_adj(IN_INFO.n, vector <size_t> (IN_INFO.n, 0));
	vector <size_t> graph_col(IN_INFO.n, 0);

	size_t n_used = 0;
	size_t last_n = 0;

	size_t k3_root;

	for (size_t h = 0; h < seq.w_ind.size(); ++h){
		if (h != 0){
			if (h == IN_INFO.diameter - 3){
				graph_adj[n_used][k3_root] = seq.mul[h - 1];
				graph_adj[k3_root][n_used] = seq.mul[h - 1];
			} else {
				graph_adj[n_used][last_n] = seq.mul[h - 1];
				graph_adj[last_n][n_used] = seq.mul[h - 1];
			}
			
		}

		last_n = n_used;

		if (h == l1){
			k3_root = n_used;
		}

		if (h == 0 || h == IN_INFO.diameter - 4 || h == seq.w_ind.size() - 1){
			auto& FT_ind = map_FT_W1.at(make_pair(seq.rs[h], seq.w_ind[h]))[st[h]];	

			add_fringe_tree(IN_INFO, All_FT_W1[seq.rs[h]][FT_ind], n_used, graph_adj, graph_col);
		} else if (h == l1){
			auto& FT_ind = map_FT_W3.at(make_pair(seq.rs[h], seq.w_ind[h]))[st[h]];	

			add_fringe_tree(IN_INFO, All_FT_W3[seq.rs[h]][FT_ind], n_used, graph_adj, graph_col);
		} else {
			auto& FT_ind = map_FT_W2.at(make_pair(seq.rs[h], seq.w_ind[h]))[st[h]];	

			add_fringe_tree(IN_INFO, All_FT_W2[seq.rs[h]][FT_ind], n_used, graph_adj, graph_col);
		}
		++n_used;
	}

	++total_num;

	if (total_num > num_limit){
		_stop_gen = true;
		return;
	}

	output_SDF(IN_INFO, graph_adj, graph_col, outputfilename);
}

void generate(const input_info& IN_INFO,
	map_rv& seq,
	unordered_map <pair <unsigned short, size_t>, vector <size_t>>& map_FT_W1,
	unordered_map <pair <unsigned short, size_t>, vector <size_t>>& map_FT_W2,
	unordered_map <pair <unsigned short, size_t>, vector <size_t>>& map_FT_W3,
	vector_2D <FringeTree>& All_FT_W1,
	vector_2D <FringeTree>& All_FT_W2,
	vector_2D <FringeTree>& All_FT_W3,
	size_t l1,
	string outputfilename
){
	// A function to generate the indices of the used fringe trees.
	vector <int> st(seq.w_ind.size() + 1, 0);
	prepare_SDF(IN_INFO, st, seq, map_FT_W1, map_FT_W2, map_FT_W3, All_FT_W1, All_FT_W2, All_FT_W3, l1, outputfilename);
}

// A function to find feasible pairs
void combine_check(const input_info& IN_INFO,
	vector <unordered_set <root_status>>& all_rs_W1,
	vector <unordered_set <root_status>>& all_rs_W4,
	vector <unordered_map <root_status, map <size_t, pair <map_rv, size_t>>>>& map_W1,
	vector <unordered_map <root_status, map <size_t, pair <map_rv, size_t>>>>& map_W4,
	vector_2D <PathTrie>& all_trie_node,
	vector_2D <map_rv>& all_seq
){
	size_t M = IN_INFO.num_kind_atoms;
	size_t count = 0;
	size_t count_num = 0;

	size_t l3 = IN_INFO.k3;
	size_t l2_min = IN_INFO.k3;
	size_t l2_max = (IN_INFO.diameter - 6) / 2;
	size_t l1_min = IN_INFO.diameter - 5 - l2_max;
	size_t l1_max = IN_INFO.diameter - 5 - l3;
	size_t l = IN_INFO.diameter - 4;
	size_t trie_l_ind = IN_INFO.k1 + 1;

	_start_time = globalTimeKeeper::tt.diff();
	_stop = false;

	for (size_t l1 = l1_min; l1 <= l1_max; ++l1){
		if (_stop) break;
		size_t l2 = l - 1 - l1;
		for (auto& rs4 : all_rs_W4[l1]){
			if (_stop) break;
			auto& set_w1 = map_W4[l1].at(rs4);
			for (auto& rs1 : all_rs_W1[l3]){
				auto& _set_w2 = map_W1[l3].at(rs1);
				if (_stop) break;
				for (size_t k = 1; k <= 3; ++k){
					if (_stop) break;
					size_t a1 = rs4.color;
					size_t d1 = rs4.deg;
					size_t v1 = rs4.val;

					size_t a2 = rs1.color;
					size_t d2 = rs1.deg;
					size_t v2 = rs1.val;

					if (v1 + k <= IN_INFO.atoms[a1].VALENCE && v2 + k <= IN_INFO.atoms[a2].VALENCE){

						set <resource_vector_1D> set_w2;
						for (auto& _rv : _set_w2){
							auto _w2 = traverse(_rv.first, all_trie_node[l3], IN_INFO.rv1D_size);
							set_w2.insert(_w2);
						}

						for (auto& _w1 : set_w1){
							if (globalTimeKeeper::tt.diff() - _start_time > time_limit){

								// time limit reached.
								// cout << "Time over, limit is " << time_limit << "s." << endl;
								_stop = true;
							}
							if (_stop) break;

							auto w1 = traverse(_w1.first, all_trie_node[trie_l_ind], IN_INFO.rv1D_size);
							++w1[rv1D_ind_CVE_in(a1, a2, k)];
							++w1[rv1D_ind_bc_in(IN_INFO, d1, d2, k)];

							if (
								TestLessEqAB(w1, IN_INFO.rv1D)
							){
								auto _w2 = SUMinus3D(IN_INFO.rv1D, w1, false);

								if (set_w2.find(_w2) != set_w2.end()){
									++count;
									size_t ind_w2 = find_trie_index(_w2, all_trie_node[l3]);
									auto& _w2 = map_W1[l3].at(rs1).at(ind_w2);
									// cout << "count = " << count << " _w1.first = " << _w1.first << endl;
									size_t tmp = _w1.second.second * _w2.second;

									// if (l2 == l3){
									// 	if (l1 == l2 + 1){
									// 		count_num += (tmp + 5) / 6;
									// 	} else {
									// 		count_num += (tmp + 1) / 2;
									// 	}
									// } else {
									// 	if (l1 == l2 + 1){
									// 		count_num += (tmp + 5) / 6;
									// 	} else {
									// 		count_num += tmp;
									// 	}
									// }
									count_num += (tmp + 5) / 6;

									auto& seq1 = _w1.second.first;
									auto& seq2 = _w2.first;

									if (count <= num_limit){
										all_seq[l1].emplace_back(seq1, seq2, k, false);
									} 
								}
							}

						}
					} // if (v1 + k <= IN_INFO.atoms[a1].VALENCE && v2 + k <= IN_INFO.atoms[a2].VALENCE)
				} // for k
			} // for rs1
		} // for rs4
		// cout << "l1 = " << l1 << " count = " << count << endl;
	} // for l1

	size_t _size = 0;
	for (auto& temp : all_seq){
		_size += temp.size();
	}

	cout << "Number of feasible pairs = " << count << endl;
	cout << "A lower bound on the number of graphs = " << count_num << endl;
	cout << "Number of generated graphs = " << _size << endl;

}

// A function to generate SDF format files
void Gen_Graph(const input_info& IN_INFO,
	vector_2D <map_rv>& all_seq,
	unordered_map <pair <unsigned short, size_t>, vector <size_t>>& map_FT_W1,
	unordered_map <pair <unsigned short, size_t>, vector <size_t>>& map_FT_W2,
	unordered_map <pair <unsigned short, size_t>, vector <size_t>>& map_FT_W3,
	vector_2D <FringeTree>& All_FT_W1,
	vector_2D <FringeTree>& All_FT_W2,
	vector_2D <FringeTree>& All_FT_W3,
	string outputfilename
){
	ofstream output(outputfilename, ios::out);
	output.close();

	size_t l3 = IN_INFO.k3;
	size_t l2_min = IN_INFO.k3;
	size_t l2_max = (IN_INFO.diameter - 6) / 2;
	size_t l1_min = IN_INFO.diameter - 5 - l2_max;
	size_t l1_max = IN_INFO.diameter - 5 - l3;

	for (size_t l1 = l1_min; l1 <= l1_max; ++l1){
		for (auto& seq : all_seq[l1]){
			generate(IN_INFO, seq, map_FT_W1, map_FT_W2, map_FT_W3, All_FT_W1, All_FT_W2, All_FT_W3, l1, outputfilename);
		}
	}
}


#endif /* INCLUDES_FRINGE_TREE_HPP_ */

