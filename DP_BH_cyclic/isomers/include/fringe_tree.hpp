// The main header file of the code.
// This files contains functions used to generate fringe trees, the sets W1^(h),
// W2^(h), W3^(h) and W4^(h), finding feasible pairs and outputting in SDF
// format.

#ifndef INCLUDES_FRINGE_TREE_HPP_
#define INCLUDES_FRINGE_TREE_HPP_

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#include "data_structures.hpp"
#include "tools.hpp"

using namespace std;

// double _start_time;
size_t time_limit = 3600;  // time limit used for one step of merging

bool _stop = false;

size_t pair_limit = 0;  // limit of the number of feasible pairs per component
unsigned long long UB_limit = 10000000;

size_t fringe_tree_sum = 0;
size_t fringe_tree_vectors_sum = 0;
double fringe_tree_vectors_time = 0;
size_t end_subtree_vectors_sum = 0;
double end_subtree_vectors_time = 0;
size_t rooted_core_subtrees_sum = 0;
double rooted_core_subtrees_time = 0;
size_t bi_rooted_core_subtrees_sum = 0;
double bi_rooted_core_subtrees_time = 0;
size_t feasible_pairs_sum = 0;
double feasible_pairs_time = 0;

size_t DAG_size_v = 0;
size_t DAG_size_e = 0;
double DAG_construct_time = 0;
double DAG_traverse_path_time = 0;
double DAG_traverse_graph_time = 0;
double DAG_traverse_num_time = 0;

/**
 * A recursive function to generate all possible assignments
 * of resources to a BMT in lexicographically descending order
 */
void extendTree(const input_info& IN_INFO,  // resource
                FringeTree& T,              // the FT in build
                const size_t& next_index, const size_t& residual_valence,
                vector<FringeTree>& Tau,  // the set of all generated trees
                const size_t& root_degree) {
  size_t dmax = IN_INFO.dmax;
  if (next_index >= dmax + 1) return;
  size_t parent_index = 1;
  MultCol& last_assignment = T.seq[next_index - 1];
  Mult last_mult = last_assignment.first;
  Color last_color = last_assignment.second;
  Color parent_color = T.seq[parent_index].second;

  Color max_color = IN_INFO.atoms.size() - 1;
  if (next_index == 2) {
    last_mult = 3;
    last_color = max_color;
  }

  // In order to get left-heavy trees, for the next assignment
  // we only consider MultCol >= last assignment
  for (Mult next_mult = last_mult; next_mult > 0; --next_mult) {
    if (next_mult <= residual_valence) {
      for (Color color_it = max_color; color_it >= 0; --color_it) {
        Color next_color = color_it;
        if (next_mult == last_mult && next_color > last_color) continue;
        if (next_mult == 0 && next_color >= 0) continue;

        if (next_mult <= IN_INFO.atoms[next_color].VALENCE) {
          MultCol next_assignment = MultCol(next_mult, next_color);
          T.seq[next_index] = next_assignment;

          if (next_mult > 0) ++T.num_verts;
          if (next_index == 2) ++T.height;

          if (next_index <= dmax) {
            extendTree(IN_INFO, T, next_index + 1, residual_valence - next_mult,
                       Tau, root_degree);
          }

          Tau.push_back(T);
          if (next_index == 2) --T.height;

          MultCol empty_assignment(0, 0);
          T.seq[next_index] = empty_assignment;
          if (next_mult > 0) --T.num_verts;
        }
      }  // for color_it
    }    // if
  }      // for next_mult
}

/**
 * An initializer function that assigns color a to the root of
 * a BMT, then iterates
 */
void rootedFringeTree(const input_info& IN_INFO, Color a,
                      vector<FringeTree>& Tau, const size_t root_degree) {
  size_t M = IN_INFO.num_kind_atoms;
  size_t dmax = IN_INFO.dmax;
  FringeTree T(M, dmax);
  T.height = 0;
  T.seq[0] = MultCol(0, a);
  T.num_verts = 1;

  Color max_color = IN_INFO.atoms.size() - 1;
  Mult max_val;
  if (IN_INFO.atoms[a].VALENCE < 3 + root_degree) {
    max_val = IN_INFO.atoms[a].VALENCE - root_degree;
  } else {
    max_val = 3;
  }
  for (Mult m = max_val; m >= 1; --m) {
    for (Color b = max_color; b >= 0; --b) {
      if (m <= IN_INFO.atoms[b].VALENCE) {
        MultCol first_assignment(m, b);
        size_t residual_valence = IN_INFO.atoms[b].VALENCE - m;
        size_t next_index = 2;
        ++T.num_verts;
        T.seq[1] = first_assignment;
        T.root_val = m;
        T.height = 1;
        extendTree(IN_INFO, T, next_index, residual_valence, Tau, root_degree);
        Tau.push_back(T);
        T.root_val = 0;
        --T.num_verts;
      }
    }  // for b
  }    // for m
  T.seq[1] = MultCol(0, 0);
  T.height = 0;
  Tau.push_back(T);
}

/**
 * The above is the old code to get output of step 1 of Fringe-tree
 * This part should be modified by new trie and degree based bounds in phase 1
 */

// Generating Fringe-trees for a given color a: Step 2 //
vector<FringeTree> GenAllFringeTrees(const input_info& IN_INFO,
                                     const size_t& index_a,
                                     vector<FringeTree>& AllFringeTrees,
                                     const size_t root_degree,
                                     const bool height_check = true) {
  size_t m = AllFringeTrees.size();
  vector<FringeTree> AllReqFringeTrees;

  for (size_t i = 0; i < m; ++i) {
    AllFringeTrees[i].calc_rv(IN_INFO, root_degree);
  }
  for (size_t i = 0; i < m; ++i) {
    size_t num_vert_Ti = AllFringeTrees[i].num_verts;
    resource_vector& ri = AllFringeTrees[i].rv;

    // printseq(IN_INFO, AllFringeTrees[i].seq);
    // cout << "num_vert_Ti = " << num_vert_Ti << "
    // AllFringeTrees[i].seq[1].first + root_degree = " <<
    // AllFringeTrees[i].seq[1].first + root_degree << endl;

    if (num_vert_Ti <= 4 &&
        ((height_check && num_vert_Ti >= 3) || (!height_check)) &&
        TestLessEqAB(ri, IN_INFO.rv) &&
        AllFringeTrees[i].seq[1].first + root_degree <=
            IN_INFO.atoms[index_a].VALENCE) {
      FringeTree Tip(AllFringeTrees[i]);
      if (Tip.num_verts > 1) {
        Tip.root_degree = root_degree + 1;
      } else {
        Tip.root_degree = root_degree;
      }
      AllReqFringeTrees.push_back(Tip);
    }

    for (size_t j = i; j < m; ++j) {
      size_t num_vert_Tj = AllFringeTrees[j].num_verts;
      size_t k_i = AllFringeTrees[i].seq[1].first;
      size_t k_j = AllFringeTrees[j].seq[1].first;

      if (k_i + k_j + root_degree <= IN_INFO.atoms[index_a].VALENCE) {
        for (size_t h = j; h < m; ++h) {
          size_t num_vert_Th = AllFringeTrees[h].num_verts;
          size_t k_h = AllFringeTrees[h].seq[1].first;

          if ((k_i + k_j + k_h + root_degree <=
               IN_INFO.atoms[index_a].VALENCE) &&
              num_vert_Ti + num_vert_Tj + num_vert_Th - 2 <= 8 &&
              ((height_check &&
                (num_vert_Ti >= 3 || num_vert_Tj >= 3 || num_vert_Th >= 3)) ||
               (!height_check)) &&
              (num_vert_Ti > 1 && num_vert_Tj > 1 && num_vert_Th > 1)) {
            vector<MultCol> TiTjTh;
            for (const auto& s : AllFringeTrees[i].seq) {
              TiTjTh.push_back(s);
            }
            vector<size_t> P = {j, h};
            for (const auto& p : P) {
              for (size_t q = 1; q < AllFringeTrees[p].seq.size(); ++q) {
                TiTjTh.push_back(AllFringeTrees[p].seq[q]);
              }
            }
            FringeTree Tijkp;
            Tijkp.seq = TiTjTh;
            Tijkp.num_verts = num_vert_Ti + num_vert_Tj + num_vert_Th - 2;

            Tijkp.calc_rv(IN_INFO, root_degree);
            Tijkp.root_degree = root_degree + 3;
            Tijkp.root_val = k_i + k_j + k_h;
            Tijkp.height = AllFringeTrees[i].height;
            if (AllFringeTrees[j].height > Tijkp.height)
              Tijkp.height = AllFringeTrees[j].height;
            if (AllFringeTrees[h].height > Tijkp.height)
              Tijkp.height = AllFringeTrees[h].height;

            if (TestLessEqAB(Tijkp.rv, IN_INFO.rv)) {
              AllReqFringeTrees.push_back(Tijkp);
            }
          }
        }  // for h

        if (num_vert_Ti + num_vert_Tj - 1 <= 6 &&
            ((height_check && (num_vert_Ti >= 3 || num_vert_Tj >= 3)) ||
             (!height_check)) &&
            (num_vert_Ti > 1 && num_vert_Tj > 1)) {
          vector<MultCol> Tij;
          for (const auto& s : AllFringeTrees[i].seq) {
            Tij.push_back(s);
          }
          for (size_t jp = 1; jp < AllFringeTrees[j].seq.size(); ++jp) {
            Tij.push_back(AllFringeTrees[j].seq[jp]);
          }
          FringeTree Tijp;
          Tijp.seq = Tij;
          Tijp.num_verts = num_vert_Ti + num_vert_Tj - 1;

          Tijp.calc_rv(IN_INFO, root_degree);
          Tijp.root_degree = root_degree + 2;
          Tijp.root_val = k_i + k_j;
          Tijp.height = AllFringeTrees[i].height;
          if (AllFringeTrees[j].height > Tijp.height)
            Tijp.height = AllFringeTrees[j].height;

          if (TestLessEqAB(Tijp.rv, IN_INFO.rv)) {
            AllReqFringeTrees.push_back(Tijp);
          }
        }
      }  //  if k_i + k_j
    }    // for j
  }      // for i
  return AllReqFringeTrees;
}

// Gene
void GenW1(
    const input_info& IN_INFO, vector<unordered_set<root_status>>& all_rs_W1,
    unordered_map<pair<unsigned short, size_t>, vector<size_t>>& map_FT_W1,
    vector<unordered_map<root_status, map<size_t, bool>>>& map_W1,
    vector_2D<FringeTree>& All_FT_W1, vector_2D<PathTrie>& all_trie_node) {
  size_t M = IN_INFO.num_kind_atoms;
  map_W1[0].clear();

  for (unsigned short a = 0; a < M; ++a) {
    vector<FringeTree> Tau = {};

    rootedFringeTree(IN_INFO, a, Tau, 1);
    All_FT_W1[a] = GenAllFringeTrees(IN_INFO, a, Tau, 1, true);

    for (size_t i = 0; i < All_FT_W1[a].size(); ++i) {
      auto& T = All_FT_W1[a][i];
      T.rv.bl = 1;

      root_status rs(a, T.root_val,
                     T.root_degree);  // here T.root_degree means the real
                                      // degree, have added the [+2]
      all_rs_W1[0].emplace(a, T.root_val, T.root_degree);

      resource_vector_1D rv1D = rv2rv1D(T.rv, IN_INFO);
      size_t ind = find_trie_index(rv1D, all_trie_node[0]);
      // for (auto& atom : IN_INFO.atoms) {
      //   cout << atom.NAME << ": " << atom.LABEL << endl;
      // }

      if (_DEBUG) {
        cout << "W1 ind = " << ind << endl;
        print(T, IN_INFO);
      }

      map_FT_W1[make_pair(a, ind)].push_back(i);

      if (map_W1[0][rs].find(ind) == map_W1[0][rs].end()) {
        map_W1[0][rs][ind] = false;
      }
    }
    // if (_DEBUG) cout << "a = " << a <<" All_FT_W1[a].size() = " <<
    // All_FT_W1[a].size() << endl;

    if (_EXP) fringe_tree_sum += All_FT_W1[a].size();
  }

  size_t _size = 0;
  for (auto& _rs : all_rs_W1[0]) {
    _size += map_W1[0].at(_rs).size();
  }
  if (_EXP) fringe_tree_vectors_sum += _size;
}

void GenW2_0(
    const input_info& IN_INFO, unordered_set<half_path_status>& all_rs_W2,
    unordered_map<pair<unsigned short, size_t>, vector<size_t>>& map_FT_W2,
    unordered_map<half_path_status, map<size_t, bool>>& map_W2,
    unordered_map<half_path_status, vector<size_t>>& map_W2_ind,
    vector_2D<FringeTree>& All_FT_W2, vector_2D<PathTrie>& all_trie_node) {
  size_t M = IN_INFO.num_kind_atoms;
  map_W2.clear();

  for (unsigned short a = 0; a < M; ++a) {
    vector<FringeTree> Tau = {};

    rootedFringeTree(IN_INFO, a, Tau, 2);
    All_FT_W2[a] = GenAllFringeTrees(IN_INFO, a, Tau, 2, false);

    for (size_t i = 0; i < All_FT_W2[a].size(); ++i) {
      auto& T = All_FT_W2[a][i];
      root_status rs(a, T.root_val,
                     T.root_degree);  // here T.root_degree means the real
                                      // degree, have added the [+2]
      half_path_status hps(rs, rs);
      all_rs_W2.emplace(rs, rs);

      resource_vector_1D rv1D = rv2rv1D(T.rv, IN_INFO);
      size_t ind = find_trie_index(rv1D, all_trie_node[0]);

      if (_DEBUG) {
        cout << "W2_0 ind = " << ind << endl;
        print(T, IN_INFO);
      }

      map_FT_W2[make_pair(a, ind)].push_back(i);

      if (map_W2[hps].find(ind) == map_W2[hps].end()) {
        map_W2[hps][ind] = false;
      }
    }
    // if (_DEBUG) cout << "a = " << a <<" All_FT_W2[a].size() = " <<
    // All_FT_W2[a].size() << endl;

    if (_EXP) fringe_tree_sum += All_FT_W2[a].size();
  }

  for (auto& hps : all_rs_W2) {
    map_W2_ind[hps].clear();
    for (auto& tmp : map_W2.at(hps)) {
      map_W2_ind[hps].push_back(tmp.first);
    }
  }

  size_t _size = 0;
  for (auto& _rs : all_rs_W2) {
    _size += map_W2.at(_rs).size();
  }
  if (_EXP) fringe_tree_vectors_sum += _size;
}

void GenW2_0_core(
    const input_info& IN_INFO,
    unordered_set<root_status_height>& all_rs_W2_core,
    unordered_map<pair<unsigned short, size_t>, vector<size_t>>& map_FT_W2_core,
    unordered_map<root_status_height, map<size_t, bool>>& map_W2_core,
    vector_2D<FringeTree>& All_FT_W2_core, vector_2D<PathTrie>& all_trie_node,
    size_t delta_i = 2) {
  size_t M = IN_INFO.num_kind_atoms;
  map_W2_core.clear();

  for (unsigned short a = 0; a < M; ++a) {
    vector<FringeTree> Tau = {};

    rootedFringeTree(IN_INFO, a, Tau, delta_i);
    All_FT_W2_core[a] = GenAllFringeTrees(IN_INFO, a, Tau, delta_i, false);

    for (size_t i = 0; i < All_FT_W2_core[a].size(); ++i) {
      auto& T = All_FT_W2_core[a][i];
      root_status_height rsh(a, T.root_val, T.root_degree,
                             T.height);  // here T.root_degree means the real
                                         // degree, have added the [+2]
      all_rs_W2_core.emplace(rsh);
      // cout << "rsh.color = "  << rsh.color << " rsh.deg = " << rsh.deg << "
      // rsh.val = " <<rsh.val << " rsh.height = " << rsh.height << endl;

      resource_vector_1D rv1D = rv2rv1D(T.rv, IN_INFO);
      size_t ind = find_trie_index(rv1D, all_trie_node[0]);

      if (_DEBUG) {
        cout << "W2_0_core ind = " << ind << endl;
        print(T, IN_INFO);
      }

      map_FT_W2_core[make_pair(a, ind)].push_back(i);

      if (map_W2_core[rsh].find(ind) == map_W2_core[rsh].end()) {
        map_W2_core[rsh][ind] = false;
      }
    }
    // if (_DEBUG) cout << "a = " << a << " M = " << M << "
    // All_FT_W2_core[a].size() = " << All_FT_W2_core[a].size() << endl;

    if (_EXP) fringe_tree_sum += All_FT_W2_core[a].size();
  }

  size_t _size = 0;
  for (auto& _rs : all_rs_W2_core) {
    _size += map_W2_core.at(_rs).size();
  }
  if (_EXP) fringe_tree_vectors_sum += _size;
}

void GenW3_0(
    const input_info& IN_INFO, vector<unordered_set<root_status>>& all_rs_W3,
    unordered_map<pair<unsigned short, size_t>, vector<size_t>>& map_FT_W3,
    vector<unordered_map<root_status, map<size_t, bool>>>& map_W3,
    vector_2D<FringeTree>& All_FT_W3, vector_2D<PathTrie>& all_trie_node,
    size_t color, size_t delta_i = 3) {
  size_t M = IN_INFO.num_kind_atoms;
  map_W3[0].clear();

  _start_time = globalTimeKeeper::tt.diff();

  for (unsigned short a = 0; a < M; ++a) {
    if (color < M && a != color) continue;
    vector<FringeTree> Tau = {};

    rootedFringeTree(IN_INFO, a, Tau, delta_i);
    All_FT_W3[a] = GenAllFringeTrees(IN_INFO, a, Tau, delta_i, false);

    for (size_t i = 0; i < All_FT_W3[a].size(); ++i) {
      auto& T = All_FT_W3[a][i];
      root_status rs(a, T.root_val,
                     T.root_degree);  // here T.root_degree means the real
                                      // degree, have added the [+2]
      all_rs_W3[0].emplace(a, T.root_val, T.root_degree);

      resource_vector_1D rv1D = rv2rv1D(T.rv, IN_INFO);
      size_t ind = find_trie_index(rv1D, all_trie_node[0]);

      if (_DEBUG) {
        cout << "W3 ind = " << ind << endl;
        print(T, IN_INFO);
      }

      map_FT_W3[make_pair(a, ind)].push_back(i);

      if (map_W3[0][rs].find(ind) == map_W3[0][rs].end()) {
        map_W3[0][rs][ind] = false;
      }
    }
    // if (_DEBUG) cout << "a = " << a <<" All_FT_W3[a].size() = " <<
    // All_FT_W3[a].size() << endl;

    if (_EXP) fringe_tree_sum += All_FT_W3[a].size();
  }

  size_t _size = 0;
  for (auto& _rs : all_rs_W3[0]) {
    _size += map_W3[0].at(_rs).size();
  }
  if (_EXP) fringe_tree_vectors_sum += _size;
}

void merge_W1(const input_info& IN_INFO,
              vector<unordered_set<root_status>>& all_rs_W1,
              unordered_set<half_path_status>& all_rs_W2,
              vector<unordered_map<root_status, map<size_t, bool>>>& map_W1,
              unordered_map<half_path_status, map<size_t, bool>>& map_W2,
              unordered_map<half_path_status, vector<size_t>>& map_W2_ind,
              vector_2D<PathTrie>& all_trie_node,
              size_t len  // max core height - 1
) {
  size_t M = IN_INFO.num_kind_atoms;

  for (size_t l = 1; l <= len; ++l) {
    _start_time = globalTimeKeeper::tt.diff();
    _stop = false;

    size_t l1 = l - 1;
    size_t l2 = 0;

    map_W1[l].clear();

    size_t vector_size = 0;

    for (auto& rhps : all_rs_W2) {
      if (_stop) break;
      for (auto& rs : all_rs_W1[l1]) {
        if (_stop) break;
        for (size_t k = 1; k <= 3; ++k) {
          if (_stop) break;
          size_t a1 = rs.color;
          size_t a2 = rhps.lrs.color;
          size_t d1 = rs.deg;
          size_t d2 = rhps.lrs.deg;
          size_t m1 = rs.val;
          size_t m2 = rhps.lrs.val;

          atom_degree ad1(a1, d1);
          atom_degree ad2(a2, d2);

          if (m1 + k <= IN_INFO.atoms[a1].VALENCE &&
              m2 + k + 1 <= IN_INFO.atoms[a2].VALENCE) {
            auto& set_left = map_W1[l1].at(rs);
            auto& set_right = map_W2_ind.at(rhps);

            rv_edge e_tmp("in", ad2, ad1, k);  // direction: (a2,d2)->(a1,d1)
            size_t e_tmp_ind = IN_INFO.rv_edge_map.at(e_tmp);

            size_t index_size = map_W2_ind.at(rhps).size();
            size_t start_index = (index_size - 1) * (l - 1) / len;
            size_t ind = start_index;
            bool flag = true;

            while (ind != start_index || flag) {
              flag = false;
              if (ind >= index_size) {
                ind -= index_size;
              }
              if (_stop) break;
              auto& _w2 = map_W2_ind[rhps][ind];
              ++ind;
              if (ind == index_size) {
                ind = 0;
              }

              auto w2 = traverse(_w2, all_trie_node[l2], IN_INFO.rv1D_size);
              auto& _w2_second = map_W2.at(rhps).at(_w2);
              for (auto& _w1 : set_left) {
                if (globalTimeKeeper::tt.diff() - _start_time > time_limit) {
                  // time limit reached.
                  // cout << "len = " << l << ". Time over, limit is " <<
                  // time_limit << "s." << endl;
                  _stop = true;
                }
                if (UB_limit != 0 && vector_size > UB_limit) {
                  _stop = true;
                }
                if (_stop) break;

                auto w1 =
                    traverse(_w1.first, all_trie_node[l1], IN_INFO.rv1D_size);
                auto new_rv = SUMinus3D(w1, w2, true);

                ++new_rv[e_tmp_ind];

                if (TestLessEqAB(
                        new_rv,
                        IN_INFO.rv1D)  // w1 + w2 + 1_mu + 1_gamma <= x_star
                ) {
                  root_status _rs(rhps.rrs);
                  _rs.val += k;

                  all_rs_W1[l].emplace(_rs);
                  size_t ind = find_trie_index(new_rv, all_trie_node[l]);

                  // auto& seq1_vector = _w1.second.first;
                  // auto& seq2_vector = _w2_second.first;

                  if (map_W1[l][_rs].find(ind) == map_W1[l][_rs].end()) {
                    ++vector_size;
                    // vector <map_rv> tmp;
                    // tmp.clear();
                    map_W1[l][_rs][ind] = false;
                  }

                  // auto& _seq = map_W1[l][_rs][ind].first;
                  // auto& _num = map_W1[l][_rs][ind].second;

                  // for (auto& seq1 : seq1_vector){
                  // 	for (auto& seq2 : seq2_vector){
                  // 		if (seq_limit == 0 || _seq.size() < seq_limit){
                  // 			_seq.emplace_back(seq1, seq2, k, true);
                  // 		}
                  // 	} // for seq2
                  // }// for seq1
                  // _num += _w1.second.second * _w2_second.second;
                }  // if (TestLessEqAB)
              }    // for _w1
            }      // while
          }        // if
        }          // for k
      }            // for rs
    }              // for rhps

    if (_DEBUG)
      cout << "end subtrees l = " << l << " len = " << len
           << " vector_size = " << vector_size
           << " time = " << globalTimeKeeper::tt.diff() - _start_time << endl;
    end_subtree_vectors_sum += vector_size;
  }  // for l
}

void merge_W3(const input_info& IN_INFO,
              vector<unordered_set<root_status>>& all_rs_W1,
              vector<unordered_set<root_status>>& all_rs_W3,
              vector<unordered_map<root_status, map<size_t, bool>>>& map_W1,
              vector<unordered_map<root_status, map<size_t, bool>>>& map_W3,
              vector_2D<PathTrie>& all_trie_node,
              size_t l_min,  //
              size_t l_max,  // for tree_component, l_min=l_max = height(T)-k*,
                             // else lmin=1, lmax=
              size_t delta_i =
                  2  // >= 2, the number of neighbour of the core vertex in core
) {
  size_t M = IN_INFO.num_kind_atoms;

  for (size_t l = l_min; l <= l_max; ++l) {
    if (l == 0) continue;

    _start_time = globalTimeKeeper::tt.diff();
    _stop = false;

    size_t l2 = 0;
    size_t l1 = l - 1;

    map_W3[l].clear();

    size_t vector_size = 0;

    for (auto& rs1 : all_rs_W1[l1]) {
      if (_stop) break;
      for (auto& rs3 : all_rs_W3[l2]) {
        if (_stop) break;
        for (size_t k = 1; k <= 3; ++k) {
          if (_stop) break;
          size_t a1 = rs1.color;
          size_t a2 = rs3.color;
          size_t d1 = rs1.deg;
          size_t d2 = rs3.deg;
          size_t m1 = rs1.val;
          size_t m2 = rs3.val;

          atom_degree ad1(a1, d1);
          atom_degree ad2(a2, d2);  // direction ad2->ad1

          // if (DEBUG) cout << "a1 = " << a1 << " d1 = " << d1 << " m1 = " <<
          // m1 << " a2 = " << a2 << " d2 = " << d2 << " m2 = " << m2 << " k = "
          // << k << endl;

          if (m1 + k <= IN_INFO.atoms[a1].VALENCE &&
              m2 + k + delta_i <= IN_INFO.atoms[a2].VALENCE) {
            auto& set_left = map_W1[l1].at(rs1);
            auto& set_right = map_W3[l2].at(rs3);

            rv_edge e_tmp("in", ad2, ad1, k);  // direction: (a2,d2)->(a1,d1)
            size_t e_tmp_ind = IN_INFO.rv_edge_map.at(e_tmp);

            for (auto& _w1 : set_left) {
              if (_stop) break;

              auto w1 =
                  traverse(_w1.first, all_trie_node[l1], IN_INFO.rv1D_size);
              for (auto& _w2 : set_right) {
                if (globalTimeKeeper::tt.diff() - _start_time > time_limit) {
                  // time limit reached.
                  // cout << "len = " << l << ". Time over, limit is " <<
                  // time_limit << "s." << endl;
                  _stop = true;
                }
                if (UB_limit != 0 && vector_size > UB_limit) {
                  _stop = true;
                }
                if (_stop) break;

                auto w2 =
                    traverse(_w2.first, all_trie_node[l2], IN_INFO.rv1D_size);

                auto new_rv = SUMinus3D(w1, w2, true);

                ++new_rv[e_tmp_ind];

                if (TestLessEqAB(new_rv,
                                 IN_INFO.rv1D)  // w1 + w2 + 1_gamma <= x_star
                ) {
                  root_status _rs(rs3);
                  _rs.val += k;

                  all_rs_W3[l].emplace(_rs);
                  size_t ind = find_trie_index(new_rv, all_trie_node[l]);

                  // auto& seq1_vector = _w1.second.first;
                  // auto& seq2_vector = _w2.second.first;

                  if (map_W3[l][_rs].find(ind) == map_W3[l][_rs].end()) {
                    ++vector_size;
                    // vector <map_rv> tmp;
                    // tmp.clear();
                    map_W3[l][_rs][ind] = false;
                  }

                  // auto& _seq = map_W3[l][_rs][ind].first;
                  // auto& _num = map_W3[l][_rs][ind].second;

                  // for (auto& seq1 : seq1_vector){
                  // 	for (auto& seq2 : seq2_vector){
                  // 		if (seq_limit == 0 || _seq.size() < seq_limit){
                  // 			_seq.emplace_back(seq2, seq1, k, false);
                  // 		}
                  // 	} // for seq2
                  // }// for seq1
                  // _num += _w1.second.second * _w2.second.second;
                }  // if (TestLessEqAB)
              }    // for _w2
            }      // for _w1
          }        // if
        }          // for k
      }            // for rs3
    }              // for rs1

    if (_DEBUG)
      cout << "rooted core-subtrees l = " << l << " l_min = " << l_min
           << " l_max = " << l_max << " vector_size = " << vector_size
           << " time = " << globalTimeKeeper::tt.diff() - _start_time << endl;
  }
}

void merge_W3_tree(
    const input_info& IN_INFO, unordered_set<half_path_status>& all_rs_W2,
    vector<unordered_set<root_status>>& all_rs_W3,
    unordered_map<half_path_status, map<size_t, bool>>& map_W2,
    vector<unordered_map<root_status, map<size_t, bool>>>& map_W3,
    vector_2D<PathTrie>& all_trie_node,
    size_t l_min,  //
    size_t l_max,  // for tree_component, l_min=l_max = height(T)-k*, else
                   // lmin=1, lmax=
    size_t delta_i =
        2  // >= 2, the number of neighbour of the core vertex in core
) {
  size_t M = IN_INFO.num_kind_atoms;

  for (size_t l = l_min; l <= l_max; ++l) {
    if (l == 0) continue;

    _start_time = globalTimeKeeper::tt.diff();
    _stop = false;

    size_t l2 = l - 1;
    size_t l1 = 0;

    map_W3[l].clear();

    // cout << "l = "  << l << " l_min = " << l_min << " l_max = " << l_max  <<
    // endl;

    size_t vector_size = 0;

    for (auto& rs1 : all_rs_W2) {
      if (_stop) break;
      for (auto& rs3 : all_rs_W3[l2]) {
        if (_stop) break;
        for (size_t k = 1; k <= 3; ++k) {
          if (_stop) break;
          size_t a1 = rs1.lrs.color;
          size_t a2 = rs3.color;
          size_t d1 = rs1.lrs.deg;
          size_t d2 = rs3.deg;
          size_t m1 = rs1.lrs.val;
          size_t m2 = rs3.val;

          atom_degree ad1(a1, d1);
          atom_degree ad2(a2, d2);  // direction ad2->ad1

          // if (_DEBUG) cout << "a1 = " << a1 << " d1 = " << d1 << " m1 = " <<
          // m1 << " a2 = " << a2 << " d2 = " << d2 << " m2 = " << m2 << " k = "
          // << k << endl;

          if (((l2 > 0) && m1 + k + 1 <= IN_INFO.atoms[a1].VALENCE &&
               m2 + k <= IN_INFO.atoms[a2].VALENCE) ||
              ((l2 == 0) && m1 + k + 1 <= IN_INFO.atoms[a1].VALENCE &&
               m2 + k + delta_i <= IN_INFO.atoms[a2].VALENCE)) {
            auto& set_left = map_W2.at(rs1);
            auto& set_right = map_W3[l2].at(rs3);

            // if (_DEBUG) cout << "A a1 = " << a1 << " d1 = " << d1 << " m1 = "
            // << m1 << " a2 = " << a2 << " d2 = " << d2 << " m2 = " << m2 << "
            // k = " << k << endl;
            rv_edge e_tmp("in", ad2, ad1, k);  // direction: (a2,d2)->(a1,d1)
            size_t e_tmp_ind = IN_INFO.rv_edge_map.at(e_tmp);

            for (auto& _w1 : set_left) {
              if (_stop) break;

              auto w1 =
                  traverse(_w1.first, all_trie_node[l1], IN_INFO.rv1D_size);
              for (auto& _w2 : set_right) {
                if (globalTimeKeeper::tt.diff() - _start_time > time_limit) {
                  // time limit reached.
                  // cout << "len = " << l << ". Time over, limit is " <<
                  // time_limit << "s." << endl;
                  _stop = true;
                }
                if (UB_limit != 0 && vector_size > UB_limit) {
                  _stop = true;
                }
                if (_stop) break;

                auto w2 =
                    traverse(_w2.first, all_trie_node[l2], IN_INFO.rv1D_size);

                auto new_rv = SUMinus3D(w1, w2, true);

                ++new_rv[e_tmp_ind];

                if (TestLessEqAB(new_rv,
                                 IN_INFO.rv1D)  // w1 + w2 + 1_gamma <= x_star
                ) {
                  root_status _rs(rs1.lrs);
                  _rs.val += k;

                  all_rs_W3[l].emplace(_rs);
                  size_t ind = find_trie_index(new_rv, all_trie_node[l]);

                  // auto& seq1_vector = _w1.second.first;
                  // auto& seq2_vector = _w2.second.first;

                  if (map_W3[l][_rs].find(ind) == map_W3[l][_rs].end()) {
                    ++vector_size;
                    // vector <map_rv> tmp;
                    // tmp.clear();
                    map_W3[l][_rs][ind] = false;
                  }

                  // auto& _seq = map_W3[l][_rs][ind].first;
                  // auto& _num = map_W3[l][_rs][ind].second;

                  // for (auto& seq1 : seq1_vector){
                  // 	for (auto& seq2 : seq2_vector){
                  // 		if (seq_limit == 0 || _seq.size() < seq_limit){
                  // 			_seq.emplace_back(seq2, seq1, k, false);
                  // 		}
                  // 	} // for seq2
                  // }// for seq1
                  // _num += _w1.second.second * _w2.second.second;
                }  // if (TestLessEqAB)
              }    // for _w2
            }      // for _w1
          }        // if
        }          // for k
      }            // for rs3
    }              // for rs1

    if (_DEBUG)
      cout << "rooted core-subtrees l = " << l << " l_min = " << l_min
           << " l_max = " << l_max << " vector_size = " << vector_size
           << " time = " << globalTimeKeeper::tt.diff() - _start_time << endl;
  }
}

// A function to find feasible pairs
void combine_check_tree_component(
    const input_info& IN_INFO, vector<unordered_set<root_status>>& all_rs_W1,
    vector<unordered_set<root_status>>& all_rs_W3,
    vector<unordered_map<root_status, map<size_t, bool>>>& map_W1,
    vector<unordered_map<root_status, map<size_t, bool>>>& map_W3,
    vector_2D<PathTrie>& all_trie_node,
    // vector <map_rv>& all_seq_h,
    // vector <map_rv>& all_seq_tree,
    vector<feasible_pair_W>& all_fp_W, size_t l1, size_t l3,
    root_status& target_rs,
    // const size_t& core_height,
    // size_t& count_num_h,
    size_t& count_num_tree) {
  size_t M = IN_INFO.num_kind_atoms;
  size_t count = 0;
  // size_t count_num = 0;

  all_fp_W.clear();

  // size_t l1 = 0;
  // size_t l3 = l - 1;

  _start_time = globalTimeKeeper::tt.diff();
  _stop = false;

  for (auto& rs3 : all_rs_W3[l3]) {
    if (_stop) break;
    auto& set_w1 = map_W3[l3].at(rs3);
    if (l3 == 0 && (rs3.color != target_rs.color || rs3.deg != target_rs.deg)) {
      continue;
    }
    for (auto& rs1 : all_rs_W1[l1]) {
      auto& _set_w2 = map_W1[l1].at(rs1);
      if (_stop) break;
      for (size_t k = 1; k <= 3; ++k) {
        if (_stop) break;
        size_t a1 = rs3.color;
        size_t d1 = rs3.deg;
        size_t m1 = rs3.val;

        size_t a2 = rs1.color;
        size_t d2 = rs1.deg;
        size_t m2 = rs1.val;

        atom_degree ad1(a1, d1);
        atom_degree ad2(a2, d2);

        if (m1 + k <= IN_INFO.atoms[a1].VALENCE &&
            m2 + k <= IN_INFO.atoms[a2].VALENCE &&
            ((l3 != 0) || (l3 == 0 && m1 + k == target_rs.val))) {
          rv_edge e_tmp("in", ad1, ad2, k);
          size_t e_tmp_ind = IN_INFO.rv_edge_map.at(e_tmp);

          set<resource_vector_1D> set_w2;
          for (auto& _rv : _set_w2) {
            auto _w2 =
                traverse(_rv.first, all_trie_node[l1], IN_INFO.rv1D_size);
            set_w2.insert(_w2);
          }

          for (auto& _w1 : set_w1) {
            // if (globalTimeKeeper::tt.diff() - _start_time > time_limit){

            // 	// time limit reached.
            // 	// cout << "Time over, limit is " << time_limit << "s." << endl;
            // 	_stop = true;
            // }
            if (_stop) break;

            auto w1 = traverse(_w1.first, all_trie_node[l3], IN_INFO.rv1D_size);

            ++w1[e_tmp_ind];

            if (TestLessEqAB(w1, IN_INFO.rv1D)) {
              auto w2 = SUMinus3D(IN_INFO.rv1D, w1, false);

              if (set_w2.find(w2) != set_w2.end()) {
                ++count;
                size_t ind_w2 = find_trie_index(w2, all_trie_node[l1]);
                auto& _w2 = map_W1[l1].at(rs1).at(ind_w2);
                // cout << "count = " << count << " _w1.first = " << _w1.first
                // << endl; size_t tmp = _w1.second.second * _w2.second;

                // auto& seq1 = _w1.second.first;
                // auto& seq2 = _w2.first;

                // count_num_tree += tmp;

                // auto& seq1_vector = _w1.second.first;
                // auto& seq2_vector = _w2.first;

                // for (auto& seq1 : seq1_vector){
                // 	for (auto& seq2 : seq2_vector){
                // 		if (pair_limit == 0 || count <= pair_limit){
                // 			all_seq_tree.emplace_back(seq1, seq2, k,
                // false);
                // 		}
                // 	} // for seq2
                // }// for seq1
                _w1.second = true;
                _w2 = true;

                // cout << "rs3.color = " << rs3.color << " _w1.first = " <<
                // _w1.first << endl;
                all_fp_W.emplace_back(rs3, rs1, l3, _w1.first, l1, ind_w2, k);
              }
            }

          }  // for _w1
        }    // if (v1 + k <= IN_INFO.atoms[a1].VALENCE && v2 + k <=
             // IN_INFO.atoms[a2].VALENCE)
      }      // for k
    }        // for rs1
  }          // for rs4
  // if (_DEBUG) cout << "l1 = " << l1 << " count = " << count << endl;

  // size_t _size = 0;
  // for (auto& temp : all_seq){
  // 	_size += temp.size();
  // }
  // _size = all_seq.size();

  // if (_DEBUG) cout << "Number of feasible pairs (h = core_height) = " <<
  // count_h << endl; if (_DEBUG) cout << "Number of feasible pairs (others) = "
  // << count - count_h << endl;

  if (_DEBUG) cout << "Number of feasible pairs = " << count << endl;
  // cout << "A lower bound on the number of graphs = " << count_num << endl;
  // cout << "Number of generated graphs = " << _size << endl;

  if (_EXP) feasible_pairs_sum += count;
}

//*** IDO begin ***//

void check_tree_component(
    const input_info& IN_INFO,
    // vector <unordered_set <root_status>>& all_rs_W1,
    vector<unordered_set<root_status>>& all_rs_W3,
    // vector <unordered_map <root_status, map <size_t, bool>>>& map_W1,
    vector<unordered_map<root_status, map<size_t, bool>>>& map_W3,
    vector_2D<PathTrie>& all_trie_node,
    // vector <map_rv>& all_seq_h,
    // vector <map_rv>& all_seq_tree,
    vector<DAG_node>& all_fp_W,
    // size_t l1,
    size_t l3, root_status& target_rs) {
  // const size_t& core_height,
  // size_t& count_num_h,
  // size_t& count_num_tree){

  for (auto& rs3 : all_rs_W3[l3]) {
    if (target_rs == rs3) {
      auto& set_w3 = map_W3[l3].at(rs3);
      for (auto& _w3 : set_w3) {
        auto w3 = traverse(_w3.first, all_trie_node[l3], IN_INFO.rv1D_size);
        if (TestEqAB(w3, IN_INFO.rv1D)) {
          // DAG_nodeに加える
          all_fp_W.emplace_back(rs3, l3, _w3.first);
          // DAG
          _w3.second = true;
        }
      }
    }
  }
}

//*** IDO end ***//

void merge_W1_DAG(const input_info& IN_INFO,
                  vector<unordered_set<root_status>>& all_rs_W1,
                  unordered_set<half_path_status>& all_rs_W2,
                  vector<unordered_map<root_status, map<size_t, bool>>>& map_W1,
                  unordered_map<half_path_status, map<size_t, bool>>& map_W2,
                  unordered_map<half_path_status, vector<size_t>>& map_W2_ind,
                  vector_2D<PathTrie>& all_trie_node,
                  unordered_map<DAG_node, vector<DAG_edge>>& map_DAG_W1,
                  size_t len  // max core height - 1
) {
  size_t M = IN_INFO.num_kind_atoms;

  // for  (size_t l = 1; l <= len; ++l){
  for (size_t l = len; l >= 1; --l) {
    _start_time = globalTimeKeeper::tt.diff();
    _stop = false;

    size_t l1 = l - 1;
    size_t l2 = 0;

    // map_W1[l].clear();

    size_t vector_size = 0;

    for (auto& rhps : all_rs_W2) {
      if (_stop) break;
      for (auto& rs : all_rs_W1[l1]) {
        if (_stop) break;
        for (size_t k = 1; k <= 3; ++k) {
          if (_stop) break;
          size_t a1 = rs.color;
          size_t a2 = rhps.lrs.color;
          size_t d1 = rs.deg;
          size_t d2 = rhps.lrs.deg;
          size_t m1 = rs.val;
          size_t m2 = rhps.lrs.val;

          atom_degree ad1(a1, d1);
          atom_degree ad2(a2, d2);

          if (m1 + k <= IN_INFO.atoms[a1].VALENCE &&
              m2 + k + 1 <= IN_INFO.atoms[a2].VALENCE) {
            auto& set_left = map_W1[l1].at(rs);
            auto& set_right = map_W2_ind.at(rhps);

            rv_edge e_tmp("in", ad2, ad1, k);  // direction: (a2,d2)->(a1,d1)
            size_t e_tmp_ind = IN_INFO.rv_edge_map.at(e_tmp);

            size_t index_size = map_W2_ind.at(rhps).size();
            size_t start_index = (index_size - 1) * (l - 1) / len;
            size_t ind = start_index;
            bool flag = true;

            while (ind != start_index || flag) {
              flag = false;
              if (ind >= index_size) {
                ind -= index_size;
              }
              if (_stop) break;
              auto& _w2 = map_W2_ind[rhps][ind];
              ++ind;
              if (ind == index_size) {
                ind = 0;
              }

              auto w2 = traverse(_w2, all_trie_node[l2], IN_INFO.rv1D_size);
              auto& _w2_second = map_W2.at(rhps).at(_w2);
              for (auto& _w1 : set_left) {
                if (globalTimeKeeper::tt.diff() - _start_time > time_limit) {
                  // time limit reached.
                  // cout << "len = " << l << ". Time over, limit is " <<
                  // time_limit << "s." << endl;
                  _stop = true;
                }
                if (UB_limit != 0 && vector_size > UB_limit) {
                  _stop = true;
                }
                if (_stop) break;

                auto w1 =
                    traverse(_w1.first, all_trie_node[l1], IN_INFO.rv1D_size);
                auto new_rv = SUMinus3D(w1, w2, true);

                ++new_rv[e_tmp_ind];

                if (TestLessEqAB(
                        new_rv,
                        IN_INFO.rv1D)  // w1 + w2 + 1_mu + 1_gamma <= x_star
                ) {
                  root_status _rs(rhps.rrs);
                  _rs.val += k;

                  // all_rs_W1[l].emplace(_rs);
                  size_t ind = find_trie_index(new_rv, all_trie_node[l]);

                  // auto& seq1_vector = _w1.second.first;
                  // auto& seq2_vector = _w2_second.first;

                  if (map_W1[l][_rs].find(ind) != map_W1[l][_rs].end()) {
                    if (map_W1[l][_rs][ind]) {
                      _w1.second = true;
                      DAG_node tmp(_rs, l, ind);
                      map_DAG_W1[tmp].emplace_back(rs, l1, _w1.first, a2, _w2,
                                                   k);
                    }
                  }

                }  // if (TestLessEqAB)
              }    // for _w1
            }      // while
          }        // if
        }          // for k
      }            // for rs
    }              // for rhps

  }  // for l
}

void merge_W3_DAG(
    const input_info& IN_INFO, vector<unordered_set<root_status>>& all_rs_W1,
    vector<unordered_set<root_status>>& all_rs_W3,
    vector<unordered_map<root_status, map<size_t, bool>>>& map_W1,
    vector<unordered_map<root_status, map<size_t, bool>>>& map_W3,
    vector_2D<PathTrie>& all_trie_node,
    unordered_map<DAG_node, vector<DAG_edge>>& map_DAG_W3,
    size_t l_min,  //
    size_t l_max,  // for tree_component, l_min=l_max = height(T)-k*, else
                   // lmin=1, lmax=
    size_t delta_i =
        2  // >= 2, the number of neighbour of the core vertex in core
) {
  size_t M = IN_INFO.num_kind_atoms;

  // for (size_t l = l_min; l <= l_max; ++l){
  for (size_t l = l_max; l >= l_min; --l) {
    if (l == 0) break;

    _start_time = globalTimeKeeper::tt.diff();
    _stop = false;

    size_t l2 = 0;
    size_t l1 = l - 1;

    // map_W3[l].clear();

    size_t vector_size = 0;

    for (auto& rs1 : all_rs_W1[l1]) {
      if (_stop) break;
      for (auto& rs3 : all_rs_W3[l2]) {
        if (_stop) break;
        for (size_t k = 1; k <= 3; ++k) {
          if (_stop) break;
          size_t a1 = rs1.color;
          size_t a2 = rs3.color;
          size_t d1 = rs1.deg;
          size_t d2 = rs3.deg;
          size_t m1 = rs1.val;
          size_t m2 = rs3.val;

          atom_degree ad1(a1, d1);
          atom_degree ad2(a2, d2);  // direction ad2->ad1

          // if (DEBUG) cout << "a1 = " << a1 << " d1 = " << d1 << " m1 = " <<
          // m1 << " a2 = " << a2 << " d2 = " << d2 << " m2 = " << m2 << " k = "
          // << k << endl;

          if (m1 + k <= IN_INFO.atoms[a1].VALENCE &&
              m2 + k + delta_i <= IN_INFO.atoms[a2].VALENCE) {
            auto& set_left = map_W1[l1].at(rs1);
            auto& set_right = map_W3[l2].at(rs3);

            rv_edge e_tmp("in", ad2, ad1, k);  // direction: (a2,d2)->(a1,d1)
            size_t e_tmp_ind = IN_INFO.rv_edge_map.at(e_tmp);

            for (auto& _w1 : set_left) {
              if (_stop) break;

              auto w1 =
                  traverse(_w1.first, all_trie_node[l1], IN_INFO.rv1D_size);
              for (auto& _w2 : set_right) {
                if (globalTimeKeeper::tt.diff() - _start_time > time_limit) {
                  // time limit reached.
                  // cout << "len = " << l << ". Time over, limit is " <<
                  // time_limit << "s." << endl;
                  _stop = true;
                }
                if (UB_limit != 0 && vector_size > UB_limit) {
                  _stop = true;
                }
                if (_stop) break;

                auto w2 =
                    traverse(_w2.first, all_trie_node[l2], IN_INFO.rv1D_size);

                auto new_rv = SUMinus3D(w1, w2, true);

                ++new_rv[e_tmp_ind];

                if (TestLessEqAB(new_rv,
                                 IN_INFO.rv1D)  // w1 + w2 + 1_gamma <= x_star
                ) {
                  root_status _rs(rs3);
                  _rs.val += k;

                  // all_rs_W3[l].emplace(_rs);
                  size_t ind = find_trie_index(new_rv, all_trie_node[l]);

                  // auto& seq1_vector = _w1.second.first;
                  // auto& seq2_vector = _w2.second.first;

                  if (map_W3[l][_rs].find(ind) != map_W3[l][_rs].end()) {
                    if (map_W3[l][_rs][ind]) {
                      _w1.second = true;
                      DAG_node tmp(_rs, l, ind);
                      map_DAG_W3[tmp].emplace_back(rs1, l1, _w1.first, a2,
                                                   _w2.first, k);
                    }
                  }
                }  // if (TestLessEqAB)
              }    // for _w2
            }      // for _w1
          }        // if
        }          // for k
      }            // for rs3
    }              // for rs1
  }
}

void merge_W3_tree_DAG(
    const input_info& IN_INFO, unordered_set<half_path_status>& all_rs_W2,
    vector<unordered_set<root_status>>& all_rs_W3,
    unordered_map<half_path_status, map<size_t, bool>>& map_W2,
    vector<unordered_map<root_status, map<size_t, bool>>>& map_W3,
    vector_2D<PathTrie>& all_trie_node,
    unordered_map<DAG_node, vector<DAG_edge>>& map_DAG_W3,
    size_t l_min,  //
    size_t l_max,  // for tree_component, l_min=l_max = height(T)-k*, else
                   // lmin=1, lmax=
    size_t delta_i =
        2  // >= 2, the number of neighbour of the core vertex in core
) {
  size_t M = IN_INFO.num_kind_atoms;

  // for (size_t l = l_min; l <= l_max; ++l){
  for (size_t l = l_max; l >= l_min; --l) {
    // cout << "l = " << l << endl;
    if (l == 0) break;

    _start_time = globalTimeKeeper::tt.diff();
    _stop = false;

    size_t l2 = l - 1;
    size_t l1 = 0;

    // map_W3[l].clear();

    // cout << "l = "  << l << " l_min = " << l_min << " l_max = " << l_max  <<
    // endl;

    size_t vector_size = 0;

    for (auto& rs1 : all_rs_W2) {
      if (_stop) break;
      for (auto& rs3 : all_rs_W3[l2]) {
        if (_stop) break;
        for (size_t k = 1; k <= 3; ++k) {
          if (_stop) break;
          size_t a1 = rs1.lrs.color;
          size_t a2 = rs3.color;
          size_t d1 = rs1.lrs.deg;
          size_t d2 = rs3.deg;
          size_t m1 = rs1.lrs.val;
          size_t m2 = rs3.val;

          atom_degree ad1(a1, d1);
          atom_degree ad2(a2, d2);  // direction ad2->ad1

          // if (_DEBUG) cout << "a1 = " << a1 << " d1 = " << d1 << " m1 = " <<
          // m1 << " a2 = " << a2 << " d2 = " << d2 << " m2 = " << m2 << " k = "
          // << k << endl;

          if (((l2 > 0) && m1 + k + 1 <= IN_INFO.atoms[a1].VALENCE &&
               m2 + k <= IN_INFO.atoms[a2].VALENCE) ||
              ((l2 == 0) && m1 + k + 1 <= IN_INFO.atoms[a1].VALENCE &&
               m2 + k + delta_i <= IN_INFO.atoms[a2].VALENCE)) {
            auto& set_left = map_W2.at(rs1);
            auto& set_right = map_W3[l2].at(rs3);

            // if (_DEBUG) cout << "A a1 = " << a1 << " d1 = " << d1 << " m1 = "
            // << m1 << " a2 = " << a2 << " d2 = " << d2 << " m2 = " << m2 << "
            // k = " << k << endl;
            rv_edge e_tmp("in", ad2, ad1, k);  // direction: (a2,d2)->(a1,d1)
            size_t e_tmp_ind = IN_INFO.rv_edge_map.at(e_tmp);

            for (auto& _w1 : set_left) {
              if (_stop) break;

              auto w1 =
                  traverse(_w1.first, all_trie_node[l1], IN_INFO.rv1D_size);
              for (auto& _w2 : set_right) {
                if (globalTimeKeeper::tt.diff() - _start_time > time_limit) {
                  // time limit reached.
                  // cout << "len = " << l << ". Time over, limit is " <<
                  // time_limit << "s." << endl;
                  _stop = true;
                }
                if (UB_limit != 0 && vector_size > UB_limit) {
                  _stop = true;
                }
                if (_stop) break;

                auto w2 =
                    traverse(_w2.first, all_trie_node[l2], IN_INFO.rv1D_size);

                auto new_rv = SUMinus3D(w1, w2, true);

                ++new_rv[e_tmp_ind];

                if (TestLessEqAB(new_rv,
                                 IN_INFO.rv1D)  // w1 + w2 + 1_gamma <= x_star
                ) {
                  root_status _rs(rs1.lrs);
                  _rs.val += k;

                  // all_rs_W3[l].emplace(_rs);
                  size_t ind = find_trie_index(new_rv, all_trie_node[l]);

                  // auto& seq1_vector = _w1.second.first;
                  // auto& seq2_vector = _w2.second.first;

                  if (map_W3[l][_rs].find(ind) != map_W3[l][_rs].end()) {
                    if (map_W3[l][_rs][ind]) {
                      _w2.second = true;
                      DAG_node tmp(_rs, l, ind);
                      map_DAG_W3[tmp].emplace_back(rs3, l2, _w2.first, a1,
                                                   _w1.first, k);
                    }
                  }
                }  // if (TestLessEqAB)
              }    // for _w2
            }      // for _w1
          }        // if
        }          // for k
      }            // for rs3
    }              // for rs1
  }
}

void GenT1(const input_info& IN_INFO,
           vector<unordered_set<root_status_T>>& all_rs_T1,
           vector<unordered_map<root_status_T, map<size_t, bool>>>& map_T1,
           vector<PathTrie>& all_trie_node_path,

           // root_status_height& _rsh
           root_status& _rs) {
  // root_status _rs(_rsh.color, _rsh.val, _rsh.deg);
  root_status_height _rsh(_rs, 0);
  root_status_T rsT = make_pair(_rs, _rsh);
  all_rs_T1[0].insert(rsT);

  resource_vector null_rv(IN_INFO.num_kind_atom_degree);
  resource_vector_1D rv1D = rv2rv1D(null_rv, IN_INFO);

  size_t ind = find_trie_index(rv1D, all_trie_node_path);

  if (map_T1[0][rsT].find(ind) == map_T1[0][rsT].end()) {
    if (_EXP) bi_rooted_core_subtrees_sum += 1;
    map_rv tmp(_rsh.color, 0);
    // map_T1[0][rsT][ind] = make_pair(vector <map_tree> ({map_tree(tmp)}), 1);
    map_T1[0][rsT][ind] = false;
  }
}

void GenT2_0(const input_info& IN_INFO,
             unordered_set<root_status_height>& all_rs_T2,
             unordered_map<root_status_height, map<size_t, bool>>& map_T2,
             unordered_map<root_status_height,
                           map<size_t, vector<DAG_node_height>>>& map_T2_W3,
             unordered_map<root_status_height, vector<size_t>>& map_T2_ind,
             vector<PathTrie>& all_trie_node_path,

             unordered_set<root_status_height>& all_rs_W2_core,
             unordered_map<root_status_height, map<size_t, bool>>& map_W2_core,
             vector<unordered_set<root_status>>& all_rs_W3,
             vector<unordered_map<root_status, map<size_t, bool>>>& map_W3,
             vector_2D<PathTrie>& all_trie_node,

             size_t chUB,

             bool exact = false

) {
  size_t max_h = 0;
  if (chUB >= IN_INFO.k_star) {
    max_h = chUB - IN_INFO.k_star;
  } else {
    max_h = 0;
  }
  size_t min_h = 0;
  size_t vector_size = 0;

  for (size_t h = min_h; h <= max_h; ++h) {
    if (h == 0) {
      for (auto& rsh : all_rs_W2_core) {
        if (rsh.height > chUB) {
          continue;
        }
        // cout << "rsh.color = "  << rsh.color << " rsh.deg = " << rsh.deg << "
        // rsh.val = " <<rsh.val << " rsh.height = " << rsh.height << endl;
        for (auto& _w : map_W2_core.at(rsh)) {
          auto w = traverse(_w.first, all_trie_node[0], IN_INFO.rv1D_size);

          if (exact && IN_INFO.rv1D != w) {
            continue;
          }

          all_rs_T2.insert(rsh);

          size_t ind = find_trie_index(w, all_trie_node_path);

          if (map_T2[rsh].find(ind) == map_T2[rsh].end()) {
            ++vector_size;
            // vector <map_tree> tmp;
            // tmp.clear();
            // for (auto& _seq : _w.second.first){
            // 	tmp.emplace_back(_seq, rsh.height);
            // }
            vector<DAG_node_height> tmp;
            tmp.emplace_back(rsh, 0, _w.first);
            map_T2[rsh][ind] = false;
            map_T2_W3[rsh][ind] = tmp;
          } else {
            // for (auto& _seq: _w.second.first){
            // 	if (seq_limit == 0 || map_T2[rsh][ind].first.size() <
            // seq_limit){ map_T2[rsh][ind].first.emplace_back(_seq,
            // rsh.height);
            // 	}
            // }
            map_T2_W3[rsh][ind].emplace_back(rsh, 0, _w.first);
            // map_T2[rsh][ind].second += _w.second.second;
          }
        }
      }
    } else {
      for (auto& rs : all_rs_W3[h]) {
        for (auto& _w : map_W3[h].at(rs)) {
          auto w = traverse(_w.first, all_trie_node[h], IN_INFO.rv1D_size);

          if (exact && IN_INFO.rv1D != w) {
            continue;
          }

          root_status_height rsh(rs, h + IN_INFO.k_star);

          all_rs_T2.emplace(rs, h + IN_INFO.k_star);

          size_t ind = find_trie_index(w, all_trie_node_path);

          if (map_T2[rsh].find(ind) == map_T2[rsh].end()) {
            ++vector_size;
            // vector <map_tree> tmp;
            // tmp.clear();
            // for (auto& _seq : _w.second.first){
            // 	tmp.emplace_back(_seq, rsh.height);
            // }
            vector<DAG_node_height> tmp;
            tmp.emplace_back(rsh, h, _w.first);
            map_T2[rsh][ind] = false;
            map_T2_W3[rsh][ind] = tmp;
          } else {
            // for (auto& _seq: _w.second.first){
            // 	if (seq_limit == 0 || map_T2[rsh][ind].first.size() <
            // seq_limit){ map_T2[rsh][ind].first.emplace_back(_seq,
            // rsh.height);
            // 	}
            // }
            map_T2_W3[rsh][ind].emplace_back(rsh, h, _w.first);
            // map_T2[rsh][ind].second += _w.second.second;
          }
        }
      }
    }
  }

  if (_DEBUG) cout << "sum of size of core-subtrees = " << vector_size << endl;
  if (_EXP) rooted_core_subtrees_sum += vector_size;

  for (auto& rsh : all_rs_T2) {
    map_T2_ind[rsh].clear();
    for (auto& tmp : map_T2.at(rsh)) {
      map_T2_ind[rsh].push_back(tmp.first);
    }
  }
}

void merge_T1(const input_info& IN_INFO,
              vector<unordered_set<root_status_T>>& all_rs_T1,
              unordered_set<root_status_height>& all_rs_T2,
              vector<unordered_map<root_status_T, map<size_t, bool>>>& map_T1,
              unordered_map<root_status_height, map<size_t, bool>>& map_T2,
              unordered_map<root_status_height, vector<size_t>>& map_T2_ind,
              vector<PathTrie>& all_trie_node_path, size_t len) {
  size_t M = IN_INFO.num_kind_atoms;

  for (size_t l = 1; l <= len; ++l) {
    _start_time = globalTimeKeeper::tt.diff();
    _stop = false;

    size_t l1 = l - 1;
    size_t l2 = 0;

    map_T1[l].clear();

    size_t vector_size = 0;

    for (auto& rhps : all_rs_T2) {
      if (_stop) break;
      for (auto& rs : all_rs_T1[l1]) {
        if (_stop) break;
        for (size_t k = 1; k <= 3; ++k) {
          if (_stop) break;
          size_t a1 = rs.second.color;
          size_t a2 = rhps.color;
          size_t d1 = rs.second.deg;
          size_t d2 = rhps.deg;
          size_t m1 = rs.second.val;
          size_t m2 = rhps.val;
          size_t h1 = rs.second.height;
          size_t h2 = rhps.height;

          atom_degree ad1(a1, d1);
          atom_degree ad2(a2, d2);  // no direction in core

          size_t h = h1;
          if (h2 > h) h = h2;

          if ((l1 > 0 && m1 + k <= IN_INFO.atoms[a1].VALENCE &&
               m2 + k + 1 <= IN_INFO.atoms[a2].VALENCE) ||
              (l1 == 0 && k == rs.first.val &&
               m2 + k + 1 <= IN_INFO.atoms[a2].VALENCE)) {
            auto& set_left = map_T1[l1].at(rs);
            auto& set_right = map_T2_ind.at(rhps);

            rv_edge e_tmp;
            if (IN_INFO.atom_deg_map.at(ad1) >= IN_INFO.atom_deg_map.at(ad2)) {
              e_tmp = rv_edge("co", ad1, ad2, k);
            } else {
              e_tmp = rv_edge("co", ad2, ad1, k);
            }
            size_t e_tmp_ind = IN_INFO.rv_edge_map.at(e_tmp);

            size_t index_size = map_T2_ind.at(rhps).size();
            size_t start_index = (index_size - 1) * (l - 1) / len;
            size_t ind = start_index;
            bool flag = true;

            while (ind != start_index || flag) {
              flag = false;
              if (ind >= index_size) {
                ind -= index_size;
              }
              if (_stop) break;
              auto& _w2 = map_T2_ind[rhps][ind];
              // if (_DEBUG) cout << "ind = " <<ind << " index_size = " <<
              // index_size  << endl;
              ++ind;
              if (ind == index_size) {
                ind = 0;
              }

              auto w2 = traverse(_w2, all_trie_node_path, IN_INFO.rv1D_size);
              auto& _w2_second = map_T2.at(rhps).at(_w2);
              for (auto& _w1 : set_left) {
                if (globalTimeKeeper::tt.diff() - _start_time > time_limit) {
                  // time limit reached.
                  // cout << "len = " << l << ". Time over, limit is " <<
                  // time_limit << "s." << endl;
                  _stop = true;
                }
                if (UB_limit != 0 && vector_size > UB_limit) {
                  _stop = true;
                }
                if (_stop) break;

                auto w1 =
                    traverse(_w1.first, all_trie_node_path, IN_INFO.rv1D_size);
                auto new_rv = SUMinus3D(w1, w2, true);

                ++new_rv[e_tmp_ind];

                if (TestLessEqAB(
                        new_rv,
                        IN_INFO.rv1D)  // w1 + w2 + 1_mu + 1_gamma <= x_star
                ) {
                  root_status_T _rsT = make_pair(rs.first, rhps);
                  _rsT.second.val += k;
                  if (rs.second.height > _rsT.second.height) {
                    _rsT.second.height = rs.second.height;
                  }

                  all_rs_T1[l].emplace(_rsT);
                  size_t ind = find_trie_index(new_rv, all_trie_node_path);

                  // auto& seq1_vector = _w1.second.first;
                  // auto& seq2_vector = _w2_second.first;

                  if (map_T1[l][_rsT].find(ind) == map_T1[l][_rsT].end()) {
                    ++vector_size;
                    // vector <map_tree> tmp;
                    // tmp.clear();
                    map_T1[l][_rsT][ind] = false;
                  }

                  // auto& _seq = map_T1[l][_rsT][ind].first;
                  // auto& _num = map_T1[l][_rsT][ind].second;

                  // for (auto& seq1 : seq1_vector){
                  // 	for (auto& seq2 : seq2_vector){
                  // 		if (seq_limit == 0 || _seq.size() < seq_limit){
                  // 			_seq.emplace_back(seq1, seq2, k, h,
                  // true);
                  // 		}
                  // 	} // for seq2
                  // }// for seq1
                  // _num += _w1.second.second * _w2_second.second;
                }  // if (TestLessEqAB)
              }    // for _w1
            }      // while
          }        // if
        }          // for k
      }            // for rs
    }              // for rhps

    if (_DEBUG)
      cout << "bi-rooted core-subtrees l = " << l << " len = " << len
           << " vector_size = " << vector_size
           << " time = " << globalTimeKeeper::tt.diff() - _start_time << endl;
    if (_EXP) bi_rooted_core_subtrees_sum += vector_size;
  }  // for l
}

// A function to find feasible pairs
void combine_check_edge_component(
    const input_info& IN_INFO, vector<unordered_set<root_status_T>>& all_rs_T1,
    vector<unordered_map<root_status_T, map<size_t, bool>>>& map_T1,
    vector<PathTrie>& all_trie_node_path,
    // vector <map_tree>& all_seq_h,
    // vector <map_tree>& all_seq,
    vector<feasible_pair_T>& all_fp_T, size_t l1, size_t l2,
    root_status& target_rs1, root_status& target_rs2, const size_t& core_height,
    size_t& count_num_h, size_t& count_num) {
  size_t M = IN_INFO.num_kind_atoms;
  size_t count = 0;
  size_t count_h = 0;
  // size_t count_num = 0;

  all_fp_T.clear();

  _start_time = globalTimeKeeper::tt.diff();
  _stop = false;

  for (auto& rsT1 : all_rs_T1[l1]) {
    if (_stop) break;
    if (rsT1.first != target_rs1) continue;
    auto& rs1 = rsT1.second;
    auto& set_w1 = map_T1[l1].at(rsT1);
    for (auto& rsT2 : all_rs_T1[l2]) {
      if (_stop) break;
      if (rsT2.first != target_rs2) continue;
      auto& rs2 = rsT2.second;
      auto& _set_w2 = map_T1[l2].at(rsT2);
      for (size_t k = 1; k <= 3; ++k) {
        if (_stop) break;
        size_t a1 = rs1.color;
        size_t d1 = rs1.deg;
        size_t m1 = rs1.val;
        size_t h1 = rs1.height;

        size_t a2 = rs2.color;
        size_t d2 = rs2.deg;
        size_t m2 = rs2.val;
        size_t h2 = rs2.height;

        atom_degree ad1(a1, d1);
        atom_degree ad2(a2, d2);

        size_t h = h1;
        if (h2 > h) h = h2;

        if (l1 == 0 && m1 != k) continue;
        if (l2 == 0 && m2 != k) continue;

        if (m1 + k <= IN_INFO.atoms[a1].VALENCE &&
            m2 + k <= IN_INFO.atoms[a2].VALENCE &&
            (h >= IN_INFO.chLB && h <= IN_INFO.chUB)) {
          // if (_DEBUG) cout << "rsT1.first.color = " << rsT1.first.color << "
          // rsT1.first.deg = " << rsT1.first.deg << " rsT1.first.val = " <<
          // rsT1.first.val << endl; if (_DEBUG) cout << "rsT2.first.color = "
          // << rsT2.first.color << " rsT2.first.deg = " << rsT2.first.deg << "
          // rsT2.first.val = " << rsT2.first.val << endl;

          // if (_DEBUG) cout << "a1 = " << a1 << " d1 = " << d1 << " m1 = " <<
          // m1 << " h1 = " << h1 << endl; if (_DEBUG) cout << "a2 = " << a2 <<
          // " d2 = " << d2 << " m2 = " << m2 << " h2 = "  << h2 <<" k = " << k
          // << endl;

          rv_edge e_tmp;
          if (IN_INFO.atom_deg_map.at(ad1) >= IN_INFO.atom_deg_map.at(ad2)) {
            e_tmp = rv_edge("co", ad1, ad2, k);
          } else {
            e_tmp = rv_edge("co", ad2, ad1, k);
          }
          size_t e_tmp_ind = IN_INFO.rv_edge_map.at(e_tmp);

          set<resource_vector_1D> set_w2;
          for (auto& _rv : _set_w2) {
            auto _w2 =
                traverse(_rv.first, all_trie_node_path, IN_INFO.rv1D_size);
            set_w2.insert(_w2);
          }

          for (auto& _w1 : set_w1) {
            // if (globalTimeKeeper::tt.diff() - _start_time > time_limit){

            // 	// time limit reached.
            // 	// cout << "Time over, limit is " << time_limit << "s." << endl;
            // 	_stop = true;
            // }
            if (_stop) break;

            auto w1 =
                traverse(_w1.first, all_trie_node_path, IN_INFO.rv1D_size);

            ++w1[e_tmp_ind];

            // if (_DEBUG) print_rv1D(IN_INFO, IN_INFO.rv1D);
            // if (_DEBUG) print_rv1D(IN_INFO, w1);

            if (TestLessEqAB(w1, IN_INFO.rv1D)) {
              auto w2 = SUMinus3D(IN_INFO.rv1D, w1, false);

              if (set_w2.find(w2) != set_w2.end()) {
                ++count;
                size_t ind_w2 = find_trie_index(w2, all_trie_node_path);
                auto& _w2 = map_T1[l2].at(rsT2).at(ind_w2);

                // size_t tmp = _w1.second.second * _w2.second;

                // auto& seq1 = _w1.second.first;
                // auto& seq2 = _w2.first;

                // if (h == core_height){
                // 	// count_num_h += (tmp + 1) / 2;
                // 	count_num_h += tmp;
                // } else {
                // 	// count_num += (tmp + 1) / 2;
                // 	count_num += tmp;
                // }

                // auto& seq1_vector = _w1.second.first;
                // auto& seq2_vector = _w2.first;

                // if (h == core_height){
                // 	++count_h;
                // }

                // for (auto& seq1 : seq1_vector){
                // 	for (auto& seq2 : seq2_vector){
                // 		if (pair_limit == 0 || count <= pair_limit){
                // 			if (h == core_height){
                // 				all_seq_h.emplace_back(seq1,
                // seq2, k, h, false); 			} else {
                // all_seq.emplace_back(seq1, seq2, k, h, false);
                // 			}
                // 		}
                // 	} // for seq2
                // }// for seq1

                _w1.second = true;
                _w2 = true;

                all_fp_T.emplace_back(rsT1, rsT2, l1, _w1.first, l2, ind_w2, k);
              }
            }
          }
        }  // if
      }    // for k
    }      // for rsT2
  }        // for rsT1

  // size_t _size = 0;
  // for (auto& temp : all_seq){
  // 	_size += temp.size();
  // }
  // _size = all_seq.size();

  // if (_DEBUG) cout << "Number of feasible pairs (h = core_height) = " <<
  // count_h << endl; if (_DEBUG) cout << "Number of feasible pairs (others) = "
  // << count - count_h << endl;
  if (_DEBUG) cout << "Number of feasible pairs = " << count << endl;
  // cout << "A lower bound on the number of graphs = " << count_num << endl;
  // cout << "Number of generated graphs = " << _size << endl;
  if (_EXP) feasible_pairs_sum += count;
}

void merge_T1_DAG(
    const input_info& IN_INFO, vector<unordered_set<root_status_T>>& all_rs_T1,
    unordered_set<root_status_height>& all_rs_T2,
    vector<unordered_map<root_status_T, map<size_t, bool>>>& map_T1,
    unordered_map<root_status_height, map<size_t, bool>>& map_T2,
    unordered_map<root_status_height, map<size_t, vector<DAG_node_height>>>&
        map_T2_W3,
    unordered_map<root_status_height, vector<size_t>>& map_T2_ind,
    vector<PathTrie>& all_trie_node_path,
    unordered_map<DAG_node_T, vector<DAG_edge_T>>& map_DAG_T, size_t len) {
  size_t M = IN_INFO.num_kind_atoms;

  // for (size_t l = 1; l <= len; ++l){
  for (size_t l = len; l >= 1; --l) {
    _start_time = globalTimeKeeper::tt.diff();
    _stop = false;

    size_t l1 = l - 1;
    size_t l2 = 0;

    // map_T1[l].clear();

    size_t vector_size = 0;

    for (auto& rhps : all_rs_T2) {
      if (_stop) break;
      for (auto& rs : all_rs_T1[l1]) {
        if (_stop) break;
        for (size_t k = 1; k <= 3; ++k) {
          if (_stop) break;
          size_t a1 = rs.second.color;
          size_t a2 = rhps.color;
          size_t d1 = rs.second.deg;
          size_t d2 = rhps.deg;
          size_t m1 = rs.second.val;
          size_t m2 = rhps.val;
          size_t h1 = rs.second.height;
          size_t h2 = rhps.height;

          atom_degree ad1(a1, d1);
          atom_degree ad2(a2, d2);  // no direction in core

          size_t h = h1;
          if (h2 > h) h = h2;

          if ((l1 > 0 && m1 + k <= IN_INFO.atoms[a1].VALENCE &&
               m2 + k + 1 <= IN_INFO.atoms[a2].VALENCE) ||
              (l1 == 0 && k == rs.first.val &&
               m2 + k + 1 <= IN_INFO.atoms[a2].VALENCE)) {
            auto& set_left = map_T1[l1].at(rs);
            auto& set_right = map_T2_ind.at(rhps);

            rv_edge e_tmp;
            if (IN_INFO.atom_deg_map.at(ad1) >= IN_INFO.atom_deg_map.at(ad2)) {
              e_tmp = rv_edge("co", ad1, ad2, k);
            } else {
              e_tmp = rv_edge("co", ad2, ad1, k);
            }
            size_t e_tmp_ind = IN_INFO.rv_edge_map.at(e_tmp);

            size_t index_size = map_T2_ind.at(rhps).size();
            size_t start_index = (index_size - 1) * (l - 1) / len;
            size_t ind = start_index;
            bool flag = true;

            while (ind != start_index || flag) {
              flag = false;
              if (ind >= index_size) {
                ind -= index_size;
              }
              if (_stop) break;
              auto& _w2 = map_T2_ind[rhps][ind];
              // if (_DEBUG) cout << "ind = " <<ind << " index_size = " <<
              // index_size  << endl;
              ++ind;
              if (ind == index_size) {
                ind = 0;
              }

              auto w2 = traverse(_w2, all_trie_node_path, IN_INFO.rv1D_size);
              auto& _w2_second = map_T2.at(rhps).at(_w2);
              for (auto& _w1 : set_left) {
                if (globalTimeKeeper::tt.diff() - _start_time > time_limit) {
                  // time limit reached.
                  // cout << "len = " << l << ". Time over, limit is " <<
                  // time_limit << "s." << endl;
                  _stop = true;
                }
                if (UB_limit != 0 && vector_size > UB_limit) {
                  _stop = true;
                }
                if (_stop) break;

                auto w1 =
                    traverse(_w1.first, all_trie_node_path, IN_INFO.rv1D_size);
                auto new_rv = SUMinus3D(w1, w2, true);

                ++new_rv[e_tmp_ind];

                if (TestLessEqAB(
                        new_rv,
                        IN_INFO.rv1D)  // w1 + w2 + 1_mu + 1_gamma <= x_star
                ) {
                  root_status_T _rsT = make_pair(rs.first, rhps);
                  _rsT.second.val += k;
                  if (rs.second.height > _rsT.second.height) {
                    _rsT.second.height = rs.second.height;
                  }

                  // all_rs_T1[l].emplace(_rsT);
                  size_t ind = find_trie_index(new_rv, all_trie_node_path);

                  if (map_T1[l][_rsT].find(ind) != map_T1[l][_rsT].end()) {
                    if (map_T1[l][_rsT][ind]) {
                      _w1.second = true;
                      _w2_second = true;
                      DAG_node_T tmp(_rsT, l, ind);
                      map_DAG_T[tmp].emplace_back(rs, l1, _w1.first, rhps, _w2,
                                                  k);
                    }
                  }
                }  // if (TestLessEqAB)
              }    // for _w1
            }      // while
          }        // if
        }          // for k
      }            // for rs
    }              // for rhps

  }  // for l
}

void GenT2_0_DAG(
    const input_info& IN_INFO, unordered_set<root_status_height>& all_rs_T2,
    unordered_map<root_status_height, map<size_t, bool>>& map_T2,
    unordered_map<root_status_height, map<size_t, vector<DAG_node_height>>>&
        map_T2_W3,
    unordered_map<root_status_height, vector<size_t>>& map_T2_ind,
    vector<PathTrie>& all_trie_node_path,

    unordered_set<root_status_height>& all_rs_W2_core,
    unordered_map<root_status_height, map<size_t, bool>>& map_W2_core,
    vector<unordered_set<root_status>>& all_rs_W3,
    vector<unordered_map<root_status, map<size_t, bool>>>& map_W3,
    vector_2D<PathTrie>& all_trie_node,

    size_t chUB,

    bool exact = false

) {
  size_t max_h = 0;
  if (chUB >= IN_INFO.k_star) {
    max_h = chUB - IN_INFO.k_star;
  } else {
    max_h = 0;
  }
  size_t min_h = 0;
  size_t vector_size = 0;

  for (size_t h = min_h; h <= max_h; ++h) {
    if (h == 0) {
      for (auto& rsh : all_rs_W2_core) {
        if (rsh.height > chUB) {
          continue;
        }
        // cout << "rsh.color = "  << rsh.color << " rsh.deg = " << rsh.deg << "
        // rsh.val = " <<rsh.val << " rsh.height = " << rsh.height << endl;
        for (auto& _w : map_W2_core.at(rsh)) {
          auto w = traverse(_w.first, all_trie_node[0], IN_INFO.rv1D_size);

          if (exact && IN_INFO.rv1D != w) {
            continue;
          }

          all_rs_T2.insert(rsh);

          size_t ind = find_trie_index(w, all_trie_node_path);

          if (map_T2[rsh].find(ind) != map_T2[rsh].end()) {
            if (map_T2[rsh][ind]) {
              _w.second = true;
            }
          }
        }
      }
    } else {
      for (auto& rs : all_rs_W3[h]) {
        for (auto& _w : map_W3[h].at(rs)) {
          auto w = traverse(_w.first, all_trie_node[h], IN_INFO.rv1D_size);

          if (exact && IN_INFO.rv1D != w) {
            continue;
          }

          root_status_height rsh(rs, h + IN_INFO.k_star);

          all_rs_T2.emplace(rs, h + IN_INFO.k_star);

          size_t ind = find_trie_index(w, all_trie_node_path);

          if (map_T2[rsh].find(ind) != map_T2[rsh].end()) {
            if (map_T2[rsh][ind]) {
              _w.second = true;
            }
          }
        }
      }
    }
  }
}

#endif /* INCLUDES_FRINGE_TREE_HPP_ */
