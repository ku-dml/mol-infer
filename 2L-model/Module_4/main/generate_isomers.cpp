/* 2021/06/10
   Backtrack, isomer, for 2L-model, fc included

   Add fc, 06/09

   Change the way of storing paths to better check fc, 06/10

*/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

bool _DEBUG;  // a global boolean variable to show message for debugging
bool _EXP;

bool _LOOSE;

bool core_height_fixed = true;

#include "../include/chemical_graph.hpp"
#include "../include/cross_timer.h"
#include "../include/data_structures.hpp"
#include "../include/fringe_tree.hpp"
#include "../include/tools.hpp"
#include "../include/tree_signature.hpp"

using namespace std;

// initialize a global timekeeper
Timer globalTimeKeeper::tt;

int main(int argc, char** argv) {
  /////////////////////////////////
  // change this value  to true when debugging
  _DEBUG = false;
  _EXP = false;

  _LOOSE = true;  
  // when the variable "_LOOSE" is set to "false",
  // the code will only generate graphs that each component has the same descriptor (including fc) as the input graph
  ////////////////////////////////

  if (argc < 10) {
    cout << "Please supply: " << endl;
    cout << "1. an SDF file with one cyclic graph " << endl;
    cout << "2. an integer of time limit per iteration" << endl;
    cout << "3. an integer of vector size bound per iteration (0 if no limit)"
         << endl;
    cout << "4. an integer of number of sample graphs stored for each vector"
         << endl;
    cout << "5. an integer of global time limit" << endl;
    cout << "6. an integer of global number of path in DAG" << endl;
    cout << "7. an integer of number of wanted graphs" << endl;
    cout << "8. a file name of output" << endl;
    cout << "9. a file with information of partition" << endl;
    cout << "10. a file with fringe tree information" << endl;
    exit(-1);
  }
  globalTimeKeeper tk;

  string infile;
  try {
    infile = string(argv[1]);
  } catch (std::exception& e) {
    cout << "Caught an exception for the filename argument:" << endl;
    cout << e.what() << endl;
    cout << "The program will now exit." << endl;
    return 1;
  }

  time_limit = stoi(argv[2]);
  UB_limit = stoull(argv[3]);
  seq_limit = stoi(argv[4]);
  global_time_limit = stoi(argv[5]);
  global_path_limit = stoi(argv[6]);
  num_limit = stoi(argv[7]);

  _start_time = globalTimeKeeper::tt.diff();

  string outputfilename = argv[8];
  string partitionfilename;

  bool no_partition_info = false;
  // if (argc == 9) {
  //   no_partition_info = true;
  // } else {
    partitionfilename = argv[9];
  // }

  string csv_filename = argv[10];
  map <size_t, TreeSignature> TS_map;
  read_TS(csv_filename, TS_map);

  map <size_t, size_t> fc_map;
  fc_map.clear();
  for (auto& TS : TS_map){
    size_t id = fc_map.size();
    fc_map.emplace(TS.first, id);
  }

  Graph h = read_graph_sdf(infile);
  Graph g;
  H_suppressed_convert(h, g);
  vector<Vertex> core_set;
  vector<Vertex> internal_set;
  core_set.clear();
  // calcCoreVertexSet(g, core_set);
  calcInternalVertexSet(g, internal_set);
  // size_t g_n = calcEffectiveVertexNum(g);
  size_t g_n = 1000;
  size_t g_m = calcNumOfEdges(g);
  size_t core_height = 0;
  size_t bc = 0;

  double _ts0 = tk.tt.diff();

  vector<Vertex> base_vertices;
  vector_2D<Vertex> base_edges;
  vector<size_t> chLB_v;
  vector<size_t> chUB_v;
  vector<size_t> chLB_e;
  vector<size_t> chUB_e;

  vector<bool> fixed_v;
  vector<bool> fixed_e;

  vector_2D <size_t> fringe_tree_indices_v;
  vector_2D <size_t> fringe_tree_indices_e;

  if (no_partition_info) {
    get_partition(g, core_set, base_vertices, base_edges);
    fixed_v = vector<bool>(base_vertices.size(), false);
    fixed_e = vector<bool>(base_edges.size(), false);
    // get_partition(g, core_set, internal_set, base_vertices, base_edges);
  } else {
    get_partition(base_vertices, base_edges, chLB_v, chUB_v, chLB_e, chUB_e,
                  fixed_v, fixed_e, fringe_tree_indices_v, fringe_tree_indices_e, core_set, partitionfilename);
  }

  if (_DEBUG) {
    cout << "DEBUG on !!!" << endl;
  }
  if (_DEBUG) {
    cout << "|V_B| = " << base_vertices.size()
         << " |E_B| = " << base_edges.size() << endl;
  }

  cout.unsetf(ios::floatfield);
  cout.precision(3);

  vector <Atom> g_atom_set;
  g_atom_set.clear();
  get_atom_set(g_atom_set, g);

  vector <unsigned short> fc_all(TS_map.size(), 0);

  vector<component> set_tree_component;
  for (size_t i = 0; i < base_vertices.size(); ++i) {
    auto& v = base_vertices[i];
    component _cp;
    _cp.no_feasible_pair = false;
    _cp.no_feasible_pair_loose = false;
    _cp.with_core_height = false;
    // input_info IN_INFO;
    _cp.IN_INFO.num_fc = TS_map.size();
    _cp.IN_INFO.fc_map = fc_map;
    get_descriptors_tree(_cp.IN_INFO, g, core_set, internal_set, v, g_atom_set, TS_map);

    for (size_t i = 0; i < _cp.IN_INFO.num_fc; ++i){
      fc_all[i] += _cp.IN_INFO.rv.fc[i];
    }

    if (no_partition_info) {
      _cp.IN_INFO.chLB = 0;
      _cp.IN_INFO.chUB = _cp.IN_INFO.height;
      // _cp.IN_INFO.chUB = 100;
    } else {
      _cp.IN_INFO.chLB = chLB_v[i];
      _cp.IN_INFO.chUB = chUB_v[i];
      _cp.IN_INFO.fringe_tree_indices = fringe_tree_indices_v[i];
    }

    if (_cp.IN_INFO.height > core_height) {
      core_height = _cp.IN_INFO.height;
    }

    if (fixed_v[i]) {
      _cp.no_feasible_pair = true;
    }

    bc += _cp.IN_INFO.rv.bl;
    _cp.comp_graph.tree_component = true;

    _cp.IN_INFO.calculate_Gamma();

    set_tree_component.push_back(_cp);
    if (_DEBUG) cout << "v = " << v << endl;
  }

  vector<component> set_edge_component;
  for (size_t i = 0; i < base_edges.size(); ++i) {
    auto& e = base_edges[i];
    component _cp;
    _cp.no_feasible_pair = false;
    _cp.no_feasible_pair_loose = false;
    _cp.with_core_height = false;
    _cp.IN_INFO.num_fc = TS_map.size();
    _cp.IN_INFO.fc_map = fc_map;
    get_descriptors_edge(_cp.IN_INFO, g, core_set, internal_set, e, g_atom_set, TS_map);

    for (size_t i = 0; i < _cp.IN_INFO.num_fc; ++i){
      fc_all[i] += _cp.IN_INFO.rv.fc[i];
    }

    if (no_partition_info) {
      _cp.IN_INFO.chLB = 0;
      _cp.IN_INFO.chUB = _cp.IN_INFO.height;
      // _cp.IN_INFO.chUB = 100;
    } else {
      _cp.IN_INFO.chLB = chLB_e[i];
      _cp.IN_INFO.chUB = chUB_e[i];
      _cp.IN_INFO.fringe_tree_indices = fringe_tree_indices_e[i];
    }

    if (_cp.IN_INFO.height > core_height) {
      core_height = _cp.IN_INFO.height;
    }

    if (fixed_e[i]) {
      _cp.no_feasible_pair = true;
    }

    bc += _cp.IN_INFO.rv.bl;
    _cp.comp_graph.tree_component = false;

    _cp.IN_INFO.calculate_Gamma();

    set_edge_component.push_back(_cp);
    // if (_DEBUG) {
    //     cout << "e = ";
    //     for (auto& u : e){
    //         cout << u << " ";
    //     }
    //     cout<<endl;
    // }
  }
  for (size_t i = 0; i < base_vertices.size(); ++i){
    set_tree_component[i].IN_INFO.gen_rv1D_loose(fc_all);
  }
  for (size_t i = 0; i < base_edges.size(); ++i){
    set_edge_component[i].IN_INFO.gen_rv1D_loose(fc_all);
  }
  if (_EXP) {
    unordered_set<string> Lambda;
    unordered_set<string> Gamma_int;
    // unordered_set<string> Gamma_ex;
    // unordered_set<string> Gamma_co;

    for (auto& cp : set_tree_component) {
      auto& IN_INFO = cp.IN_INFO;
      for (auto& atom : IN_INFO.atoms) {
        Lambda.insert(atom.NAME);
      }
      size_t M = IN_INFO.num_kind_atom_degree;
      // for (size_t i = 0; i < M; ++i) {
      //   for (size_t j = 0; j <= i; ++j) {
      //     for (size_t k = 1; k <= 3; ++k) {
      //       if (IN_INFO.rv.ec_co[i][j][k] != 0) {
      //         string tmp = "";
      //         atom_degree ad1 = IN_INFO.atom_deg[i];
      //         tmp += IN_INFO.atoms[ad1.color].NAME + to_string(ad1.deg);
      //         atom_degree ad2 = IN_INFO.atom_deg[j];
      //         tmp += IN_INFO.atoms[ad2.color].NAME + to_string(ad2.deg);
      //         tmp += to_string(k);
      //         Gamma_co.insert(tmp);
      //       }
      //     }
      //   }
      // }
      for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j <= i; ++j) {
          for (size_t k = 1; k <= 3; ++k) {
            if (IN_INFO.rv.ec_int[i][j][k] != 0) {
              string tmp = "";
              atom_degree ad1 = IN_INFO.atom_deg[i];
              tmp += IN_INFO.atoms[ad1.color].NAME + to_string(ad1.deg);
              atom_degree ad2 = IN_INFO.atom_deg[j];
              tmp += IN_INFO.atoms[ad2.color].NAME + to_string(ad2.deg);
              tmp += to_string(k);
              Gamma_int.insert(tmp);
            }
          }
        }
      }
      // for (size_t i = 0; i < M; ++i) {
      //   for (size_t j = 0; j < M; ++j) {
      //     for (size_t k = 1; k <= 3; ++k) {
      //       if (IN_INFO.rv.ec_ex[i][j][k] != 0) {
      //         string tmp = "";
      //         atom_degree ad1 = IN_INFO.atom_deg[i];
      //         tmp += IN_INFO.atoms[ad1.color].NAME + to_string(ad1.deg);
      //         atom_degree ad2 = IN_INFO.atom_deg[j];
      //         tmp += IN_INFO.atoms[ad2.color].NAME + to_string(ad2.deg);
      //         tmp += to_string(k);
      //         Gamma_ex.insert(tmp);
      //       }
      //     }
      //   }
      // }
    }
    for (auto& cp : set_edge_component) {
      auto& IN_INFO = cp.IN_INFO;
      for (auto& atom : IN_INFO.atoms) {
        Lambda.insert(atom.NAME);
      }
      size_t M = IN_INFO.num_kind_atom_degree;
      // for (size_t i = 0; i < M; ++i) {
      //   for (size_t j = 0; j <= i; ++j) {
      //     for (size_t k = 1; k <= 3; ++k) {
      //       if (IN_INFO.rv.ec_co[i][j][k] != 0) {
      //         string tmp = "";
      //         atom_degree ad1 = IN_INFO.atom_deg[i];
      //         tmp += IN_INFO.atoms[ad1.color].NAME + to_string(ad1.deg);
      //         atom_degree ad2 = IN_INFO.atom_deg[j];
      //         tmp += IN_INFO.atoms[ad2.color].NAME + to_string(ad2.deg);
      //         tmp += to_string(k);
      //         Gamma_co.insert(tmp);
      //         // cout << tmp <<  endl;
      //       }
      //     }
      //   }
      // }
      for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j <= i; ++j) {
          for (size_t k = 1; k <= 3; ++k) {
            if (IN_INFO.rv.ec_int[i][j][k] != 0) {
              string tmp = "";
              atom_degree ad1 = IN_INFO.atom_deg[i];
              tmp += IN_INFO.atoms[ad1.color].NAME + to_string(ad1.deg);
              atom_degree ad2 = IN_INFO.atom_deg[j];
              tmp += IN_INFO.atoms[ad2.color].NAME + to_string(ad2.deg);
              tmp += to_string(k);
              Gamma_int.insert(tmp);
            }
          }
        }
      }
      // for (size_t i = 0; i < M; ++i) {
      //   for (size_t j = 0; j < M; ++j) {
      //     for (size_t k = 1; k <= 3; ++k) {
      //       if (IN_INFO.rv.ec_ex[i][j][k] != 0) {
      //         string tmp = "";
      //         atom_degree ad1 = IN_INFO.atom_deg[i];
      //         tmp += IN_INFO.atoms[ad1.color].NAME + to_string(ad1.deg);
      //         atom_degree ad2 = IN_INFO.atom_deg[j];
      //         tmp += IN_INFO.atoms[ad2.color].NAME + to_string(ad2.deg);
      //         tmp += to_string(k);
      //         Gamma_ex.insert(tmp);
      //       }
      //     }
      //   }
      // }
    }
	cout << "# vertices = " << calcEffectiveVertexNum(g) << endl;
    cout << "|V_B| = " << base_vertices.size() << endl;
    cout << "|E_B| = " << base_edges.size() << endl;
    cout << "|Gamma^int| = " << Gamma_int.size() << endl;
    // cout << "|Gamma^ex| = " << Gamma_ex.size() << endl;
    // cout << "|Gamma^co| = " << Gamma_co.size() << endl;
    cout << "|Lambda| = " << Lambda.size() << endl;
    cout << "bc = " << bc << endl;
    cout << "ch = " << core_height << endl;
  }

  sort(set_tree_component.begin(), set_tree_component.end(), cmp_v_component);
  sort(set_edge_component.begin(), set_edge_component.end(), cmp_e_component);

  size_t num_len = max(set_tree_component.size(), set_edge_component.size());
  size_t global_path_num = 0;

  for (size_t ind_comp = 0; ind_comp < num_len; ++ind_comp) {
    if (ind_comp < set_tree_component.size()) {
      auto& cp = set_tree_component[ind_comp];

      auto& IN_INFO = cp.IN_INFO;
      auto& comp_graph = cp.comp_graph;

      if (_MEMORY_OUT || _GLOBAL_TIME_OUT) {
        cp.no_feasible_pair = true;
        cp.no_feasible_pair_loose = true;
      }

      if (cp.no_feasible_pair) {
        if (IN_INFO.height == core_height) {
          cp.with_core_height = true;
        }
        if (cp.with_core_height) {
          comp_graph.num_h = 1;
          comp_graph.num = 0;
        } else {
          comp_graph.num = 1;
          comp_graph.num_h = 0;
        }
        // continue;
      } else {
        size_t ch = 0;
        if (IN_INFO.height >= IN_INFO.k_star) {
          ch = IN_INFO.height - IN_INFO.k_star;
        }
        size_t M = IN_INFO.num_kind_atoms;
        if (_DEBUG)
          cout << "ch = " << ch << " IN_INFO.height = " << IN_INFO.height
               << " v = " << IN_INFO.base_vertices[0]
               << " delta_i = " << IN_INFO.delta_1 << endl;
        // print_rv1D(IN_INFO, IN_INFO.rv1D);
        vector_2D<PathTrie> all_trie_node(ch + 1);
        for (size_t i = 0; i <= ch; ++i) {
          PathTrie tmp(0);
          all_trie_node[i].emplace_back(tmp);
        }

        //
        // W1: end fringe tree, always +1
        // W2: internal fringe tree,always +2
        // W3: fringe tree attached to the core vertex, +delta_i+1
        // W2_core: fringe tree attached to the core vertex, +delta_i
        vector<unordered_set<root_status>> all_rs_W1(ch);
        vector<unordered_map<root_status, map<size_t, bool>>> map_W1(ch);

        unordered_set<half_path_status> all_rs_W2;
        unordered_map<half_path_status, map<size_t, bool>> map_W2;
        unordered_map<half_path_status, vector<size_t>> map_W2_ind;

        vector<unordered_set<root_status>> all_rs_W3(ch + 1);
        vector<unordered_map<root_status, map<size_t, bool>>> map_W3(ch + 1);

        unordered_set<root_status_height> all_rs_W2_core;
        unordered_map<root_status_height, map<size_t, bool>> map_W2_core;

        // unordered_map <DAG_node, vector <DAG_edge>> map_DAG_W1;  //IDO
        unordered_map<DAG_node, vector<DAG_edge>> map_DAG_W3;

        vector<DAG_node> all_fp_W;
        vector<DAG_node> all_fp_W_loose;

        comp_graph.initialize(M + 1);

        // size_t chl1 = (ch - 1) / 2;
        size_t chl1 = ch - 1;  // IDO
        size_t chl2 = ch - 1 - chl1;

        if (ch > 0) {
          comp_graph.simple_tree = false;

          double time_tmp_1 = globalTimeKeeper::tt.diff();

          GenW1(IN_INFO, all_rs_W1, comp_graph.map_FT_W1, map_W1,
                comp_graph.All_FT_W1, all_trie_node, TS_map);
          GenW2_0(IN_INFO, all_rs_W2, comp_graph.map_FT_W2, map_W2, map_W2_ind,
                  comp_graph.All_FT_W2, all_trie_node, TS_map);
          GenW3_0(IN_INFO, all_rs_W3, comp_graph.map_FT_W3, map_W3,
                  comp_graph.All_FT_W3, all_trie_node, TS_map, IN_INFO.target_rs1.color,
                  IN_INFO.delta_1 + 1);

          double time_tmp_2 = globalTimeKeeper::tt.diff();

          // size_t chl1 = (ch - 1) / 2;
          merge_W1(IN_INFO, all_rs_W1, all_rs_W2, map_W1, map_W2, map_W2_ind,
                   all_trie_node, chl1);  // IDO

          double time_tmp_3 = globalTimeKeeper::tt.diff();

          merge_W3(IN_INFO, all_rs_W1, all_rs_W3, map_W1, map_W3, all_trie_node,
                   ch, ch, IN_INFO.delta_1);  // IDO
          // size_t chl2 = ch - 1 - chl1;
          // merge_W3_tree(IN_INFO, all_rs_W2, all_rs_W3, map_W2, map_W3,
          // all_trie_node, 0, chl2, IN_INFO.delta_1);

          double time_tmp_4 = globalTimeKeeper::tt.diff();

          fringe_tree_vectors_time += time_tmp_2 - time_tmp_1;
          end_subtree_vectors_time += time_tmp_3 - time_tmp_2;
          // rooted_core_subtrees_time += time_tmp_4 - time_tmp_3;
          end_subtree_vectors_time += time_tmp_4 - time_tmp_3;
        } else {
          comp_graph.simple_tree = true;

          double time_tmp_1 = globalTimeKeeper::tt.diff();
          
          GenW2_0_core(IN_INFO, all_rs_W2_core, comp_graph.map_FT_W2_core,
                       map_W2_core, comp_graph.All_FT_W2_core, all_trie_node,
                       TS_map, IN_INFO.delta_1);

          double time_tmp_2 = globalTimeKeeper::tt.diff();

          fringe_tree_vectors_time += time_tmp_2 - time_tmp_1;
        }

        //////////////////////////////////////////////////
        if (ch > 0) {
          double time_tmp_1 = globalTimeKeeper::tt.diff();

          vector<map_rv> all_seq_tmp;
          vector<map_rv> all_seq_loose_tmp;
          size_t num_tmp = 0;
          // combine_check_tree_component(IN_INFO, all_rs_W1, all_rs_W3, map_W1,
          // map_W3, all_trie_node,
          //     all_seq_tmp, chl1, chl2, IN_INFO.target_rs1, core_height,
          //     num_tmp);

          // combine_check_tree_component(IN_INFO, all_rs_W1, all_rs_W3, map_W1,
          // map_W3, all_trie_node,
          //     all_fp_W, chl1, chl2, IN_INFO.target_rs1, num_tmp);

          check_tree_component(IN_INFO, all_rs_W3, map_W3, all_trie_node,
                               all_fp_W, all_fp_W_loose, ch, IN_INFO.target_rs1);  // IDO

          // check all rs(a,d,m) of w3 and find same rs(a,d,m) of target rs1
          // (IN_INFO.target_rs1)

          double time_tmp_2 = globalTimeKeeper::tt.diff();

          // merge_W1_DAG(IN_INFO, all_rs_W1, all_rs_W2, map_W1, map_W2,
          // map_W2_ind, all_trie_node, map_DAG_W1, chl1);

          merge_W3_DAG(IN_INFO, all_rs_W1, all_rs_W3, map_W1, map_W3,
                       all_trie_node, map_DAG_W3, ch, ch,
                       IN_INFO.delta_1);  // IDO

          double time_tmp_3 = globalTimeKeeper::tt.diff();

          // merge_W3_tree_DAG(IN_INFO, all_rs_W2, all_rs_W3, map_W2, map_W3,
          // all_trie_node, map_DAG_W3, 0, chl2, IN_INFO.delta_1);

          merge_W1_DAG(IN_INFO, all_rs_W1, all_rs_W2, map_W1, map_W2,
                       map_W2_ind, all_trie_node, map_DAG_W3, chl1);

          double time_tmp_4 = globalTimeKeeper::tt.diff();

          feasible_pairs_time += time_tmp_2 - time_tmp_1;
          DAG_construct_time += time_tmp_3 - time_tmp_2;
          DAG_construct_time += time_tmp_4 - time_tmp_3;

          DAG_size_v += map_DAG_W3.size();  // IDO
          // for (auto& tmp : map_DAG_W1){   //IDO
          //     DAG_size_e += tmp.second.size();
          // }
          for (auto& tmp : map_DAG_W3) {
            DAG_size_e += tmp.second.size();
          }

          if (_DEBUG) {
            set<resource_vector_1D> xxx;
            cout << "DAG about v-component" << endl;
            int num = 1;
            for (auto& tmp : map_DAG_W3) {
              cout << num << ": ";
              auto w1 = traverse(tmp.first.ind, all_trie_node[tmp.first.l],
                                 IN_INFO.rv1D_size);
              xxx.insert(w1);
              for (auto& edge : tmp.second) {
                auto w2 = traverse(edge.node.ind, all_trie_node[edge.node.l],
                                   IN_INFO.rv1D_size);
                xxx.insert(w2);
                cout << "(" << tmp.first.rs.color << "," << tmp.first.rs.val
                     << "," << tmp.first.rs.deg << "," << tmp.first.ind << ","
                     << tmp.first.l << ") (";
                cout << edge.node.rs.color << "," << edge.node.rs.val << ","
                     << edge.node.rs.deg << "," << edge.node.ind << ","
                     << edge.node.l << ") (";
                cout << edge.w_ind << "," << edge.k << ")" << endl;
                cout << "head(" << tmp.first.ind << "," << tmp.first.l
                     << ") :" << endl;
                print_rv1D(IN_INFO, w1);
                print_rv1D_2(IN_INFO, w1);
                cout << "tail(" << edge.node.ind << "," << edge.node.l
                     << ") :" << endl;
                print_rv1D(IN_INFO, w2);
                print_rv1D_2(IN_INFO, w2);
                num += 1;
              }
            }
            cout << "kind of vector = " << xxx.size() << endl;
            // for (auto& w : xxx) {
            //   print_rv1D_2(IN_INFO, w);
            // }
            xxx.clear();
          }
		  if (_DEBUG) {
			  cout << "# feasible = " << feasible_pairs_sum << endl;
			  cout << "time for feasible = " << feasible_pairs_time << endl;
			  cout << "size of DAG, |V| = " << DAG_size_v << endl;
			  cout << "size of DAG, |E| = " << DAG_size_e << endl;
			  cout << "time to construct DAG = " << DAG_construct_time << endl;
		  }

          // cout << "all_fp_W.size() = " << all_fp_W.size() << endl;

          // get all paths by traversing DAG
          double time_tmp_5_1 = globalTimeKeeper::tt.diff();

          size_t path_num = 0;

          for (auto& fp : all_fp_W) {
            if (_GLOBAL_TIME_OUT) continue;

            // vector <map_rv> all_seq_W1;
            vector<map_rv> all_seq_W3;

            // all_seq_W1.clear();
            // all_seq_W3.clear();
            size_t all_seq_size = 0;

            map_rv seq;
            time_count_path = globalTimeKeeper::tt.diff();
            Gen_Path_W_in_T_path_count(
                map_DAG_W3, comp_graph.map_FT_W1, comp_graph.map_FT_W2,
                comp_graph.map_FT_W3, fp, seq, all_seq_size);
            // Gen_Path_W_path(map_DAG_W1, comp_graph.map_FT_W2,
            // comp_graph.map_FT_W1, fp.node2, seq, all_seq_W1);
            path_num += all_seq_size;  // all_seq_W1.size() * all_seq_W3.size();
          }

          double time_tmp_5 = globalTimeKeeper::tt.diff();
          DAG_traverse_path_time += time_tmp_5 - time_tmp_5_1;

         // cout << "# paths = " << path_num << endl;
          //cout << "time to traverse all paths from DAG = "
           //    << DAG_traverse_path_time << endl;

          // get numbers by traversing DAG
          num_tmp = 0;

          for (auto& fp : all_fp_W) {
            if (_GLOBAL_TIME_OUT) continue;

            size_t num_tmp_W3 = 0;
            // size_t num_tmp_W1 = 0;

            Gen_Path_W_in_T_num(map_DAG_W3, comp_graph.map_FT_W1,
                                comp_graph.map_FT_W2, comp_graph.map_FT_W3, fp,
                                1, num_tmp_W3);
            // Gen_Path_W_num(map_DAG_W1, comp_graph.map_FT_W2,
            // comp_graph.map_FT_W1, fp.node2, 1, num_tmp_W1);
            num_tmp += num_tmp_W3;  // num_tmp_W3 * num_tmp_W1;
          }

          double time_tmp_5_2 = globalTimeKeeper::tt.diff();
          DAG_traverse_num_time += time_tmp_5_2 - time_tmp_5;
		  if (_DEBUG) {
			  cout << "time to get exact number of graphs from DAG = "
				  << DAG_traverse_num_time << endl;
			  cout << "# graphs (v-comp) = " << num_tmp << endl;
		  }

          // get all graphs by traversing DAG
          for (auto& fp : all_fp_W) {
            if (seq_limit != 0 && all_seq_tmp.size() >= seq_limit) {
              break;
            }

            // vector <map_rv> all_seq_W1;
            vector<map_rv> all_seq_W3;

            // all_seq_W1.clear();
            all_seq_W3.clear();

            map_rv seq;
            Gen_Path_W_in_T(map_DAG_W3, comp_graph.map_FT_W1,
                            comp_graph.map_FT_W2, comp_graph.map_FT_W3, fp, seq,
                            all_seq_W3, global_path_num);
            // Gen_Path_W(map_DAG_W1, comp_graph.map_FT_W2,
            // comp_graph.map_FT_W1, fp.node2, seq, all_seq_W1);

            // for (auto& seq1 : all_seq_W3){
            //     for (auto& seq2 : all_seq_W1){
            //         map_rv map_rv_tmp(seq1, false);
            //         all_seq_tmp.emplace_back(map_rv_tmp, seq2, fp.k, true);
            //         if (seq_limit != 0 && all_seq_tmp.size() >= seq_limit){
            //             break;
            //         }
            //     }
            //     if (seq_limit != 0 && all_seq_tmp.size() >= seq_limit){
            //         break;
            //     }
            // }

            all_seq_tmp.insert(all_seq_tmp.end(), all_seq_W3.begin(),
                               all_seq_W3.end());
          }

          if (_LOOSE){
            for (auto& fp : all_fp_W_loose){
              if (seq_limit != 0 && all_seq_loose_tmp.size() >= seq_limit) {
                break;
              }
              vector<map_rv> all_seq_W3;
              all_seq_W3.clear();

              map_rv seq;
              Gen_Path_W_in_T(map_DAG_W3, comp_graph.map_FT_W1,
                              comp_graph.map_FT_W2, comp_graph.map_FT_W3, fp, seq,
                              all_seq_W3, global_path_num);
              all_seq_loose_tmp.insert(all_seq_loose_tmp.end(), all_seq_W3.begin(),
                               all_seq_W3.end());
            }
          }

          if (IN_INFO.height == core_height) {
            cp.with_core_height = true;
            comp_graph.num_h = num_tmp;
            if (all_seq_tmp.size() > 0) {
              for (auto& w : all_seq_tmp) {
                comp_graph.all_seq_h.emplace_back(w, IN_INFO.height);
              }

              if (comp_graph.num_h < comp_graph.all_seq_h.size()){
                comp_graph.num_h = comp_graph.all_seq_h.size();
              }
        
            } else {
              cp.no_feasible_pair = true;
              comp_graph.num_h = 1;
            }
            comp_graph.num = 0;
            // cout << "comp_graph.all_seq_h.size() = " <<
            // comp_graph.all_seq_h.size() << endl;

            if (_LOOSE){
              if (all_seq_loose_tmp.size() > 0) {
                for (auto& w : all_seq_loose_tmp) {
                  comp_graph.all_seq_h_loose.emplace_back(w, IN_INFO.height);
                }
              } else if (all_seq_loose_tmp.size() + all_seq_tmp.size() == 0){ // change, 0610
                cp.no_feasible_pair_loose = true;
              }
            }
          } else {
            comp_graph.num = num_tmp;
            if (all_seq_tmp.size() > 0) {
              for (auto& w : all_seq_tmp) {
                comp_graph.all_seq.emplace_back(w, IN_INFO.height);
              }

              if (comp_graph.num < comp_graph.all_seq.size()){
                comp_graph.num = comp_graph.all_seq.size();
              }

            } else {
              cp.no_feasible_pair = true;
              comp_graph.num = 1;
            }
            comp_graph.num_h = 0;
            // cout << "comp_graph.all_seq.size() = " <<
            // comp_graph.all_seq.size() << endl;

            if (_LOOSE){
              if (all_seq_loose_tmp.size() > 0) {
                for (auto& w : all_seq_loose_tmp) {
                  comp_graph.all_seq_loose.emplace_back(w, IN_INFO.height);
                }
              } else if (all_seq_loose_tmp.size() + all_seq_tmp.size() == 0){ // change, 0610
                cp.no_feasible_pair_loose = true;
              }
            }
          }

          double time_tmp_6 = globalTimeKeeper::tt.diff();

          DAG_traverse_graph_time += time_tmp_6 - time_tmp_5_2;

         // cout << "time to traverse to get all graphs from DAG = "
          //     << DAG_traverse_graph_time << endl;

        } else {
          vector<PathTrie> all_trie_node_path;
          PathTrie tmp(0);
          all_trie_node_path.emplace_back(tmp);
          unordered_set<root_status_height> all_rs_T2;
          unordered_map<root_status_height, map<size_t, bool>> map_T2;
          unordered_map<root_status_height,
                        map<size_t, vector<DAG_node_height>>>
              map_T2_W3;
          unordered_map<root_status_height, vector<size_t>> map_T2_ind;

          double time_tmp_1 = globalTimeKeeper::tt.diff();

          GenT2_0(IN_INFO, all_rs_T2, map_T2, map_T2_W3, map_T2_ind,
                  all_trie_node_path, all_rs_W2_core, map_W2_core, all_rs_W3,
                  map_W3, all_trie_node, IN_INFO.height, true);

          double time_tmp_2 = globalTimeKeeper::tt.diff();

          // rooted_core_subtrees_time += time_tmp_2 - time_tmp_1;

          for (auto& rsh : all_rs_T2) {
            for (auto& w : map_T2_W3.at(rsh)) {
              if (rsh.height == core_height) {
                for (auto& tmp_node : w.second) {
                  auto tmp_pair = make_pair(rsh, tmp_node.ind);
                  for (auto& tmp : comp_graph.map_FT_W2_core.at(tmp_pair)) {
                    auto& FT = comp_graph.All_FT_W2_core[rsh.color][tmp];
                    if (FT.height == rsh.height) {
                      map_rv rv_tmp(rsh, tmp);
                      comp_graph.all_seq_h.emplace_back(rv_tmp, rsh.height);
                      ++comp_graph.num_h;
                    }
                  }
                }
              } else {
                for (auto& tmp_node : w.second) {
                  auto tmp_pair = make_pair(rsh, tmp_node.ind);
                  for (auto& tmp : comp_graph.map_FT_W2_core.at(tmp_pair)) {
                    auto& FT = comp_graph.All_FT_W2_core[rsh.color][tmp];
                    if (FT.height == rsh.height) {
                      map_rv rv_tmp(rsh, tmp);
                      comp_graph.all_seq.emplace_back(rv_tmp, rsh.height);
                      ++comp_graph.num;
                    }
                  }
                }
              }
            }
          }

          if (comp_graph.all_seq.size() == 0 && comp_graph.all_seq_h.size() == 0){
            cp.no_feasible_pair = true;
          }

          if (_LOOSE){
            all_rs_T2.clear();
            map_T2.clear();
            map_T2_W3.clear();
            map_T2_ind.clear();
            GenT2_0(IN_INFO, all_rs_T2, map_T2, map_T2_W3, map_T2_ind,
                  all_trie_node_path, all_rs_W2_core, map_W2_core, all_rs_W3,
                  map_W3, all_trie_node, IN_INFO.height, true, true);

            for (auto& rsh : all_rs_T2) {
              for (auto& w : map_T2_W3.at(rsh)) {
                if (rsh.height == core_height) {
                  for (auto& tmp_node : w.second) {
                    auto tmp_pair = make_pair(rsh, tmp_node.ind);
                    for (auto& tmp : comp_graph.map_FT_W2_core.at(tmp_pair)) {
                      auto& FT = comp_graph.All_FT_W2_core[rsh.color][tmp];
                      if (FT.height == rsh.height) {
                        map_rv rv_tmp(rsh, tmp);
                        comp_graph.all_seq_h_loose.emplace_back(rv_tmp, rsh.height);
                      }
                    }
                  }
                } else {
                  for (auto& tmp_node : w.second) {
                    auto tmp_pair = make_pair(rsh, tmp_node.ind);
                    for (auto& tmp : comp_graph.map_FT_W2_core.at(tmp_pair)) {
                      auto& FT = comp_graph.All_FT_W2_core[rsh.color][tmp];
                      if (FT.height == rsh.height) {
                        map_rv rv_tmp(rsh, tmp);
                        comp_graph.all_seq_loose.emplace_back(rv_tmp, rsh.height);
                      }
                    }
                  }
                }
              }
            }
            // change, 0610
            if (comp_graph.all_seq_loose.size() + comp_graph.all_seq.size() == 0 
              && 
              comp_graph.all_seq_h_loose.size() + comp_graph.all_seq_h.size() == 0
            ){
              cp.no_feasible_pair_loose = true;
            }
          }

        }
      }
    }

    if (ind_comp < set_edge_component.size()) {
      auto& cp = set_edge_component[ind_comp];

      auto& IN_INFO = cp.IN_INFO;
      auto& comp_graph = cp.comp_graph;

      if (_MEMORY_OUT || _GLOBAL_TIME_OUT) {
        cp.no_feasible_pair = true;
        cp.no_feasible_pair_loose = true;
      }

      if (cp.no_feasible_pair) {
        if (IN_INFO.height == core_height) {
          cp.with_core_height = true;
        }
        if (cp.with_core_height) {
          comp_graph.num_h = 1;
          comp_graph.num = 0;
        } else {
          comp_graph.num = 1;
          comp_graph.num_h = 0;
        }
        // continue;
      } else {
        size_t chlb = 0;
        if (IN_INFO.chLB >= IN_INFO.k_star) {
          chlb = IN_INFO.chLB - IN_INFO.k_star;
        }
        size_t chub = 0;
        if (IN_INFO.chUB >= IN_INFO.k_star) {
          chub = IN_INFO.chUB - IN_INFO.k_star;
        }

        size_t M = IN_INFO.num_kind_atoms;
        size_t len = IN_INFO.base_vertices.size() - 1;
        size_t l2 = (len - 1) / 2;
        size_t l1 = len - 1 - l2;

        if (_DEBUG)
          cout << "IN_INFO.base_vertices[0] = " << IN_INFO.base_vertices[0]
               << " IN_INFO.base_vertices[IN_INFO.base_vertices.size() - 1] = "
               << IN_INFO.base_vertices[IN_INFO.base_vertices.size() - 1]
               << endl;
        if (_DEBUG)
          cout << "chlb = " << chlb << " IN_INFO.chLB = " << IN_INFO.chLB
               << " chub = " << chub << " IN_INFO.chUB = " << IN_INFO.chUB
               << endl;
        if (_DEBUG) cout << "l1 = " << l1 << " l2 = " << l2 << endl;
        // print_rv1D(IN_INFO, IN_INFO.rv1D);
        vector_2D<PathTrie> all_trie_node(chub + 1);
        for (size_t i = 0; i <= chub; ++i) {
          PathTrie tmp(0);
          all_trie_node[i].emplace_back(tmp);
        }

        //
        // W1: end fringe tree, always +1
        // W2: internal fringe tree,always +2
        // W3: fringe tree attached to the core vertex, +delta_i+1
        // W2_core: fringe tree attached to the core vertex, +delta_i
        vector<unordered_set<root_status>> all_rs_W1(chub);
        vector<unordered_map<root_status, map<size_t, bool>>> map_W1(chub);

        unordered_set<half_path_status> all_rs_W2;
        unordered_map<half_path_status, map<size_t, bool>> map_W2;
        unordered_map<half_path_status, vector<size_t>> map_W2_ind;

        vector<unordered_set<root_status>> all_rs_W3(chub + 1);
        vector<unordered_map<root_status, map<size_t, bool>>> map_W3(chub + 1);

        unordered_set<root_status_height> all_rs_W2_core;
        unordered_map<root_status_height, map<size_t, bool>> map_W2_core;

        unordered_map<DAG_node, vector<DAG_edge>> map_DAG_W;
        unordered_map<DAG_node_T, vector<DAG_edge_T>> map_DAG_T;

        vector<feasible_pair_T> all_fp_T;

        comp_graph.initialize(M + 1);

        if (chub > 0) {
          double time_tmp_1 = globalTimeKeeper::tt.diff();

          GenW1(IN_INFO, all_rs_W1, comp_graph.map_FT_W1, map_W1,
                comp_graph.All_FT_W1, all_trie_node, TS_map);
          GenW2_0(IN_INFO, all_rs_W2, comp_graph.map_FT_W2, map_W2, map_W2_ind,
                  comp_graph.All_FT_W2, all_trie_node, TS_map);
          GenW3_0(IN_INFO, all_rs_W3, comp_graph.map_FT_W3, map_W3,
                  comp_graph.All_FT_W3, all_trie_node, TS_map, M, 3);

          double time_tmp_2 = globalTimeKeeper::tt.diff();

          merge_W1(IN_INFO, all_rs_W1, all_rs_W2, map_W1, map_W2, map_W2_ind,
                   all_trie_node, chub - 1);

          double time_tmp_3 = globalTimeKeeper::tt.diff();

          merge_W3(IN_INFO, all_rs_W1, all_rs_W3, map_W1, map_W3, all_trie_node,
                   1, chub, 2);

          double time_tmp_4 = globalTimeKeeper::tt.diff();

          fringe_tree_vectors_time += time_tmp_2 - time_tmp_1;
          end_subtree_vectors_time += time_tmp_3 - time_tmp_2;
          rooted_core_subtrees_time += time_tmp_4 - time_tmp_3;
        }

        double time_tmp_1 = globalTimeKeeper::tt.diff();

        GenW2_0_core(IN_INFO, all_rs_W2_core, comp_graph.map_FT_W2_core,
                     map_W2_core, comp_graph.All_FT_W2_core, all_trie_node, TS_map, 2);

        double time_tmp_2 = globalTimeKeeper::tt.diff();
        fringe_tree_vectors_time += time_tmp_2 - time_tmp_1;

        vector<PathTrie> all_trie_node_path;
        PathTrie tmp(0);
        all_trie_node_path.emplace_back(tmp);

        vector<unordered_set<root_status_T>> all_rs_T1(l1 + 1);
        vector<unordered_map<root_status_T, map<size_t, bool>>> map_T1(l1 + 1);

        GenT1(IN_INFO, all_rs_T1, map_T1, all_trie_node_path,
              IN_INFO.target_rs1);
        GenT1(IN_INFO, all_rs_T1, map_T1, all_trie_node_path,
              IN_INFO.target_rs2);

        unordered_set<root_status_height> all_rs_T2;
        unordered_map<root_status_height, map<size_t, bool>> map_T2;
        unordered_map<root_status_height, map<size_t, vector<DAG_node_height>>>
            map_T2_W3;
        unordered_map<root_status_height, vector<size_t>> map_T2_ind;

        time_tmp_1 = globalTimeKeeper::tt.diff();

        GenT2_0(IN_INFO, all_rs_T2, map_T2, map_T2_W3, map_T2_ind,
                all_trie_node_path, all_rs_W2_core, map_W2_core, all_rs_W3,
                map_W3, all_trie_node, IN_INFO.chUB, false);

        time_tmp_2 = globalTimeKeeper::tt.diff();
        rooted_core_subtrees_time += time_tmp_2 - time_tmp_1;

        merge_T1(IN_INFO, all_rs_T1, all_rs_T2, map_T1, map_T2, map_T2_ind,
                 all_trie_node_path, l1);

        double time_tmp_3 = globalTimeKeeper::tt.diff();
        bi_rooted_core_subtrees_time += time_tmp_3 - time_tmp_2;

        combine_check_edge_component(
            IN_INFO, all_rs_T1, map_T1, all_trie_node_path, all_fp_T, l1, l2,
            IN_INFO.target_rs1, IN_INFO.target_rs2, core_height,
            comp_graph.num_h, comp_graph.num);

        double time_tmp_4 = globalTimeKeeper::tt.diff();
        feasible_pairs_time += time_tmp_4 - time_tmp_3;

        merge_T1_DAG(IN_INFO, all_rs_T1, all_rs_T2, map_T1, map_T2, map_T2_W3,
                     map_T2_ind, all_trie_node_path, map_DAG_T, l1);

        GenT2_0_DAG(IN_INFO, all_rs_T2, map_T2, map_T2_W3, map_T2_ind,
                    all_trie_node_path, all_rs_W2_core, map_W2_core, all_rs_W3,
                    map_W3, all_trie_node, IN_INFO.chUB, false);

        double time_tmp_4_2 = globalTimeKeeper::tt.diff();

        if (chub > 0) {
          double time_tmp_5 = globalTimeKeeper::tt.diff();

          merge_W3_DAG(IN_INFO, all_rs_W1, all_rs_W3, map_W1, map_W3,
                       all_trie_node, map_DAG_W, 1, chub, 2);

          double time_tmp_6 = globalTimeKeeper::tt.diff();

          merge_W1_DAG(IN_INFO, all_rs_W1, all_rs_W2, map_W1, map_W2,
                       map_W2_ind, all_trie_node, map_DAG_W, chub - 1);

          double time_tmp_7 = globalTimeKeeper::tt.diff();

          end_subtree_vectors_time += time_tmp_7 - time_tmp_6;
          rooted_core_subtrees_time += time_tmp_6 - time_tmp_5;
        }

        double time_tmp_8 = globalTimeKeeper::tt.diff();

        DAG_construct_time += time_tmp_8 - time_tmp_4;
        DAG_size_v += map_DAG_W.size() + map_DAG_T.size();
        size_t DAG_size_e_W = 0;
        size_t DAG_size_e_T = 0;

        for (auto& tmp : map_DAG_W) {
          DAG_size_e_W += tmp.second.size();
        }
        for (auto& tmp : map_DAG_T) {
          DAG_size_e_T += tmp.second.size();
        }
        DAG_size_e += DAG_size_e_T + DAG_size_e_W;

        // DAG about v-component in e-component information

        if (_DEBUG) {
          set<resource_vector_1D> xxx;
          cout << "DAG about v-component in e-component" << endl;
          int num = 1;
          for (auto& tmp : map_DAG_W) {
            cout << num << ": ";
            auto w1 = traverse(tmp.first.ind, all_trie_node[tmp.first.l],
                               IN_INFO.rv1D_size);
            xxx.insert(w1);
            for (auto& edge : tmp.second) {
              auto w2 = traverse(edge.node.ind, all_trie_node[edge.node.l],
                                 IN_INFO.rv1D_size);
              xxx.insert(w2);
              cout << "(" << tmp.first.rs.color << "," << tmp.first.rs.val
                   << "," << tmp.first.rs.deg << "," << tmp.first.ind << ","
                   << tmp.first.l << ") (";
              cout << edge.node.rs.color << "," << edge.node.rs.val << ","
                   << edge.node.rs.deg << "," << edge.node.ind << ","
                   << edge.node.l << ") (";
              cout << edge.w_ind << "," << edge.k << ")" << endl;
              cout << "head(" << tmp.first.ind << "," << tmp.first.l
                   << ") :" << endl;
              print_rv1D(IN_INFO, w1);
              print_rv1D_2(IN_INFO, w1);
              cout << "tail(" << edge.node.ind << "," << edge.node.l
                   << ") :" << endl;
              print_rv1D(IN_INFO, w2);
              print_rv1D_2(IN_INFO, w2);
              num += 1;
            }
          }
          cout << "kind of vector = " << xxx.size() << endl;
          xxx.clear();
        }

        // DAG about e-component information

        if (_DEBUG) {
          set<resource_vector_1D> xxx;
          cout << "DAG about e-component" << endl;
          int num = 1;
          for (auto& tmp : map_DAG_T) {  // output everthing
            auto w1 =
                traverse(tmp.first.ind, all_trie_node_path, IN_INFO.rv1D_size);
            xxx.insert(w1);
            for (auto& edge : tmp.second) {
              auto w2 = traverse(edge.node.ind, all_trie_node_path,
                                 IN_INFO.rv1D_size);
              auto w3 =
                  traverse(edge.T2_ind, all_trie_node_path, IN_INFO.rv1D_size);
              xxx.insert(w2);
              cout << num << ": ";
              cout << "(" << tmp.first.rsT.first.color << ","
                   << tmp.first.rsT.first.val << "," << tmp.first.rsT.first.deg
                   << "," << tmp.first.rsT.second.color << ","
                   << tmp.first.rsT.second.val << ","
                   << tmp.first.rsT.second.deg << ","
                   << tmp.first.rsT.second.height << "," << tmp.first.ind << ","
                   << tmp.first.l << ") (";
              cout << edge.node.rsT.first.color << ","
                   << edge.node.rsT.first.val << "," << edge.node.rsT.first.deg
                   << "," << edge.node.rsT.second.color << ","
                   << edge.node.rsT.second.val << ","
                   << edge.node.rsT.second.deg << ","
                   << edge.node.rsT.second.height << "," << edge.node.ind << ","
                   << edge.node.l << ") (";
              cout << edge.T2_ind << "," << edge.k << "," << edge.rsh.color
                   << "," << edge.rsh.val << "," << edge.rsh.deg << ","
                   << edge.rsh.height << ")" << endl;
              cout << "head(" << tmp.first.ind << "," << tmp.first.l
                   << ") :" << endl;
              print_rv1D(IN_INFO, w1);
              print_rv1D_2(IN_INFO, w1);
              cout << "tail(" << edge.node.ind << "," << edge.node.l
                   << ") :" << endl;
              print_rv1D(IN_INFO, w2);
              print_rv1D_2(IN_INFO, w2);
              cout << "arc:" << endl;
              print_rv1D(IN_INFO, w3);
              print_rv1D_2(IN_INFO, w3);
              num += 1;
            }
          }
          cout << "kind of vector = " << xxx.size() << endl;

          // for (auto& w : xxx) {
          //   print_rv1D_2(IN_INFO, w);
          // }
          xxx.clear();
        }

       /* cout << "# feasible pairs = " << feasible_pairs_sum << endl;
        cout << "time for feasible pairs = " << feasible_pairs_time << endl;
        cout << "size of DAG of core part, |V| = " << map_DAG_T.size() << endl;
        cout << "size of DAG of core part, |E| = " << DAG_size_e_T << endl;
        cout << "size of DAG of non-core part, |V| = " << map_DAG_W.size()
             << endl;
        cout << "size of DAG of non-core part, |E| = " << DAG_size_e_W << endl;
        cout << "time to construct DAG of core part = "
             << time_tmp_4_2 - time_tmp_4 << endl;
        cout << "time to construct DAG of non-core part = "
             << time_tmp_8 - time_tmp_4_2 << endl;*/

        // get all paths from DAG_T

        size_t path_num_core = 0;

        for (auto& fp : all_fp_T) {
          if (_GLOBAL_TIME_OUT) continue;
          if (fp.dir != DAG_direction_status::EXACT_possible) continue;
          vector<map_tree> all_seq_T1;
          vector<map_tree> all_seq_T2;

          all_seq_T1.clear();
          all_seq_T2.clear();

          map_tree seq;

          Gen_Path_T_only_path(map_DAG_T, map_DAG_W, map_T2_W3,
                               comp_graph.map_FT_W1, comp_graph.map_FT_W2,
                               comp_graph.map_FT_W3, comp_graph.map_FT_W2_core,
                               fp.node1, seq, all_seq_T1);
          Gen_Path_T_only_path(map_DAG_T, map_DAG_W, map_T2_W3,
                               comp_graph.map_FT_W1, comp_graph.map_FT_W2,
                               comp_graph.map_FT_W3, comp_graph.map_FT_W2_core,
                               fp.node2, seq, all_seq_T2);
          path_num_core += all_seq_T1.size() * all_seq_T2.size();
        }

        double time_tmp_9_3 = globalTimeKeeper::tt.diff();
        DAG_traverse_path_time += time_tmp_9_3 - time_tmp_8;

        //cout << "# paths of DAG of core part = " << path_num_core << endl;
        //cout << "time to traverse all paths from DAG of only core part = "
         //    << time_tmp_9_3 - time_tmp_8 << endl;

        // get all paths from DAG

        size_t path_num = 0;

        for (auto& fp : all_fp_T) {
          if (_GLOBAL_TIME_OUT) continue;
          if (fp.dir != DAG_direction_status::EXACT_possible) continue;
          vector<map_tree> all_seq_T1;
          vector<map_tree> all_seq_T2;

          all_seq_T1.clear();
          all_seq_T2.clear();

          map_tree seq;

          Gen_Path_T_path(map_DAG_T, map_DAG_W, map_T2_W3, comp_graph.map_FT_W1,
                          comp_graph.map_FT_W2, comp_graph.map_FT_W3,
                          comp_graph.map_FT_W2_core, fp.node1, seq, all_seq_T1);
          Gen_Path_T_path(map_DAG_T, map_DAG_W, map_T2_W3, comp_graph.map_FT_W1,
                          comp_graph.map_FT_W2, comp_graph.map_FT_W3,
                          comp_graph.map_FT_W2_core, fp.node2, seq, all_seq_T2);
          path_num += all_seq_T1.size() * all_seq_T2.size();
        }

        double time_tmp_9_2 = globalTimeKeeper::tt.diff();
        DAG_traverse_path_time += time_tmp_9_2 - time_tmp_9_3;

        //cout << "# paths of DAG of both core and non-core part = " << path_num
         //    << endl;
        //cout << "time to traverse all paths from DAG of both core and non-core "
         //       "part = "
          //   << time_tmp_9_2 - time_tmp_9_3 << endl;

        // get number from DAG

        for (auto& fp : all_fp_T) {
          if (fp.dir != DAG_direction_status::EXACT_possible) continue;
          map_tree seq;

          size_t num_tmp_T1 = 0;
          size_t num_tmp_T2 = 0;

          Gen_Path_T_num(map_DAG_T, map_DAG_W, map_T2_W3, comp_graph.map_FT_W1,
                         comp_graph.map_FT_W2, comp_graph.map_FT_W3,
                         comp_graph.map_FT_W2_core, fp.node1, 1, num_tmp_T1);
          Gen_Path_T_num(map_DAG_T, map_DAG_W, map_T2_W3, comp_graph.map_FT_W1,
                         comp_graph.map_FT_W2, comp_graph.map_FT_W3,
                         comp_graph.map_FT_W2_core, fp.node2, 1, num_tmp_T2);

          if (max(fp.node1.rsT.second.height, fp.node2.rsT.second.height) ==
              core_height) {
            comp_graph.num_h += num_tmp_T1 * num_tmp_T2;
          } else {
            comp_graph.num += num_tmp_T1 * num_tmp_T2;
          }
        }

        double time_tmp_9_4 = globalTimeKeeper::tt.diff();
        DAG_traverse_num_time += time_tmp_9_4 - time_tmp_9_2;

      //  cout << "# graphs (e-comp) = " << comp_graph.num + comp_graph.num_h
        //     << endl;
        //cout << "time to get exact number of graphs from DAG = "
          //   << DAG_traverse_num_time << endl;
        // cout << "all_fp_T.size() = " << all_fp_T.size() << endl;

        for (auto& fp : all_fp_T) {
          if (seq_limit != 0 && comp_graph.all_seq.size() >= seq_limit &&
              comp_graph.all_seq_h.size() >= seq_limit) {
            break;
          }

          if (!_LOOSE && fp.dir != DAG_direction_status::EXACT_possible) continue;

          vector<map_tree> all_seq_T1;
          vector<map_tree> all_seq_T2;

          all_seq_T1.clear();
          all_seq_T2.clear();

          map_tree seq;

          Gen_Path_T(map_DAG_T, map_DAG_W, map_T2_W3, comp_graph.map_FT_W1,
                     comp_graph.map_FT_W2, comp_graph.map_FT_W3,
                     comp_graph.map_FT_W2_core, fp.node1, seq, all_seq_T1, 
                     global_path_num);
          Gen_Path_T(map_DAG_T, map_DAG_W, map_T2_W3, comp_graph.map_FT_W1,
                     comp_graph.map_FT_W2, comp_graph.map_FT_W3,
                     comp_graph.map_FT_W2_core, fp.node2, seq, all_seq_T2, 
                     global_path_num);

          size_t h =
              max(fp.node1.rsT.second.height, fp.node2.rsT.second.height);

          for (auto& seq1 : all_seq_T1) {
            for (auto& seq2 : all_seq_T2) {
              if (h == core_height) {
                // if (fp.dir == DAG_direction_status::EXACT_possible)
                  comp_graph.all_seq_h.emplace_back(
                      seq1, seq2, fp.k, max(seq1.height, seq2.height), false);
                // if (_LOOSE)
                //   comp_graph.all_seq_h_loose.emplace_back(
                //       seq1, seq2, fp.k, max(seq1.height, seq2.height), false);
                // cout << "seq1.trees.size() = " << seq1.trees.size() << "
                // seq1.mul.size() = " << seq1.mul.size() << endl; seq1.print();
                // cout << "seq2.trees.size() = " << seq2.trees.size() << "
                // seq2.mul.size() = " << seq2.mul.size() << endl; seq2.print();
                if (seq_limit != 0 &&
                    comp_graph.all_seq_h.size() >= seq_limit) {
                  break;
                }
              } else {
                // if (fp.dir == DAG_direction_status::EXACT_possible)
                  comp_graph.all_seq.emplace_back(
                      seq1, seq2, fp.k, max(seq1.height, seq2.height), false);
                // if (_LOOSE)
                //   comp_graph.all_seq_loose.emplace_back(
                //       seq1, seq2, fp.k, max(seq1.height, seq2.height), false);

                // cout << "seq1.trees.size() = " << seq1.trees.size() << "
                // seq1.mul.size() = " << seq1.mul.size() << endl; seq1.print();
                // cout << "seq2.trees.size() = " << seq2.trees.size() << "
                // seq2.mul.size() = " << seq2.mul.size() << endl; seq2.print();
                if (seq_limit != 0 && comp_graph.all_seq.size() >= seq_limit) {
                  break;
                }
              }
            }
            if (seq_limit != 0 && comp_graph.all_seq.size() >= seq_limit &&
                comp_graph.all_seq_h.size() >= seq_limit) {
              break;
            }
          }
        }

        // add, 0610
        if (_LOOSE){
          for (auto& fp : all_fp_T) {
            if (seq_limit != 0 && comp_graph.all_seq.size() >= seq_limit &&
                comp_graph.all_seq_h.size() >= seq_limit) {
              break;
            }

            if (!_LOOSE || (_LOOSE && fp.dir != DAG_direction_status::LOOSE_possible)) continue;

            vector<map_tree> all_seq_T1;
            vector<map_tree> all_seq_T2;

            all_seq_T1.clear();
            all_seq_T2.clear();

            map_tree seq;

            Gen_Path_T(map_DAG_T, map_DAG_W, map_T2_W3, comp_graph.map_FT_W1,
                       comp_graph.map_FT_W2, comp_graph.map_FT_W3,
                       comp_graph.map_FT_W2_core, fp.node1, seq, all_seq_T1, 
                       global_path_num);
            Gen_Path_T(map_DAG_T, map_DAG_W, map_T2_W3, comp_graph.map_FT_W1,
                       comp_graph.map_FT_W2, comp_graph.map_FT_W3,
                       comp_graph.map_FT_W2_core, fp.node2, seq, all_seq_T2, 
                       global_path_num);

            size_t h =
                max(fp.node1.rsT.second.height, fp.node2.rsT.second.height);

            for (auto& seq1 : all_seq_T1) {
              for (auto& seq2 : all_seq_T2) {
                if (h == core_height) {
                  // if (fp.dir == DAG_direction_status::EXACT_possible)
                  //   comp_graph.all_seq_h.emplace_back(
                  //       seq1, seq2, fp.k, max(seq1.height, seq2.height), false);
                  // if (_LOOSE)
                    comp_graph.all_seq_h_loose.emplace_back(
                        seq1, seq2, fp.k, max(seq1.height, seq2.height), false);
                  // cout << "seq1.trees.size() = " << seq1.trees.size() << "
                  // seq1.mul.size() = " << seq1.mul.size() << endl; seq1.print();
                  // cout << "seq2.trees.size() = " << seq2.trees.size() << "
                  // seq2.mul.size() = " << seq2.mul.size() << endl; seq2.print();
                  if (seq_limit != 0 &&
                      comp_graph.all_seq_h.size() >= seq_limit) {
                    break;
                  }
                } else {
                  // if (fp.dir == DAG_direction_status::EXACT_possible)
                  //   comp_graph.all_seq.emplace_back(
                  //       seq1, seq2, fp.k, max(seq1.height, seq2.height), false);
                  // if (_LOOSE)
                    comp_graph.all_seq_loose.emplace_back(
                        seq1, seq2, fp.k, max(seq1.height, seq2.height), false);

                  // cout << "seq1.trees.size() = " << seq1.trees.size() << "
                  // seq1.mul.size() = " << seq1.mul.size() << endl; seq1.print();
                  // cout << "seq2.trees.size() = " << seq2.trees.size() << "
                  // seq2.mul.size() = " << seq2.mul.size() << endl; seq2.print();
                  if (seq_limit != 0 && comp_graph.all_seq.size() >= seq_limit) {
                    break;
                  }
                }
              }
              if (seq_limit != 0 && comp_graph.all_seq.size() >= seq_limit &&
                  comp_graph.all_seq_h.size() >= seq_limit) {
                break;
              }
            }
          }
        }

       // cout << "time to traverse to get all graphs from DAG = "
        //     << DAG_traverse_graph_time << endl;

        double time_tmp_9 = globalTimeKeeper::tt.diff();
        DAG_traverse_graph_time += time_tmp_9 - time_tmp_9_3;

        if (IN_INFO.height == core_height) {
          cp.with_core_height = true;
        }

        if (comp_graph.all_seq_h.size() + comp_graph.all_seq.size() == 0) {
          cp.no_feasible_pair = true;
          if (cp.with_core_height) {
            comp_graph.num_h = 1;
            comp_graph.num = 0;
          } else {
            comp_graph.num = 1;
            comp_graph.num_h = 0;
          }
        } else {
          if (comp_graph.num_h < comp_graph.all_seq_h.size()){
            comp_graph.num_h = comp_graph.all_seq_h.size();
          } 
          if (comp_graph.num < comp_graph.all_seq.size()){
            comp_graph.num = comp_graph.all_seq.size();
          }
        }

        if (comp_graph.all_seq_h_loose.size() + comp_graph.all_seq_loose.size() == 0 && cp.no_feasible_pair) {
          cp.no_feasible_pair_loose = true;
        }
      }
    }
  }
  // cout << "start generating" << endl;
  size_t lower_bound = 0;
  size_t possible_num = 0;

  double time_tmp_1 = globalTimeKeeper::tt.diff();

  if (!core_height_fixed) {
    Gen_Graph_unfixed_ch(g, core_set, set_tree_component, set_edge_component,
                         g_n, g_m, possible_num, lower_bound, outputfilename);
    if (_LOOSE) Gen_Graph_unfixed_ch_loose(g, core_set, set_tree_component, set_edge_component,
                         g_n, g_m, possible_num, lower_bound, outputfilename, fc_map, fc_all, TS_map);
  } else {
    Gen_Graph(g, core_set, set_tree_component, set_edge_component, g_n, g_m,
              possible_num, lower_bound, outputfilename);
    if (_LOOSE) Gen_Graph_loose(g, core_set, set_tree_component, set_edge_component, g_n, g_m,
              possible_num, lower_bound, outputfilename, fc_map, fc_all, TS_map);
  }

  double _ts8 = tk.tt.diff();

  if (lower_bound < total_num) lower_bound = total_num;

  if (_EXP) {
    cout << "# fringe trees = " << fringe_tree_sum << endl;
    cout << "# fringe tree vectors = " << fringe_tree_vectors_sum << endl;
    cout << "time for fringe tree vectors = " << fringe_tree_vectors_time
         << endl;
    cout << "# end-subtree vectors = " << end_subtree_vectors_sum << endl;
    cout << "time for end-subtree vectors = " << end_subtree_vectors_time
         << endl;
    cout << "# rooted core-subtrees = " << rooted_core_subtrees_sum << endl;
    cout << "time for rooted core-subtrees = " << rooted_core_subtrees_time
         << endl;
    cout << "# bi-rooted core-subtrees = " << bi_rooted_core_subtrees_sum
         << endl;
    cout << "time for bi-rooted core-subtrees = "
         << bi_rooted_core_subtrees_time << endl;
    cout << "# feasible pairs = " << feasible_pairs_sum << endl;
    cout << "time for feasible pairs = " << feasible_pairs_time << endl;
    cout << "# graphs generated = " << total_num << endl;
    cout << "time for generating graphs = " << _ts8 - time_tmp_1 << endl;
    cout << "size of DAG, |V| = " << DAG_size_v << endl;
    cout << "size of DAG, |E| = " << DAG_size_e << endl;
    cout << "time to construct DAG = " << DAG_construct_time << endl;
    cout << "time to traverse all paths from DAG = " << DAG_traverse_path_time
         << endl;
    cout << "time to traverse to get all graphs from DAG = "
         << DAG_traverse_graph_time << endl;
    cout << "time to get exact number of graphs from DAG = "
         << DAG_traverse_num_time << endl;
    cout << endl;
  }

  cout << "A lower bound on the number of graphs = " << lower_bound << endl;
  //cout << "Number of possible graphs to generate = " << possible_num << endl;
  cout << "Number of generated graphs = " << total_num << endl;
  cout << "Total time : " << _ts8 - _ts0 << "s." << endl;
  // cout << "global path num : " << global_path_num << endl;
  // cout << "global_path_limit : " << global_path_limit << endl;

  //cout << "***********************************" << endl;
  return 0;
}