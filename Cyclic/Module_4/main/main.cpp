#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <set>

bool _DEBUG = false; // a global boolean variable to show message for debugging
bool _EXP = false;

#include "../include/cross_timer.h"
#include "../include/data_structures.hpp"
#include "../include/chemical_graph.hpp"
#include "../include/fringe_tree.hpp"
#include "../include/tools.hpp"

using namespace std;

// initialize a global timekeeper
Timer globalTimeKeeper::tt;

int main(int argc, char** argv){
    
    /////////////////////////////////
    // change this value  to true when debugging
    _DEBUG = false;
    _EXP = false;
    ////////////////////////////////


    if (argc < 6) {
        cout << "Please supply: " << endl;
        cout << "1. an SDF file with one cyclic graph " << endl;
        cout << "2. an integer of time limit per iteration" << endl;
        cout << "3. an integer of vector size bound per iteration (0 if no limit)" << endl;
        cout << "4. an integer of number of sample graphs stored for each vector" <<  endl;
        cout << "5. an integer of number of wanted graphs" << endl;
        cout << "6. a file name of output" << endl;
        cout << "7. a file with information of partition (optional now)" << endl;
        exit(-1);
    }
    globalTimeKeeper tk;

    string infile;
    try{
        infile = string(argv[1]);
    } catch (std::exception& e){
        cout << "Caught an exception for the filename argument:" << endl;
        cout << e.what() << endl;
        cout << "The program will now exit." << endl;
        return 1;
    }

    time_limit = stoi(argv[2]);
    UB_limit = stoull(argv[3]);
    seq_limit = stoi(argv[4]);
    num_limit = stoi(argv[5]);

    string outputfilename = argv[6];
    string partitionfilename;

    bool no_partition_info = false;
    if (argc == 7){
        no_partition_info = true;
    } else {
        partitionfilename = argv[7];
    }

    Graph h = read_graph_sdf(infile);
    Graph g;
    H_suppressed_convert(h, g);
    vector <Vertex> core_set;
    vector <Vertex> internal_set;
    calcCoreVertexSet(g, core_set);
    calcInternalVertexSet(g, internal_set);
    size_t g_n = calcEffectiveVertexNum(g);
    size_t g_m = calcNumOfEdges(g);
    size_t core_height = 0;
    size_t bc = 0;

    double _ts0 = tk.tt.diff();

    vector <Vertex> base_vertices;
    vector_2D <Vertex> base_edges;
    vector <size_t> chLB_v;
    vector <size_t> chUB_v;
    vector <size_t> chLB_e;
    vector <size_t> chUB_e;

    vector <bool> fixed_v;
    vector <bool> fixed_e;
    
    if (no_partition_info){
        get_partition(g, core_set, base_vertices, base_edges);
        fixed_v = vector <bool>(base_vertices.size(), false);
        fixed_e = vector <bool>(base_edges.size(), false);
        // get_partition(g, core_set, internal_set, base_vertices, base_edges);
    } else {
        get_partition(base_vertices, base_edges, chLB_v, chUB_v, chLB_e, chUB_e, fixed_v, fixed_e, partitionfilename);
    }
    
    if (_DEBUG){
        cout << "DEBUG on !!!" << endl;
    }
    if (_DEBUG) {
        cout << "|V_B| = " << base_vertices.size() << " |E_B| = " << base_edges.size() << endl;
    }

    cout.unsetf(ios::floatfield);
    cout.precision(3);

    vector <component> set_tree_component;
    for (size_t i = 0; i < base_vertices.size(); ++i){
        auto& v = base_vertices[i];
        component _cp;
        _cp.no_feasible_pair = false;
        _cp.with_core_height = false;
        // input_info IN_INFO;
        get_descriptors_tree(_cp.IN_INFO, g, core_set, internal_set, v);

        if (no_partition_info){
            _cp.IN_INFO.chLB = 0;
            _cp.IN_INFO.chUB = _cp.IN_INFO.height;
            // _cp.IN_INFO.chUB = 100;
        } else {
            _cp.IN_INFO.chLB = chLB_v[i];
            _cp.IN_INFO.chUB = chUB_v[i];
        }

        if (_cp.IN_INFO.height > core_height){
            core_height = _cp.IN_INFO.height;
        }

        if (fixed_v[i]){
            _cp.no_feasible_pair = true;
        }

        bc += _cp.IN_INFO.rv.bl;

        set_tree_component.push_back(_cp);
        if (_DEBUG) cout << "v = " << v << endl;
    }

    vector <component> set_edge_component;
    for (size_t i = 0; i < base_edges.size(); ++i){
        auto& e = base_edges[i];
        component _cp;
        _cp.no_feasible_pair = false;
        _cp.with_core_height = false;
        get_descriptors_edge(_cp.IN_INFO, g, core_set, internal_set, e);

        if (no_partition_info){
            _cp.IN_INFO.chLB = 0;
            _cp.IN_INFO.chUB = _cp.IN_INFO.height;
            // _cp.IN_INFO.chUB = 100;
        } else {
            _cp.IN_INFO.chLB = chLB_e[i];
            _cp.IN_INFO.chUB = chUB_e[i];
        }

        if (_cp.IN_INFO.height > core_height){
            core_height = _cp.IN_INFO.height;
        }

        if (fixed_e[i]){
            _cp.no_feasible_pair = true;
        }

        bc += _cp.IN_INFO.rv.bl;

        set_edge_component.push_back(_cp);
        // if (_DEBUG) {
        //     cout << "e = ";
        //     for (auto& u : e){
        //         cout << u << " ";
        //     }
        //     cout<<endl;
        // }
    }
    if (_EXP){
        unordered_set <string> Lambda;
        unordered_set <string> Gamma_in;
        unordered_set <string> Gamma_ex;
        unordered_set <string> Gamma_co;

        for (auto& cp : set_tree_component){
            auto& IN_INFO = cp.IN_INFO;
            for (auto& atom : IN_INFO.atoms){
                Lambda.insert(atom.NAME);
            }
            size_t M = IN_INFO.num_kind_atom_degree;
            for (size_t i = 0; i < M; ++i){
                for (size_t j = 0; j <= i; ++j){
                    for (size_t k = 1; k <= 3; ++k){
                        if (IN_INFO.rv.ec_co[i][j][k] != 0){
                            string tmp = "";
                            atom_degree ad1 = IN_INFO.atom_deg[i];
                            tmp += IN_INFO.atoms[ad1.color].NAME + to_string(ad1.deg);
                            atom_degree ad2 = IN_INFO.atom_deg[j];
                            tmp += IN_INFO.atoms[ad2.color].NAME + to_string(ad2.deg);
                            tmp += to_string(k);
                            Gamma_co.insert(tmp);
                        }
                    }
                }
            }
            for (size_t i = 0; i < M; ++i){
                for (size_t j = 0; j < M; ++j){
                    for (size_t k = 1; k <= 3; ++k){
                        if (IN_INFO.rv.ec_in[i][j][k] != 0){
                            string tmp = "";
                            atom_degree ad1 = IN_INFO.atom_deg[i];
                            tmp += IN_INFO.atoms[ad1.color].NAME + to_string(ad1.deg);
                            atom_degree ad2 = IN_INFO.atom_deg[j];
                            tmp += IN_INFO.atoms[ad2.color].NAME + to_string(ad2.deg);
                            tmp += to_string(k);
                            Gamma_in.insert(tmp);
                        }
                    }
                }
            }
            for (size_t i = 0; i < M; ++i){
                for (size_t j = 0; j < M; ++j){
                    for (size_t k = 1; k <= 3; ++k){
                        if (IN_INFO.rv.ec_ex[i][j][k] != 0){
                            string tmp = "";
                            atom_degree ad1 = IN_INFO.atom_deg[i];
                            tmp += IN_INFO.atoms[ad1.color].NAME + to_string(ad1.deg);
                            atom_degree ad2 = IN_INFO.atom_deg[j];
                            tmp += IN_INFO.atoms[ad2.color].NAME + to_string(ad2.deg);
                            tmp += to_string(k);
                            Gamma_ex.insert(tmp);
                        }
                    }
                }
            }
        }
        for (auto& cp : set_edge_component){
            auto& IN_INFO = cp.IN_INFO;
            for (auto& atom : IN_INFO.atoms){
                Lambda.insert(atom.NAME);
            }
            size_t M = IN_INFO.num_kind_atom_degree;
            for (size_t i = 0; i < M; ++i){
                for (size_t j = 0; j <= i; ++j){
                    for (size_t k = 1; k <= 3; ++k){
                        if (IN_INFO.rv.ec_co[i][j][k] != 0){
                            string tmp = "";
                            atom_degree ad1 = IN_INFO.atom_deg[i];
                            tmp += IN_INFO.atoms[ad1.color].NAME + to_string(ad1.deg);
                            atom_degree ad2 = IN_INFO.atom_deg[j];
                            tmp += IN_INFO.atoms[ad2.color].NAME + to_string(ad2.deg);
                            tmp += to_string(k);
                            Gamma_co.insert(tmp);
                            // cout << tmp <<  endl;
                        }
                    }
                }
            }
            for (size_t i = 0; i < M; ++i){
                for (size_t j = 0; j < M; ++j){
                    for (size_t k = 1; k <= 3; ++k){
                        if (IN_INFO.rv.ec_in[i][j][k] != 0){
                            string tmp = "";
                            atom_degree ad1 = IN_INFO.atom_deg[i];
                            tmp += IN_INFO.atoms[ad1.color].NAME + to_string(ad1.deg);
                            atom_degree ad2 = IN_INFO.atom_deg[j];
                            tmp += IN_INFO.atoms[ad2.color].NAME + to_string(ad2.deg);
                            tmp += to_string(k);
                            Gamma_in.insert(tmp);
                        }
                    }
                }
            }
            for (size_t i = 0; i < M; ++i){
                for (size_t j = 0; j < M; ++j){
                    for (size_t k = 1; k <= 3; ++k){
                        if (IN_INFO.rv.ec_ex[i][j][k] != 0){
                            string tmp = "";
                            atom_degree ad1 = IN_INFO.atom_deg[i];
                            tmp += IN_INFO.atoms[ad1.color].NAME + to_string(ad1.deg);
                            atom_degree ad2 = IN_INFO.atom_deg[j];
                            tmp += IN_INFO.atoms[ad2.color].NAME + to_string(ad2.deg);
                            tmp += to_string(k);
                            Gamma_ex.insert(tmp);
                        }
                    }
                }
            }
        }
        cout << "# vertices = " << g_n << endl;
        cout << "|V_B| = " << base_vertices.size() << endl;
        cout << "|E_B| = " << base_edges.size() << endl;
        cout << "|Gamma^in| = " << Gamma_in.size() << endl;
        cout << "|Gamma^ex| = " << Gamma_ex.size() << endl;
        cout << "|Gamma^co| = " << Gamma_co.size() << endl;
        cout << "|Lambda| = " << Lambda.size() << endl;
        cout << "bc = " << bc << endl;
        cout << "ch = " << core_height << endl;
    }

    for (auto& cp : set_tree_component){
        auto& IN_INFO = cp.IN_INFO;
        auto& comp_graph = cp.comp_graph;

        if (cp.no_feasible_pair){
            if (IN_INFO.height == core_height){
                cp.with_core_height = true;
                comp_graph.num = 0;
                comp_graph.num_h = 1;
            } else {
                comp_graph.num = 1;
                comp_graph.num_h = 0;
            }
            continue;
        }

        size_t ch = 0;
        if (IN_INFO.height >= IN_INFO.k_star){
            ch = IN_INFO.height - IN_INFO.k_star;
        }
        size_t M = IN_INFO.num_kind_atoms;
        if (_DEBUG) cout << "ch = " << ch << " IN_INFO.height = " << IN_INFO.height << " v = " << IN_INFO.base_vertices[0] << " delta_i = " << IN_INFO.delta_1 << endl;
        // print_rv1D(IN_INFO, IN_INFO.rv1D);
        vector_2D <PathTrie> all_trie_node(ch + 1);
        for (size_t i = 0; i <= ch; ++i){
            PathTrie tmp(0);
            all_trie_node[i].emplace_back(tmp);
        }

        // 
        // W1: end fringe tree, always +1
        // W2: internal fringe tree,always +2
        // W3: fringe tree attached to the core vertex, +delta_i+1
        // W2_core: fringe tree attached to the core vertex, +delta_i
        vector <unordered_set <root_status>> all_rs_W1(ch);
        vector <unordered_map <root_status, map <size_t, pair <vector <map_rv>, size_t>>>> map_W1(ch);

        unordered_set <half_path_status> all_rs_W2;
        unordered_map <half_path_status, map <size_t, pair <vector <map_rv>, size_t>>> map_W2;
        unordered_map <half_path_status, vector <size_t>> map_W2_ind;

        vector <unordered_set <root_status>> all_rs_W3(ch + 1);
        vector <unordered_map <root_status, map <size_t, pair <vector <map_rv>, size_t>>>> map_W3(ch + 1);

        unordered_set <root_status_height> all_rs_W2_core;
        unordered_map <root_status_height, map <size_t, pair <vector <map_rv>, size_t>>> map_W2_core;

        comp_graph.initialize(M + 1);

        size_t chl1 = (ch - 1) / 2;
        size_t chl2 = ch - 1 - chl1;

        if (ch > 0){
            double time_tmp_1 = globalTimeKeeper::tt.diff();

            GenW1(IN_INFO, all_rs_W1, comp_graph.map_FT_W1, map_W1, comp_graph.All_FT_W1, all_trie_node);
            GenW2_0(IN_INFO, all_rs_W2, comp_graph.map_FT_W2, map_W2, map_W2_ind, comp_graph.All_FT_W2, all_trie_node);
            GenW3_0(IN_INFO, all_rs_W3, comp_graph.map_FT_W3, map_W3, comp_graph.All_FT_W3, all_trie_node, IN_INFO.target_rs1.color, IN_INFO.delta_1 + 1);
            
            double time_tmp_2 = globalTimeKeeper::tt.diff();

            // size_t chl1 = (ch - 1) / 2;
            merge_W1(IN_INFO, all_rs_W1, all_rs_W2, map_W1, map_W2, map_W2_ind, all_trie_node, chl1);

            double time_tmp_3 = globalTimeKeeper::tt.diff();

            // merge_W3(IN_INFO, all_rs_W1, all_rs_W3, map_W1, map_W3, all_trie_node, ch, ch, IN_INFO.delta_1);
            // size_t chl2 = ch - 1 - chl1;
            merge_W3_tree(IN_INFO, all_rs_W2, all_rs_W3, map_W2, map_W3, all_trie_node, 0, chl2, IN_INFO.delta_1);

            double time_tmp_4 = globalTimeKeeper::tt.diff();
        
            fringe_tree_vectors_time += time_tmp_2 - time_tmp_1;
            end_subtree_vectors_time += time_tmp_3 - time_tmp_2;
            // rooted_core_subtrees_time += time_tmp_4 - time_tmp_3;
            end_subtree_vectors_time += time_tmp_4 - time_tmp_3;
        } else {
            double time_tmp_1 = globalTimeKeeper::tt.diff();

            GenW2_0_core(IN_INFO, all_rs_W2_core, comp_graph.map_FT_W2_core, map_W2_core, comp_graph.All_FT_W2_core, all_trie_node, IN_INFO.delta_1);
        
            double time_tmp_2 = globalTimeKeeper::tt.diff();

            fringe_tree_vectors_time += time_tmp_2 - time_tmp_1;
        }

        //////////////////////////////////////////////////
        if (ch > 0){
            // vector <map_rv> all_seq_h_tmp;
            double time_tmp_1 = globalTimeKeeper::tt.diff();

            vector <map_rv> all_seq_tmp;
            size_t num_tmp = 0;
            combine_check_tree_component(IN_INFO, all_rs_W1, all_rs_W3, map_W1, map_W3, all_trie_node,
                all_seq_tmp, chl1, chl2, IN_INFO.target_rs1, core_height, num_tmp);

            if (IN_INFO.height == core_height){
                cp.with_core_height = true;
                comp_graph.num_h = num_tmp;
                if (all_seq_tmp.size() > 0){
                    for (auto& w : all_seq_tmp){
                        comp_graph.all_seq_h.emplace_back(w, IN_INFO.height);
                    } 
                } else {
                    cp.no_feasible_pair = true;
                    comp_graph.num_h = 1;
                }
                comp_graph.num = 0;
                // cout << "comp_graph.all_seq_h.size() = " << comp_graph.all_seq_h.size() << endl; 
            } else {
                comp_graph.num = num_tmp;
                if (all_seq_tmp.size() > 0){
                    for (auto& w : all_seq_tmp){
                        comp_graph.all_seq.emplace_back(w, IN_INFO.height);
                    }
                } else {
                    cp.no_feasible_pair = true;
                    comp_graph.num = 1;
                }    
                comp_graph.num_h = 0;
                // cout << "comp_graph.all_seq.size() = " << comp_graph.all_seq.size() << endl;    
            }
            double time_tmp_2 = globalTimeKeeper::tt.diff();

            feasible_pairs_time += time_tmp_2 - time_tmp_1;

        } else {
            vector <PathTrie> all_trie_node_path;
            PathTrie tmp(0);
            all_trie_node_path.emplace_back(tmp);
            unordered_set <root_status_height> all_rs_T2;
            unordered_map <root_status_height, map <size_t, pair <vector <map_tree>, size_t>>> map_T2;
            unordered_map <root_status_height, vector <size_t>> map_T2_ind;

            double time_tmp_1 = globalTimeKeeper::tt.diff();

            GenT2_0(IN_INFO, all_rs_T2, map_T2, map_T2_ind, all_trie_node_path, all_rs_W2_core, map_W2_core, all_rs_W3, map_W3, all_trie_node, IN_INFO.height, true);
        
            double time_tmp_2 = globalTimeKeeper::tt.diff();

            rooted_core_subtrees_time += time_tmp_2 - time_tmp_1;

            for (auto& rsh : all_rs_T2){
                for (auto& w : map_T2.at(rsh)){
                    if (rsh.height == core_height){
                        for (auto& _seq : w.second.first){
                            comp_graph.all_seq_h.emplace_back(_seq, rsh.height);
                        }
                        comp_graph.num_h += w.second.second;
                    } else {
                        for (auto& _seq : w.second.first){
                            comp_graph.all_seq.emplace_back(_seq, rsh.height);
                        } 
                        comp_graph.num += w.second.second;
                    }
                }
            }
        }
    }

    for (auto& cp : set_edge_component){
        auto& IN_INFO = cp.IN_INFO;
        auto& comp_graph = cp.comp_graph;

        if (cp.no_feasible_pair){
            if (IN_INFO.height == core_height){
                cp.with_core_height = true;
                comp_graph.num = 0;
                comp_graph.num_h = 1;
            } else {
                comp_graph.num = 1;
                comp_graph.num_h = 0;
            }
            continue;
        }

        size_t chlb = 0;
        if (IN_INFO.chLB >= IN_INFO.k_star){
            chlb = IN_INFO.chLB - IN_INFO.k_star;
        }
        size_t chub = 0;
        if (IN_INFO.chUB >= IN_INFO.k_star){
            chub = IN_INFO.chUB - IN_INFO.k_star;
        }
        
        size_t M = IN_INFO.num_kind_atoms;
        size_t len = IN_INFO.base_vertices.size() - 1;
        size_t l2 = (len - 1) / 2;
        size_t l1 = len - 1 - l2;

        if (_DEBUG) cout << "IN_INFO.base_vertices[0] = " << IN_INFO.base_vertices[0] << " IN_INFO.base_vertices[IN_INFO.base_vertices.size() - 1] = " << IN_INFO.base_vertices[IN_INFO.base_vertices.size() - 1] << endl;
        if (_DEBUG) cout << "chlb = " << chlb << " IN_INFO.chLB = " << IN_INFO.chLB << " chub = " << chub << " IN_INFO.chUB = " << IN_INFO.chUB << endl;
        if (_DEBUG) cout << "l1 = " << l1 << " l2 = " << l2 << endl;
        // print_rv1D(IN_INFO, IN_INFO.rv1D);
        vector_2D <PathTrie> all_trie_node(chub + 1);
        for (size_t i = 0; i <= chub; ++i){
            PathTrie tmp(0);
            all_trie_node[i].emplace_back(tmp);
        }

        // 
        // W1: end fringe tree, always +1
        // W2: internal fringe tree,always +2
        // W3: fringe tree attached to the core vertex, +delta_i+1
        // W2_core: fringe tree attached to the core vertex, +delta_i
        vector <unordered_set <root_status>> all_rs_W1(chub);
        vector <unordered_map <root_status, map <size_t, pair <vector <map_rv>, size_t>>>> map_W1(chub);

        unordered_set <half_path_status> all_rs_W2;
        unordered_map <half_path_status, map <size_t, pair <vector <map_rv>, size_t>>> map_W2;
        unordered_map <half_path_status, vector <size_t>> map_W2_ind;

        vector <unordered_set <root_status>> all_rs_W3(chub + 1);
        vector <unordered_map <root_status, map <size_t, pair <vector <map_rv>, size_t>>>> map_W3(chub + 1);

        unordered_set <root_status_height> all_rs_W2_core;
        unordered_map <root_status_height, map <size_t, pair <vector <map_rv>, size_t>>> map_W2_core;

        comp_graph.initialize(M + 1);

        if (chub > 0){
            double time_tmp_1 = globalTimeKeeper::tt.diff();

            GenW1(IN_INFO, all_rs_W1, comp_graph.map_FT_W1, map_W1, comp_graph.All_FT_W1, all_trie_node);
            GenW2_0(IN_INFO, all_rs_W2, comp_graph.map_FT_W2, map_W2, map_W2_ind, comp_graph.All_FT_W2, all_trie_node);
            GenW3_0(IN_INFO, all_rs_W3, comp_graph.map_FT_W3, map_W3, comp_graph.All_FT_W3, all_trie_node, M, 3);

            double time_tmp_2 = globalTimeKeeper::tt.diff();

            merge_W1(IN_INFO, all_rs_W1, all_rs_W2, map_W1, map_W2, map_W2_ind, all_trie_node, chub - 1);

            double time_tmp_3 = globalTimeKeeper::tt.diff();

            merge_W3(IN_INFO, all_rs_W1, all_rs_W3, map_W1, map_W3, all_trie_node, 1, chub, 2);

            double time_tmp_4 = globalTimeKeeper::tt.diff();

            fringe_tree_vectors_time += time_tmp_2 - time_tmp_1;
            end_subtree_vectors_time += time_tmp_3 - time_tmp_2;
            rooted_core_subtrees_time += time_tmp_4 - time_tmp_3;
        }

        double time_tmp_1 = globalTimeKeeper::tt.diff();

        GenW2_0_core(IN_INFO, all_rs_W2_core, comp_graph.map_FT_W2_core, map_W2_core, comp_graph.All_FT_W2_core, all_trie_node, 2);
        
        double time_tmp_2 = globalTimeKeeper::tt.diff();
        fringe_tree_vectors_time += time_tmp_2 - time_tmp_1;

        vector <PathTrie> all_trie_node_path;
        PathTrie tmp(0);
        all_trie_node_path.emplace_back(tmp);

        vector <unordered_set <root_status_T>> all_rs_T1(l1 + 1);
        vector <unordered_map <root_status_T, map <size_t, pair <vector <map_tree>, size_t>>>> map_T1(l1 + 1);

        GenT1(IN_INFO, all_rs_T1, map_T1, all_trie_node_path, IN_INFO.target_rs1);
        GenT1(IN_INFO, all_rs_T1, map_T1, all_trie_node_path, IN_INFO.target_rs2);

        unordered_set <root_status_height> all_rs_T2;
        unordered_map <root_status_height, map <size_t, pair <vector <map_tree>, size_t>>> map_T2;
        unordered_map <root_status_height, vector <size_t>> map_T2_ind;

        time_tmp_1 = globalTimeKeeper::tt.diff();
        
        GenT2_0(IN_INFO, all_rs_T2, map_T2, map_T2_ind, all_trie_node_path, all_rs_W2_core, map_W2_core, all_rs_W3, map_W3, all_trie_node, IN_INFO.chUB, false);
        
        time_tmp_2 = globalTimeKeeper::tt.diff();
        rooted_core_subtrees_time += time_tmp_2 - time_tmp_1;

        merge_T1(IN_INFO, all_rs_T1, all_rs_T2, map_T1, map_T2, map_T2_ind, all_trie_node_path, l1);

        double time_tmp_3 = globalTimeKeeper::tt.diff();
        bi_rooted_core_subtrees_time += time_tmp_3 - time_tmp_2;

        combine_check_edge_component(
            IN_INFO, all_rs_T1, map_T1, all_trie_node_path, 
            comp_graph.all_seq_h, comp_graph.all_seq, 
            l1, l2, IN_INFO.target_rs1, IN_INFO.target_rs2, core_height, comp_graph.num_h, comp_graph.num);

        if (IN_INFO.height == core_height){
            cp.with_core_height = true;
        }

        if (comp_graph.all_seq_h.size() + comp_graph.all_seq.size() == 0){
            cp.no_feasible_pair = true;
            if (cp.with_core_height){
                comp_graph.num_h = 1;
                comp_graph.num = 0;
            } else {
                comp_graph.num = 1;
                comp_graph.num_h = 0;
            }
        }

        double time_tmp_4 = globalTimeKeeper::tt.diff();
        feasible_pairs_time += time_tmp_4 - time_tmp_3;        
    }
    // cout << "start generating" << endl;
    size_t lower_bound = 0;
    size_t possible_num = 0;

    double time_tmp_1 = globalTimeKeeper::tt.diff();

    Gen_Graph(g, core_set, set_tree_component, set_edge_component, g_n, g_m, possible_num, lower_bound, outputfilename);

    double _ts8 = tk.tt.diff();

    if (_EXP){
        cout << "# fringe tree vectors = " << fringe_tree_vectors_sum << endl;
        cout << "time for fringe tree vectors = " << fringe_tree_vectors_time << endl;
        cout << "# end-subtree vectors = " << end_subtree_vectors_sum << endl;
        cout << "time for end-subtree vectors = " << end_subtree_vectors_time << endl;
        cout << "# rooted core-subtrees = " << rooted_core_subtrees_sum << endl;
        cout << "time for rooted core-subtrees = " << rooted_core_subtrees_time << endl;
        cout << "# bi-rooted core-subtrees = " << bi_rooted_core_subtrees_sum << endl;
        cout << "time for bi-rooted core-subtrees = " << bi_rooted_core_subtrees_time << endl;
        cout << "# feasible pairs = " << feasible_pairs_sum << endl;
        cout << "time for feasible pairs = " << feasible_pairs_time << endl;
        cout << "# graphs generated = " << total_num << endl;
        cout << "time for generating graphs = " << _ts8 - time_tmp_1 << endl;
        cout << endl;
    }


    cout << "A lower bound on the number of graphs = " << lower_bound << endl;
    cout << "Number of possible graphs to generate = " << possible_num << endl;
    // cout << "Number of generated graphs = " << total_num << endl;
    cout << "Total time : " << _ts8 - _ts0 << "s." << endl;
            
    return 0;
}