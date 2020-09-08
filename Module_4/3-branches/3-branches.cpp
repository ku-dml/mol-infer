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

#include "../include/cross_timer.h"
#include "../include/data_structures.hpp"
#include "../include/fringe_tree.hpp"
#include "../include/debug.h"

using namespace std;

// initialize a global timekeeper
Timer globalTimeKeeper::tt;

int main(int argc, char** argv){
    if (argc < 6) {
        cout << "Please supply: " << endl;
        cout << "1. a file of resources " << endl;
        cout << "2. an integer of time limit" << endl;
        cout << "3. an integer of vector size bound per iteration" << endl;
        cout << "4. an integer of number of wanted graphs" << endl;
        cout << "5. a file name of output" << endl;
        cout << "6. (optional) a file name of a sample sdf with given resources" << endl;
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
    num_limit = stoi(argv[4]);

    // cout  <<  infile << endl;
    input_info IN_INFO;
    read_input(infile, IN_INFO);

    if (argc >= 7){
        string input_sdf_filename = string(argv[6]);
        get_F1(IN_INFO, input_sdf_filename, T_W1_F1, T_W2_F1, T_W3_F1);
    }

    size_t M = IN_INFO.num_kind_atoms;

    // IN_INFO.print();

    double _ts0 = tk.tt.diff();

    vector_2D <PathTrie> all_trie_node(IN_INFO.k1 + 2);

    for (size_t i = 0; i <= IN_INFO.k1 + 1; ++i){
        PathTrie tmp(0);
        all_trie_node[i].emplace_back(tmp);
    }

    vector <unordered_set <root_status>> all_rs_W1(IN_INFO.k1);
    unordered_set <root_status> all_rs_W1_F1;
    vector <unordered_map <root_status, map <size_t, pair <map_rv, size_t>>>> map_W1(IN_INFO.k1);
    unordered_map <root_status, unordered_set <size_t>> map_W1_F1;
    unordered_map <pair <unsigned short, size_t>, vector <size_t>> map_FT_W1;
    vector_2D <FringeTree> All_FT_W1(M + 1);
    
    GenW1(IN_INFO, all_rs_W1, all_rs_W1_F1, map_FT_W1, map_W1, map_W1_F1, All_FT_W1, all_trie_node);   // Genenrate W1(a,d,m)

    double _ts1 = tk.tt.diff();
    // cout << "Time : " << _ts1 - _ts0 << "s." << endl;

    vector <unordered_set <half_path_status>> all_rs_W2(1);
    unordered_set <half_path_status> all_hps_F1;
    vector <unordered_map <half_path_status, map <size_t, pair <map_rv, size_t>>>> map_W2(1);
    unordered_map <half_path_status, unordered_set <size_t>> map_W2_F1;
    unordered_map <half_path_status, vector <size_t>> map_W2_ind;
    unordered_map <pair <unsigned short, size_t>, vector <size_t>> map_FT_W2;
    vector_2D <FringeTree> All_FT_W2(M + 1);

    GenW2_0(IN_INFO, all_rs_W2, all_hps_F1, map_FT_W2, map_W2, map_W2_F1, map_W2_ind, All_FT_W2, all_trie_node);   // Generate W2(0)(a,d,m)
    double _ts2 = tk.tt.diff();
    // cout << "Time : " << _ts2 - _ts1 << "s." << endl;

    merge_W1(IN_INFO, all_rs_W1, all_rs_W1_F1, all_rs_W2, all_hps_F1, map_W1, map_W1_F1, map_W2, map_W2_F1, map_W2_ind, all_trie_node);   // only need call once

    double _ts3 = tk.tt.diff();
    // cout << "Time : " << _ts3 - _ts2 << "s." << endl;

    vector <unordered_set <root_status>> all_rs_W3(IN_INFO.k1 + 1);
    unordered_set <root_status> all_rs_W3_F1;
    vector <unordered_map <root_status, map <size_t, pair <map_rv, size_t>>>> map_W3(IN_INFO.k1 + 1);
    unordered_map <root_status, unordered_set <size_t>> map_W3_F1;
    unordered_map <pair <unsigned short, size_t>, vector <size_t>> map_FT_W3;
    vector_2D <FringeTree> All_FT_W3(M + 1);

    GenW3_0(IN_INFO, all_rs_W3, all_rs_W3_F1, map_FT_W3, map_W3, map_W3_F1, All_FT_W3, all_trie_node);
    double _ts4 = tk.tt.diff();
    // cout << "Time : " << _ts4 - _ts3 << "s." << endl;

    merge_W3(IN_INFO, all_rs_W1, all_rs_W3, all_rs_W3_F1, map_W1, map_W3, map_W3_F1, all_trie_node); 
    double _ts5 = tk.tt.diff();
    // cout << "Time : " << _ts5 - _ts4 << "s." << endl;

    vector <unordered_set <root_status>> all_rs_W4(IN_INFO.k1 + 1);
    vector <unordered_map <root_status, map <size_t, pair <map_rv, size_t>>>> map_W4(IN_INFO.k1 + 1);
    // cout << "Calculating the path combined by l1 and l2, denoted by W4!" <<endl;
    merge_W4(IN_INFO, all_rs_W1, all_rs_W3, all_rs_W4, map_W1, map_W3, map_W4, all_trie_node); 
    double _ts6 = tk.tt.diff();
    // cout << "Time : " << _ts6 - _ts5 << "s." << endl;

    // check number of feasible pairs
    vector_2D <map_rv> all_seq(IN_INFO.k1 + 1);
    combine_check(IN_INFO, all_rs_W1, all_rs_W4, map_W1, map_W4, all_trie_node, all_seq);
    double _ts7 = tk.tt.diff();
    // cout << "Time : " << _ts7 - _ts6 << "s." << endl;

    // string outputfilename = infile.substr(0, infile.size() - 4) + ".sdf";
    string outputfilename = argv[5];

    Gen_Graph(IN_INFO, all_seq, map_FT_W1, map_FT_W2, map_FT_W3, All_FT_W1, All_FT_W2, All_FT_W3, outputfilename);

    double _ts8 = tk.tt.diff();
    // cout << "Time : " << _ts8 - _ts7 << "s." << endl;
    cout << "Time : " << _ts8 - _ts0 << "s." << endl;
            
    return 0;
}