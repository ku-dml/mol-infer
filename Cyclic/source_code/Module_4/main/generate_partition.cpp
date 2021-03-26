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

// #include "../include/cross_timer.h"
#include "../include/data_structures.hpp"
#include "../include/chemical_graph.hpp"
// #include "../include/fringe_tree.hpp"
// #include "../include/tools.hpp"

using namespace std;

// This is the file to generate a standard partition file given a SDF file.

int main(int argc, char** argv){


    if (argc != 3) {
        cout << "Please supply: " << endl;
        cout << "1. an SDF file with one cyclic graph " << endl;
        cout << "2. a file name for the output " << endl;
        exit(-1);
    }

    string infile;
    try{
        infile = string(argv[1]);
    } catch (std::exception& e){
        cout << "Caught an exception for the filename argument:" << endl;
        cout << e.what() << endl;
        cout << "The program will now exit." << endl;
        return 1;
    }

    string outputfilename = string(argv[2]);

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

    vector <Vertex> base_vertices;
    vector_2D <Vertex> base_edges;
    
    get_partition(g, core_set, base_vertices, base_edges);

    vector <component> set_tree_component;
    for (size_t i = 0; i < base_vertices.size(); ++i){
        auto& v = base_vertices[i];
        component _cp;
        _cp.no_feasible_pair = false;
        _cp.with_core_height = false;
        // input_info IN_INFO;
        get_descriptors_tree(_cp.IN_INFO, g, core_set, internal_set, v);

        set_tree_component.push_back(_cp);
    }

    vector <component> set_edge_component;
    for (size_t i = 0; i < base_edges.size(); ++i){
        auto& e = base_edges[i];
        component _cp;
        _cp.no_feasible_pair = false;
        _cp.with_core_height = false;
        get_descriptors_edge(_cp.IN_INFO, g, core_set, internal_set, e);

        set_edge_component.push_back(_cp);
    }
    
    output_partition(g, outputfilename, base_vertices, base_edges, set_tree_component, set_edge_component);
            
    return 0;
}