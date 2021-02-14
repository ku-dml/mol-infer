// This file contains the code of the necessary functions about chemical graph.

#ifndef __CHEMICAL_GRAPH_HPP__INCLUDED
#define __CHEMICAL_GRAPH_HPP__INCLUDED

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <iomanip>
#include <stack>

#include "data_structures.hpp"

using namespace std;

bool _stop_gen = false;
size_t num_limit = 100;   //  limit of the number of graph to output

size_t total_num = 0;

typedef size_t Vertex;
struct Graph{
    string CID;                         // CID
    size_t numAtom, numBond;            // number of vertices and edges
    vector < vector < Vertex > > adj;   // adjacent vertices
    vector < string > alpha;            // atom types (C, H, O, ...)
    vector < vector < size_t > > beta;  // bond types (1, 2, ...)
    vector < bool > status;             // used in calculating branch height
};

Graph read_graph_sdf(const string& inFileName) {

    Graph read_graphs; // return value
    //open the file
    ifstream infile;
    try {
        infile.open(inFileName);
    } catch (const exception& e) {
                cerr << "Couldn't open '" << inFileName
                     << "' for reading!" << endl;
                throw(e);
            }

    if (! infile.is_open()) {
        cerr << "Error, infile not initialized!" << endl;
        throw(-1);
    }
    //read the file
    string line;
    int flag = 0;

    while (getline(infile, line)) {

        while (flag == 0){
            read_graphs.CID = line;
            // Read and discard the next two lines of the file
            getline(infile, line);
            getline(infile, line);

            int n, m;
            stringstream st;
            char *threeChars = new char[3];
            getline(infile, line);
            st << line;
            st.read(threeChars, 3);
            n = atoi(threeChars);
            st.read(threeChars, 3);
            m = atoi(threeChars);
            delete[] threeChars;

            read_graphs.numAtom = n;
            read_graphs.numBond = m;

            read_graphs.alpha.resize(n);
            read_graphs.status.resize(n);
            // a local variable to store the atom symbol
            string alpha;
            // read atom data
            for (int i = 0; i < n; ++i) {
                getline(infile, line);
                stringstream st1;
                // temporary storage for the 3D data from a line in the sdf file
                char *tdData = new char[31];
                st1 << line;
                // read 3D coordinate data in a temporary string storage
                st1.read(tdData, 31); // read the 3D information from the line
                st1 >> alpha;
                read_graphs.alpha[i] = alpha;
                read_graphs.status[i] = true;
                delete[] tdData;
            }

            // initialize a matrix for storing bond multiplicities
            read_graphs.beta.resize(n);
            for (int i = 0; i < n; ++i) {
                read_graphs.beta[i].resize(n);
                for (int j = 0; j < n; ++j) {
                    read_graphs.beta[i][j] = 0;
                }
            }

            for (int i = 0; i < n; ++i) {
                vector < Vertex > tmp;
                read_graphs.adj.push_back(tmp);
            }

            // read bond data
            for (int i = 0; i < m; ++i) {
                // local variable for vertex ids
                int v1, v2, mul;
                getline(infile, line);
                stringstream st2;
                char *vChars = new char[3];
                st2 << line;
                st2.read(vChars, 3);
                v1 = atoi(vChars); // id, 1-n of vertex
                st2.read(vChars, 3);
                v2 = atoi(vChars); // id, 1-n of vertex
                st2.read(vChars, 3);
                mul = atoi(vChars); // v1-v2 multiplicity
                delete[] vChars;

                // account for off-by-one in the indexing
                read_graphs.adj[v1-1].push_back(v2-1);
                read_graphs.adj[v2-1].push_back(v1-1);
                read_graphs.beta[v1-1][v2-1] = mul;
                read_graphs.beta[v2-1][v1-1] = mul;
            }
            flag = 1;
        }
        // The end of the graph information in the file is marked with a "$$$$"
        if (line == "$$$$" && flag == 1) {
            break;
        }
    }
    return read_graphs;
}

/**
 * Take a chemical graph H and convert it to its H-suppressed model
 * @param h: a graph possibly containing H-atoms
 * @param Graph g: the H-suppressed model of h
 * @return the number of H atoms in h
 */
int
H_suppressed_convert(const Graph &h, Graph &g){

    int numH = 0;
    // get H-vertices from the input graph
    // unordered_set <Vertex> Hvertex, nonHvertex;
    unordered_map <Vertex, Vertex> g_h;
    int count = 0; // count the number on non-H atoms
    // itereate over the atoms in graph h
    for (size_t i = 0; i < h.alpha.size(); ++i) {
        if (h.alpha[i] != "H") {
            g.alpha.push_back(h.alpha[i]);
            g.status.push_back(true);    // by Zhu
            g_h[count] = i;
            count++;
        } else {
            // increase the count of H-atoms
            ++numH;
        }
    }

    g.numAtom = g.alpha.size();
    g.numBond = 0;
    g.beta.resize((g.numAtom), vector < size_t > (g.numAtom));
    for (size_t i = 0; i < g.numAtom; ++i) {
        for (size_t j = 0; j < g.numAtom; ++j) {
            g.beta[i][j] = h.beta[g_h[i]][g_h[j]];
            if (h.beta[g_h[i]][g_h[j]] != 0) {
                g.numBond++;
            }
        }
    }
    g.numBond = g.numBond*0.5;
    g.CID = h.CID;

    g.adj.resize(g.numAtom);
    for (size_t i = 0; i < g.beta.size(); ++i) {
        for (size_t j = 0; j < g.beta.size(); ++j) {
            if (g.beta[i][j] != 0) {
                g.adj[i].push_back(Vertex(j));
            }
        }
    }
    return numH;
}

int calcEffectiveVertexNum(const Graph& g){
    int numv = 0;
    for (Vertex u = 0; u < g.numAtom; ++u) {
        if (g.status[u] == true){
            numv++;
        }
    }

    return numv;
}

int calcNumOfEdges(const Graph& g){

    int nume = 0;
    for (Vertex u = 0; u < g.numAtom; u++) {
        nume += g.adj[u].size();
    }
    return nume / 2;
}

void cutLeaf_oneLayer(Graph &g){
    if (calcEffectiveVertexNum(g) <= 2){
        return;
    }
    auto _adj = g.adj;
    for (Vertex u = 0; u < g.numAtom; ++u) {
        if (_adj[u].size() == 1){
            for (auto v : g.adj[u]) {
                g.adj[u].erase(remove(g.adj[u].begin(), g.adj[u].end(), v), g.adj[u].end());
                g.adj[v].erase(remove(g.adj[v].begin(), g.adj[v].end(), u), g.adj[v].end());
            }
            g.status[u] = false;
        }
    }
}

void calcInternalVertexSet(const Graph& g, vector <Vertex> &internal_set){

    internal_set.clear();
    Graph _g = g;
    cutLeaf_oneLayer(_g);
    cutLeaf_oneLayer(_g);

    for (Vertex u = 0; u < g.numAtom; ++u) {
        if (_g.status[u] == true){
            internal_set.push_back(u);
        }
    }
}


void calcCoreVertexSet(const Graph& g, vector <Vertex> &coreset){

    coreset.clear();
    Graph _g = g;
    bool loop = true;
    while (loop) {
        size_t numv = calcEffectiveVertexNum(_g);
        if (numv > 2){
            cutLeaf_oneLayer(_g);
            size_t numv2 = calcEffectiveVertexNum(_g);
            if (numv == numv2){
                loop = false;
            }
        } else {
            loop = false;
        }
    }

    for (Vertex u = 0; u < g.numAtom; ++u) {
        if (_g.status[u] == true){
            coreset.push_back(u);
        }
    }
}

size_t calcDegree(const Graph& g, Vertex u){
    size_t ans = 0;
    for (auto& v : g.adj[u]){
        if (g.status[v]){
            ++ans;
        }
    }
    return ans;
}

size_t calc_cs(const Graph& g){
    vector <Vertex> coreset;
    coreset.clear();
    calcCoreVertexSet(g, coreset);
    return coreset.size();
}

size_t calc_delta(const Graph& g, 
    const vector <Vertex>& coreset,
    Vertex u
){
    size_t ans = 0;
    for (auto& v: g.adj[u]){
        auto itr = find(coreset.begin(), coreset.end(), v);
        if (itr != coreset.end()){
            ++ans;
        }
    }
    return ans;
}

string get_string(input_info& IN_INFO, 
    Graph& g, 
    Vertex u, 
    Vertex parent,
    const vector <Vertex>& internal_set,
    size_t height
){
    string ans = "";
    size_t mul;
    if (parent == u){
        mul = 0;
    } else {
        mul = g.beta[parent][u];
    }
    size_t root_index = IN_INFO.find_index(g.alpha[u]);
    ans += to_string(mul) + to_string(root_index);

    vector <string> subtree_string;
    subtree_string.clear();

    for (auto& v : g.adj[u]){
        auto itr = find(internal_set.begin(), internal_set.end(), v);
        if (v != parent && itr == internal_set.end()){
            string v_string = get_string(IN_INFO, g, v, u, internal_set, height + 1);
            subtree_string.push_back(v_string);
        }
    }
    if (height == 1){    // u is not an internal vertex or a leaf
        size_t subtree_string_size = subtree_string.size();
        for (size_t i = 0; i < IN_INFO.dmax - 1 - subtree_string_size; ++i){
            string tmp_str = "00";
            subtree_string.push_back(tmp_str);
        }
    }
    if (height == 0 && subtree_string.size() == 0){
        for (size_t i = 0; i < IN_INFO.dmax; ++i){
            string tmp_str = "00";
            subtree_string.push_back(tmp_str);
        }
    }
    sort(subtree_string.rbegin(), subtree_string.rend());
    for (auto& tmp : subtree_string){
        ans += tmp;
    }
    return ans;
}

void get_F1(input_info& IN_INFO, const string& inFileName,
    unordered_set <string>& T_W1_F1, unordered_set <string>& T_W2_F1, unordered_set <string>& T_W3_F1
){
    Graph h = read_graph_sdf(inFileName);
    Graph g;
    g.adj.clear();g.alpha.clear();g.beta.clear();g.status.clear();
    H_suppressed_convert(h, g);
    vector <Vertex> internal_set;
    calcInternalVertexSet(g, internal_set);

    for (auto& u : internal_set){
        size_t neightbour_num_int = 0;
        for (auto& v : g.adj[u]){
            auto itr = find(internal_set.begin(), internal_set.end(), v);
            if (itr != internal_set.end()){
                ++neightbour_num_int;
            }
        }
        string fringrtree_string = get_string(IN_INFO, g, u, u, internal_set, 0);
        // cout << fringrtree_string << endl;
        if (neightbour_num_int == 1){
            T_W1_F1.insert(fringrtree_string);
        } else if (neightbour_num_int == 2) {
            T_W2_F1.insert(fringrtree_string);
        } else {
            T_W3_F1.insert(fringrtree_string);
        }
    }
}


void get_partition(
    const Graph& g,
    const vector <Vertex>& coreset, 
    vector <Vertex>& base_vertices,
    vector_2D <Vertex>& base_edges
){
    base_vertices.clear();
    base_edges.clear();

    for (auto& u : coreset){
        if (calc_delta(g, coreset, u) > 2){
            base_vertices.push_back(u);
        } 
    }

    if (base_vertices.size() == 0){
        for (auto& u : coreset){
            auto itr = find(base_vertices.begin(), base_vertices.end(), u);
            if (itr == base_vertices.end()){
                base_vertices.push_back(u);
                if (base_vertices.size() == 1) break;
            }
        }
    }

    unordered_set <Vertex> visited;
    vector_2D <Vertex> cycle;
    for (auto& u : base_vertices){
        visited.insert(u);
        for (auto& v : g.adj[u]){
            auto itr = find(coreset.begin(), coreset.end(), v);
            if (itr != coreset.end()){
                vector <Vertex> tmp;
                tmp.push_back(u);
                tmp.push_back(v);
                auto uu = u;
                auto vv = v;
                while (visited.find(vv) == visited.end() 
                        && 
                    find(base_vertices.begin(), base_vertices.end(), vv) == base_vertices.end()
                ){
                    visited.insert(vv);
                    for (auto& ww : g.adj[vv]){
                        if (ww != uu && find(coreset.begin(), coreset.end(), ww) != coreset.end()){
                            tmp.push_back(ww);
                            uu = vv;
                            vv = ww;
                            break;
                        }
                    }
                }
                // for (size_t i= 0; i< tmp.size(); ++i){
                //     cout <<tmp[i] << " ";
                // }
                // cout << endl;
                if (tmp.size() > 2 || visited.find(v) == visited.end()) {
                    if (tmp[0] == tmp[tmp.size() - 1]){// avoid self-loop
                        cycle.push_back(tmp);
                    } else {
                        base_edges.push_back(tmp);
                    }
                }
            }
        }
    }

    for (auto& tmp : cycle){
        base_vertices.push_back(tmp[1]);
        vector <Vertex> edge_1 = {tmp[0], tmp[1]};
        vector <Vertex> edge_2;
        for (size_t i = 1; i < tmp.size(); ++i){
            edge_2.push_back(tmp[i]);
        }
        base_edges.push_back(edge_1);
        base_edges.push_back(edge_2);
    }
}

void get_partition(
    const Graph& g,
    const vector <Vertex>& coreset, 
    const vector <Vertex>& internal_set,
    vector <Vertex>& base_vertices,
    vector_2D <Vertex>& base_edges
){
    base_vertices.clear();
    base_edges.clear();

    for (auto& u : coreset){
        if (calc_delta(g, internal_set, u) > 2){
            base_vertices.push_back(u);
        } 
    }

    if (base_vertices.size() == 0){
        for (auto& u : coreset){
            auto itr = find(base_vertices.begin(), base_vertices.end(), u);
            if (itr == base_vertices.end()){
                base_vertices.push_back(u);
                if (base_vertices.size() == 1) break;
            }
        }
    }

    unordered_set <Vertex> visited;
    vector_2D <Vertex> cycle;
    for (auto& u : base_vertices){
        visited.insert(u);
        for (auto& v : g.adj[u]){
            auto itr = find(coreset.begin(), coreset.end(), v);
            if (itr != coreset.end()){
                vector <Vertex> tmp;
                tmp.push_back(u);
                tmp.push_back(v);
                auto uu = u;
                auto vv = v;
                while (visited.find(vv) == visited.end() 
                        && 
                    find(base_vertices.begin(), base_vertices.end(), vv) == base_vertices.end()
                ){
                    visited.insert(vv);
                    for (auto& ww : g.adj[vv]){
                        if (ww != uu && find(coreset.begin(), coreset.end(), ww) != coreset.end()){
                            tmp.push_back(ww);
                            uu = vv;
                            vv = ww;
                            break;
                        }
                    }
                }
                // for (size_t i= 0; i< tmp.size(); ++i){
                //     cout <<tmp[i] << " ";
                // }
                // cout << endl;
                if (tmp.size() > 2 || visited.find(v) == visited.end()) {
                    if (tmp[0] == tmp[tmp.size() - 1]){// avoid self-loop
                        cycle.push_back(tmp);
                    } else {
                        base_edges.push_back(tmp);
                    }
                }
            }
        }
    }

    for (auto& tmp : cycle){
        base_vertices.push_back(tmp[1]);
        vector <Vertex> edge_1 = {tmp[0], tmp[1]};
        vector <Vertex> edge_2;
        for (size_t i = 1; i < tmp.size(); ++i){
            edge_2.push_back(tmp[i]);
        }
        base_edges.push_back(edge_1);
        base_edges.push_back(edge_2);
    }
}

void get_partition(
    vector <Vertex>& base_vertices,
    vector_2D <Vertex>& base_edges,
    vector <size_t>& chLB_v,
    vector <size_t>& chUB_v,
    vector <size_t>& chLB_e,
    vector <size_t>& chUB_e,
    vector <bool>& fixed_v,
    vector <bool>& fixed_e,

    string inFileName
){
    ifstream infile;
    try {
        infile.open(inFileName);
    } catch (const exception& e) {
                cerr << "Couldn't open '" << inFileName
                     << "' for reading!" << endl;
                throw(e);
            }

    if (! infile.is_open()) {
        cerr << "Error, infile not initialized!" << endl;
        throw(-1);
    }

    string line;
    stringstream st;
    
    getline(infile, line);
    st << line;
    size_t vc_size;
    st >> vc_size;
    st.clear();

    for (size_t i = 0; i < vc_size; ++i){
        stringstream st1;

        getline(infile, line);
        st1 << line;
        Vertex u;
        st1 >> u;
        base_vertices.push_back(u - 1);
        st1.clear();

        stringstream st2;
        getline(infile, line);
        st2 << line;
        size_t lb, ub;
        st2 >> lb >> ub;
        chLB_v.push_back(lb);
        chUB_v.push_back(ub);

        size_t fixed_tmp;
        if (st2 >> fixed_tmp){
            if (fixed_tmp != 0){
                fixed_v.push_back(true);
            } else {
                fixed_v.push_back(false);
            }
        } else {
            fixed_v.push_back(false);
        }

        st2.clear();
    }

    getline(infile, line);
    st << line;
    size_t ec_size;
    st >> ec_size;
    st.clear();

    for (size_t i = 0; i < ec_size; ++i){
        stringstream st1;

        getline(infile, line);
        st1 << line;

        Vertex u;
        vector <Vertex> tmp;
        while (st1 >> u){
            tmp.push_back(u - 1);
        }
        base_edges.push_back(tmp);
        st1.clear();

        stringstream st2;
        getline(infile, line);
        st2 << line;
        size_t lb, ub;
        st2 >> lb >> ub;
        chLB_e.push_back(lb);
        chUB_e.push_back(ub);

        size_t fixed_tmp;
        if (st2 >> fixed_tmp){
            if (fixed_tmp != 0){
                fixed_e.push_back(true);
            } else {
                fixed_e.push_back(false);
            }
        } else {
            fixed_e.push_back(false);
        }

        st2.clear();
    }

    infile.close();

}

atom_degree get_atom_degree(const input_info& IN_INFO,
    const Graph& g,
    Vertex u
){
    atom_degree a(IN_INFO.find_index(g.alpha[u]), calcDegree(g, u));
    return a;
}

void get_atom_set_tree(input_info& IN_INFO,
    const Graph& g, 
    const vector <Vertex>& coreset, 
    Vertex u, 
    Vertex parent,
    bool root_only = false
){
    if (IN_INFO.find_index(g.alpha[u]) == -1){
        Atom a;
        //if (DEBUG) cout << "g.alpha[u] = " << g.alpha[u] << endl;
        a.NAME = g.alpha[u];
        a.VALENCE = valence_map.at(g.alpha[u]);
        a.MASS = mass_map.at(g.alpha[u]);
        ++IN_INFO.num_kind_atoms;
        IN_INFO.atoms.push_back(a);
    }

    if (root_only) return;

    for (auto& v : g.adj[u]){
        auto itr = find(coreset.begin(), coreset.end(), v);
        if (!g.status[v] || v == parent || itr != coreset.end()) continue;
        get_atom_set_tree(IN_INFO, g, coreset, v, u);
    }
}

void get_descriptors_tree_dfs(input_info& IN_INFO,
    const Graph& g, 
    const vector <Vertex>& coreset, 
    const vector <Vertex>& internal_set,  // set of branch vertices
    Vertex u, 
    Vertex parent,
    size_t height,
    size_t& core_height
    // the  descriptors will store in IN_INFO.rv
    // need to first get the atom set
){
    if (height > core_height) core_height = height;
    if (u != parent){
        atom_degree a = get_atom_degree(IN_INFO, g, parent);
        atom_degree b = get_atom_degree(IN_INFO, g, u);
        size_t m = g.beta[parent][u];
        auto itr1 = find(internal_set.begin(), internal_set.end(), parent);
        auto itr2 = find(internal_set.begin(), internal_set.end(), u);
        if (itr1 != internal_set.end() && itr2 != internal_set.end()){
            ++IN_INFO.rv.ec_in[IN_INFO.atom_deg_map.at(a)][IN_INFO.atom_deg_map.at(b)][m];  
        } else {
            ++IN_INFO.rv.ec_ex[IN_INFO.atom_deg_map.at(a)][IN_INFO.atom_deg_map.at(b)][m];   
        }
    }

    size_t deg = calcDegree(g, u);
    if (deg > 3) IN_INFO.dmax = 4;

    for (auto& v : g.adj[u]){
        auto itr = find(coreset.begin(), coreset.end(), v);
        if (!g.status[v] || v == parent|| itr != coreset.end()) continue;
        get_descriptors_tree_dfs(IN_INFO, g, coreset, internal_set, v, u, height + 1, core_height);
    }
}

void get_descriptors_tree(input_info& IN_INFO,
    const Graph& g, 
    const vector <Vertex>& coreset, 
    const vector <Vertex>& internal_set,  // set of branch vertices
    Vertex u
    // the  descriptors will store in IN_INFO.rv
    // need to first get the atom set
){
    get_atom_set_tree(IN_INFO, g, coreset, u, u);
    sort(IN_INFO.atoms.begin(), IN_INFO.atoms.end());
    IN_INFO.initialize_atom_deg();
    IN_INFO.rv.initialize(IN_INFO.num_kind_atom_degree);
    IN_INFO.initialize_rv_edge_map();

    size_t height = 0;

    get_descriptors_tree_dfs(IN_INFO, g, coreset, internal_set, u, u, 0, height);
    if (height > IN_INFO.k_star){
        IN_INFO.rv.bl += 1;
    }

    if (height > IN_INFO.height){
        IN_INFO.height = height;
    }

    IN_INFO.rv1D = rv2rv1D(IN_INFO.rv, IN_INFO);

    IN_INFO.delta_1 = calc_delta(g, coreset, u);
    IN_INFO.base_vertices = vector <size_t>({u});
    IN_INFO.rs1_index = u;

    size_t d = calcDegree(g, u);
    size_t a = IN_INFO.find_index(g.alpha[u]);
    size_t m = 0;
    for (auto& v : g.adj[u]){
        auto itr = find(coreset.begin(), coreset.end(), v);
        if (itr == coreset.end()){
            m += g.beta[u][v];
        }
    }
    root_status rs(a, m, d);
    IN_INFO.target_rs1 = rs;

}

void get_atom_set_edge(input_info& IN_INFO,
    const Graph& g, 
    const vector <Vertex>& coreset, 
    vector <Vertex>& edge_vertex
){
    for (size_t i = 0; i < edge_vertex.size(); ++i){
        auto u = edge_vertex[i];

        if (i != 0 && i != edge_vertex.size() - 1) {
            get_atom_set_tree(IN_INFO, g, coreset, u, u);
        } else {
            get_atom_set_tree(IN_INFO, g, coreset, u, u, true);
        }
    }
}

void get_descriptors_edge(input_info& IN_INFO,
    const Graph& g, 
    const vector <Vertex>& coreset, 
    const vector <Vertex>& internal_set,  // set of branch vertices
    vector <Vertex>& edge_vertex
    // the  descriptors will store in IN_INFO.rv
    // need to first get the atom set
){
    get_atom_set_edge(IN_INFO, g, coreset, edge_vertex);
    sort(IN_INFO.atoms.begin(), IN_INFO.atoms.end());
    IN_INFO.initialize_atom_deg();
    IN_INFO.rv.initialize(IN_INFO.num_kind_atom_degree);
    IN_INFO.initialize_rv_edge_map();

    for (size_t i = 0; i < edge_vertex.size(); ++i){
        auto u = edge_vertex[i];
        atom_degree a_u = get_atom_degree(IN_INFO, g, u);
        size_t num_neighbour = 0;

        if (i != 0){
            auto v = edge_vertex[i - 1];
            atom_degree a_v = get_atom_degree(IN_INFO, g, v);
            size_t m = g.beta[u][v];
            ++IN_INFO.rv.ec_co[IN_INFO.atom_deg_map.at(a_u)][IN_INFO.atom_deg_map.at(a_v)][m];
            ++IN_INFO.rv.ec_co[IN_INFO.atom_deg_map.at(a_v)][IN_INFO.atom_deg_map.at(a_u)][m];
        }

        if (i != 0 && i != edge_vertex.size() - 1){
            size_t height = 0;

            get_descriptors_tree_dfs(IN_INFO, g, coreset, internal_set, u, u, 0, height);
            if (height > 2){
                IN_INFO.rv.bl += 1;
            }
            if (height > IN_INFO.height){
                IN_INFO.height = height;
            }
        }
    }
    IN_INFO.rv1D = rv2rv1D(IN_INFO.rv, IN_INFO);

    IN_INFO.base_vertices = edge_vertex;

    Vertex u1 = edge_vertex[0];
    size_t d1 = calcDegree(g, u1);
    size_t a1 = IN_INFO.find_index(g.alpha[u1]);
    size_t m1 = g.beta[u1][edge_vertex[1]];
    root_status rs1(a1, m1, d1);
    IN_INFO.target_rs1 = rs1;
    IN_INFO.delta_1 = d1 - 1;
    IN_INFO.rs1_index = u1;

    Vertex u2 = edge_vertex[edge_vertex.size() - 1];
    size_t d2 = calcDegree(g, u2);
    size_t a2 = IN_INFO.find_index(g.alpha[u2]);
    size_t m2 = g.beta[u2][edge_vertex[edge_vertex.size() - 2]];
    root_status rs2(a2, m2, d2);
    IN_INFO.target_rs2 = rs2;
    IN_INFO.delta_2 = d2 - 1;
    IN_INFO.rs2_index = u2;
}

void output_partition(
    const Graph& g,
    const string& outputfilename,
    const vector <Vertex>& base_vertices,
    const vector_2D <Vertex>& base_edges,
    const vector <component>& set_tree_component,
    const vector <component>& set_edge_component
){
    ofstream output(outputfilename, ios::out);

    output << base_vertices.size() << endl;
    for (size_t i = 0; i < base_vertices.size(); ++i){
        string sample_str = g.alpha[base_vertices[i]];
        output << base_vertices[i] + 1 << " # " << sample_str << "\n";
        output << "0 " << set_tree_component[i].IN_INFO.height << " 0\n";
    }

    output << base_edges.size() << endl;
    for (size_t i = 0; i < base_edges.size(); ++i){
        string sample_str = g.alpha[base_edges[i][0]];
        for (size_t j = 1; j < base_edges[i].size(); ++j){
            sample_str += to_string(g.beta[base_edges[i][j - 1]][base_edges[i][j]]) + g.alpha[base_edges[i][j]];
        }
        for (size_t j = 0; j < base_edges[i].size(); ++j){
            output << base_edges[i][j] + 1 << " ";
        }
        output << "# " << sample_str << "\n";
        output << "0 " << set_edge_component[i].IN_INFO.height << " 0\n";
    }

    output.close();
}

void output_SDF(
    size_t& n,
    size_t& m,
    vector_2D <size_t>& graph_adj,
    vector <string>& graph_col,
    string& outputfilename
){

    ofstream output(outputfilename, ios::app);

    output << total_num << "\n";
    output << "BH-cyclic" << "\n";
    output << "BH-cyclic" << "\n";
    output << std::setw(3) << n << std::setw(3) << m << "  0  0  0  0  0  0  0  0999 V2000 " << "\n";

    for (size_t i = 0; i < n; ++i){
        string atom_symbol = graph_col[i];
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

// A function used to generate SDF, especially adding a fringe tree to a vertex
void add_fringe_tree(
    const input_info& IN_INFO,
    const FringeTree& T,
    size_t& graph_ind,
    vector_2D <size_t>& graph_adj,
    vector <string>& graph_col
){
    auto& prt = _prt[IN_INFO.dmax];

    vector <size_t> ind_map(T.seq.size());
    size_t root_ind = graph_ind;
    ind_map[0] = root_ind;
    graph_col[root_ind] = IN_INFO.atoms[T.seq[0].second].NAME;

    for (size_t i = 1; i < T.seq.size(); ++i){
        // if (_DEBUG) cout << "i = " << i << " T.seq.size() = " << T.seq.size() << endl;
        if (T.seq[i].first != 0){
            ++graph_ind;
            graph_col[graph_ind] = IN_INFO.atoms[T.seq[i].second].NAME;
            // if (_DEBUG) cout <<  "i = " << i << " graph_ind = " << graph_ind << " ind_map[prt[i]] = " << ind_map[prt[i]] << " graph_col[graph_ind] = " <<  graph_col[graph_ind] << endl;
            ind_map[i] = graph_ind;
            graph_adj[graph_ind][ind_map[prt[i]]] = T.seq[i].first;
            graph_adj[ind_map[prt[i]]][graph_ind] = T.seq[i].first;
        }
    }
    // if (_DEBUG) cout << "add_fringe_tree end !!!" << endl;
}


void prepare_SDF_branch_tree(
    const input_info& IN_INFO,
    const component_result& comp_graph,
    // vector <int>& st,
    const map_rv& seq,

    size_t& graph_ind,
    vector_2D <size_t>& graph_adj,
    vector <string>& graph_col,

    map <Vertex, size_t>& base_vertex_index,

    size_t root_type,  // 1: W3, 2: W2_core
    bool tree_component
){
    size_t last_n = 0;
    for (size_t h = 0; h < seq.w_ind.size(); ++h){
        // if (_DEBUG) cout << "bt h = " <<  h  << " seq.w_ind.size() = " << seq.w_ind.size() <<  " graph_ind = " << graph_ind << " last_n = " << last_n << endl;
        if (h != 0){
            graph_adj[graph_ind][last_n] = seq.mul[h - 1];
            graph_adj[last_n][graph_ind] = seq.mul[h - 1];   
        }

        last_n = graph_ind;

        if (h == 0){
            if (tree_component){
                base_vertex_index[IN_INFO.rs1_index] = graph_ind;
            }
            if (root_type == 1){  // always first tree-component then edge-component
                // auto& FT_ind = comp_graph.map_FT_W3.at(make_pair(seq.rs[h], seq.w_ind[h]))[st[h]];
                auto& FT_ind = seq.w_ind[h];
                // cout  << " seq.rs[h] = " << seq.rs[h]  << " seq.w_ind[h] = " << seq.w_ind[h] << " FT_ind =  " << FT_ind << endl;
                add_fringe_tree(IN_INFO, comp_graph.All_FT_W3[seq.rs[h]][FT_ind], graph_ind, graph_adj, graph_col);
            } else {
                // auto& FT_ind = comp_graph.map_FT_W2_core.at(make_pair(seq.rs[h], seq.w_ind[h]))[st[h]];
                auto& FT_ind = seq.w_ind[h];
                // cout  << " seq.rs[h] = " << seq.rs[h]  << " seq.w_ind[h] = " << seq.w_ind[h] << " FT_ind =  " << FT_ind << endl;
                add_fringe_tree(IN_INFO, comp_graph.All_FT_W2_core[seq.rs[h]][FT_ind], graph_ind, graph_adj, graph_col);
            }
        } else if (h == seq.w_ind.size() - 1){
            // auto& FT_ind = comp_graph.map_FT_W1.at(make_pair(seq.rs[h], seq.w_ind[h]))[st[h]];
            auto& FT_ind = seq.w_ind[h];
            // cout  << " seq.rs[h] = " << seq.rs[h]  << " seq.w_ind[h] = " << seq.w_ind[h] << " FT_ind =  " << FT_ind << endl;
            add_fringe_tree(IN_INFO, comp_graph.All_FT_W1[seq.rs[h]][FT_ind], graph_ind, graph_adj, graph_col);
        } else {
            // auto& FT_ind = comp_graph.map_FT_W2.at(make_pair(seq.rs[h], seq.w_ind[h]))[st[h]];
            auto& FT_ind = seq.w_ind[h];
            // cout  << " seq.rs[h] = " << seq.rs[h]  << " seq.w_ind[h] = " << seq.w_ind[h] << " FT_ind =  " << FT_ind << endl;
            add_fringe_tree(IN_INFO, comp_graph.All_FT_W2[seq.rs[h]][FT_ind], graph_ind, graph_adj, graph_col);
        }

        ++graph_ind;
        // if (_DEBUG) cout << "bt h = " <<  h  << " seq.w_ind.size() = " << seq.w_ind.size() <<  " graph_ind = " << graph_ind << " last_n = " << last_n << " end!!!" << endl;
        
    }
}

void prepare_SDF_component(
    const input_info& IN_INFO,
    const component_result& comp_graph,
    const map_tree& seq,

    size_t& graph_ind,
    vector_2D <size_t>& graph_adj,
    vector <string>& graph_col,

    map <Vertex, size_t>& base_vertex_index,

    bool tree_component
){
    if (tree_component){
        auto& seq_rv = seq.trees[0];
        // vector <int> st(seq_rv.w_ind.size() + 1, 0);
        size_t root_type;
        if (IN_INFO.height > IN_INFO.k_star){
            root_type = 1;
        } else {
            root_type = 2;
        }
        // cout << "seq_rv.w_ind.size() = " << seq_rv.w_ind.size() << endl;
        prepare_SDF_branch_tree(IN_INFO, comp_graph, seq_rv, graph_ind, graph_adj, graph_col, base_vertex_index, root_type, tree_component);
    } else {
        size_t last_n = graph_ind;
        for (size_t h = 0; h < seq.trees.size(); ++h){
            // if (_DEBUG) cout << "co h = " << h << " seq.trees.size() = " << seq.trees.size() << " graph_ind = " << graph_ind << " last_n = " << last_n << endl;
            if (h != 0){
                size_t ind_tmp_1, ind_tmp_2;
                if (h == 1 && h == seq.trees.size() - 1){
                    ind_tmp_1 = base_vertex_index.at(IN_INFO.rs1_index);
                    ind_tmp_2 = base_vertex_index.at(IN_INFO.rs2_index); 
                } else if (h == 1){
                    ind_tmp_1 = base_vertex_index.at(IN_INFO.rs1_index);
                    ind_tmp_2 = graph_ind;
                } else if (h == seq.trees.size() - 1){
                    ind_tmp_1 = last_n;
                    ind_tmp_2 = base_vertex_index.at(IN_INFO.rs2_index); 
                } else {
                    ind_tmp_1 = last_n;
                    ind_tmp_2 = graph_ind;
                }
                // if (_DEBUG) cout << "ind_tmp_1 = " << ind_tmp_1 << " ind_tmp_2 = " << ind_tmp_2 << endl;
                graph_adj[ind_tmp_1][ind_tmp_2] = seq.mul[h - 1];
                graph_adj[ind_tmp_2][ind_tmp_1] = seq.mul[h - 1]; 
            }

            last_n = graph_ind;

            if (h != 0 && h != seq.trees.size() - 1){
                auto& seq_rv = seq.trees[h];
                // vector <int> st(seq_rv.w_ind.size() + 1, 0);
                size_t root_type;
                // cout << "seq_rv.w_ind.size() = " << seq_rv.w_ind.size() << endl;
                if (seq_rv.w_ind.size() > 1){
                    root_type = 1;
                } else{
                    root_type = 2;
                }

                prepare_SDF_branch_tree(IN_INFO, comp_graph, seq_rv, graph_ind, graph_adj, graph_col, base_vertex_index, root_type, tree_component);
            }
            // if (_DEBUG) cout << "co h = " << h << " seq.trees.size() = " << seq.trees.size() << " graph_ind = " << graph_ind << " last_n = " << last_n << " end!!!!!!" << endl;
            
        }
    }
}

void prepare_SDF_branch_tree_from_g(
    const Graph& g,
    vector <bool>& visited,
    Vertex u,

    size_t& graph_ind,
    vector_2D <size_t>& graph_adj,
    vector <string>& graph_col,

    map <Vertex, size_t>& base_vertex_index,

    bool tree_component
){
    stack <Vertex> _stack;
    _stack.push(u);

    map <Vertex, size_t> g_prt_map_ind;
    map <Vertex, Vertex> g_prt_map;

    base_vertex_index[u] = graph_ind;

    while (!_stack.empty()){
        Vertex _u = _stack.top();
        _stack.pop();
        
        graph_col[graph_ind] = g.alpha[_u];

        if (_u != u){
            Vertex prt_ind = g_prt_map_ind.at(_u);
            Vertex _v = g_prt_map.at(_u);
            graph_adj[graph_ind][prt_ind] = g.beta[_u][_v];
            graph_adj[prt_ind][graph_ind] = g.beta[_v][_u];
        }

        for (auto& _v : g.adj[_u]){
            if (visited[_v]){
                continue;
            }
            _stack.push(_v);

            g_prt_map[_v] = _u;
            g_prt_map_ind[_v] = graph_ind;
            visited[_v] = true;
        }

        ++graph_ind;
    }
}

void prepare_SDF_component_from_g(
    const Graph& g,
    const vector <Vertex>& coreset,
    const vector <Vertex>& base_vertices,

    size_t& graph_ind,
    vector_2D <size_t>& graph_adj,
    vector <string>& graph_col,

    map <Vertex, size_t>& base_vertex_index,

    bool tree_component
){
    vector <bool> visited(g.alpha.size(), false);
    for (auto& v : coreset){
        visited[v] = true;
    }

    if (tree_component){
        Vertex u = base_vertices[0];
        prepare_SDF_branch_tree_from_g(g, visited, u, graph_ind, graph_adj, graph_col, base_vertex_index, tree_component);
    } else {
        size_t size_tmp = base_vertices.size();
        size_t last_n = graph_ind;
        
        for (size_t h = 0; h < size_tmp; ++h){
            auto& u = base_vertices[h];
            
            if (h != 0){
                size_t ind_tmp_1, ind_tmp_2;
                if (h == 1 && h == size_tmp - 1){
                    ind_tmp_1 = base_vertex_index.at(base_vertices[0]);
                    ind_tmp_2 = base_vertex_index.at(base_vertices[size_tmp - 1]); 
                } else if (h == 1){
                    ind_tmp_1 = base_vertex_index.at(base_vertices[0]);
                    ind_tmp_2 = graph_ind;
                } else if (h == size_tmp - 1){
                    ind_tmp_1 = last_n;
                    ind_tmp_2 = base_vertex_index.at(base_vertices[size_tmp - 1]); 
                } else {
                    ind_tmp_1 = last_n;
                    ind_tmp_2 = graph_ind;
                }
                // if (_DEBUG) cout << "ind_tmp_1 = " << ind_tmp_1 << " ind_tmp_2 = " << ind_tmp_2 << endl;
                graph_adj[ind_tmp_1][ind_tmp_2] = g.beta[base_vertices[h - 1]][base_vertices[h]];
                graph_adj[ind_tmp_2][ind_tmp_1] = g.beta[base_vertices[h - 1]][base_vertices[h]]; 
            }

            last_n = graph_ind;
            if (h != 0 && h != size_tmp - 1){
                prepare_SDF_branch_tree_from_g(g, visited, u, graph_ind, graph_adj, graph_col, base_vertex_index, tree_component);
            }
        }
    }
}

void prepare_SDF(
    const Graph& g,
    const vector <Vertex>& coreset,

    const vector <component>& set_tree_component,
    const vector <component>& set_edge_component,
    size_t& g_n,
    size_t& g_m,
    const size_t& core_height_ind,

    vector <int>& st,
    string& outputfilename
){
    vector_2D <size_t> graph_adj(g_n, vector <size_t> (g_n, 0));
    vector <string> graph_col(g_n, "");
    size_t graph_ind = 0;

    map <Vertex, size_t> base_vertex_index;
    base_vertex_index.clear();

    // size_t height = 0;
    // for (size_t i = 0; i < set_tree_component.size(); ++i){
    //     auto& comp_graph = set_tree_component[i].comp_graph;
    //     if (comp_graph.all_seq[st[i]].height > height){
    //         height = comp_graph.all_seq[st[i]].height;
    //     }
    // }
    // for (size_t i = 0; i < set_edge_component.size(); ++i){
    //     auto& comp_graph = set_edge_component[i].comp_graph;
    //     if (comp_graph.all_seq[st[i]].height > height){
    //         height = comp_graph.all_seq[st[i]].height;
    //     }
    // }
    // if (height != core_height) return;

    // for (auto& _h : st){
    //     cout  << _h << " ";
    // }
    // cout <<  endl;

    for (size_t i = 0; i < set_tree_component.size(); ++i){
        auto& IN_INFO = set_tree_component[i].IN_INFO;
        auto& comp_graph = set_tree_component[i].comp_graph;

        if (set_tree_component[i].no_feasible_pair){
            prepare_SDF_component_from_g(g, coreset, IN_INFO.base_vertices, graph_ind, graph_adj, graph_col, base_vertex_index, true);
            continue;
        }

        if (i == core_height_ind){
            prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_h[st[i]], graph_ind, graph_adj, graph_col, base_vertex_index, true);
        } else if (i < core_height_ind){
            if (st[i] < comp_graph.all_seq_h.size()){
                prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_h[st[i]], graph_ind, graph_adj, graph_col, base_vertex_index, true);
            } else {
                prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq[st[i] - comp_graph.all_seq_h.size()], graph_ind, graph_adj, graph_col, base_vertex_index, true);
            }
        } else {
            prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq[st[i]], graph_ind, graph_adj, graph_col, base_vertex_index, true);
        }
    }

    size_t tree_comp_num = set_tree_component.size();
    for (size_t i = 0; i < set_edge_component.size(); ++i){
        auto& IN_INFO = set_edge_component[i].IN_INFO;
        auto& comp_graph = set_edge_component[i].comp_graph;

        if (set_edge_component[i].no_feasible_pair){
            prepare_SDF_component_from_g(g, coreset, IN_INFO.base_vertices, graph_ind, graph_adj, graph_col, base_vertex_index, false);
            continue;
        }

        if (i + tree_comp_num == core_height_ind){
            prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_h[st[i + tree_comp_num]], graph_ind, graph_adj, graph_col, base_vertex_index, false);    
        } else if (i + tree_comp_num < core_height_ind) {
            if (st[i + tree_comp_num] < comp_graph.all_seq_h.size()){
                prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_h[st[i + tree_comp_num]], graph_ind, graph_adj, graph_col, base_vertex_index, false);
            } else {
                prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq[st[i + tree_comp_num] - comp_graph.all_seq_h.size()], graph_ind, graph_adj, graph_col, base_vertex_index, false);
            }
        } else {
            prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq[st[i + tree_comp_num]], graph_ind, graph_adj, graph_col, base_vertex_index, false);
        }    
    }
    ++total_num;
    if (total_num >= num_limit) _stop_gen = true;
    output_SDF(g_n, g_m, graph_adj, graph_col, outputfilename);
    // if (_DEBUG) cout << "total_num = " << total_num << endl;
}


void generate(
    const Graph& g,
    const vector <Vertex>& coreset,

    const vector <component>& set_tree_component,
    const vector <component>& set_edge_component,
    size_t& g_n,
    size_t& g_m,
    const size_t& core_height_ind,

    size_t& possible_num,
    size_t& lower_bound,

    string& outputfilename
){

    size_t n_comp = set_tree_component.size() + set_edge_component.size();
    vector <int> st(n_comp + 1, -1);
    vector <size_t> st_max(n_comp + 1, 0);
    size_t possible_num_tmp = 1;
    size_t lower_bound_tmp = 1;

    for (size_t i = 0; i < n_comp; ++i){
        // cout << "i = " << i << endl;
        if (i < set_tree_component.size()){
            if (i == core_height_ind){
                if (set_tree_component[i].no_feasible_pair && set_tree_component[i].with_core_height){
                    st_max[i] = 1;
                } else {
                    st_max[i] = set_tree_component[i].comp_graph.all_seq_h.size();
                }
                // cout << "num_h = " << set_tree_component[i].comp_graph.num_h << endl;
                lower_bound_tmp *= set_tree_component[i].comp_graph.num_h;
            } else if (i < core_height_ind){
                if (set_tree_component[i].no_feasible_pair){
                    st_max[i] = 1;
                } else {
                    st_max[i] = set_tree_component[i].comp_graph.all_seq_h.size() + set_tree_component[i].comp_graph.all_seq.size();
                }
                // cout << "num_h + num = " << set_tree_component[i].comp_graph.num_h + set_tree_component[i].comp_graph.num << endl;
                lower_bound_tmp *= (set_tree_component[i].comp_graph.num_h + set_tree_component[i].comp_graph.num);
            } else {
                if (set_tree_component[i].no_feasible_pair && !set_tree_component[i].with_core_height){
                    st_max[i] = 1;
                } else {
                    st_max[i] = set_tree_component[i].comp_graph.all_seq.size();
                }
                // cout << "num = " << set_tree_component[i].comp_graph.num << endl;
                lower_bound_tmp *= set_tree_component[i].comp_graph.num;
            }
        } else {
            if (i == core_height_ind){
                if (set_edge_component[i - set_tree_component.size()].no_feasible_pair && set_edge_component[i - set_tree_component.size()].with_core_height){
                    st_max[i] = 1;
                } else {
                    st_max[i] = set_edge_component[i - set_tree_component.size()].comp_graph.all_seq_h.size();
                }
                // cout << "num_h = " << set_edge_component[i - set_tree_component.size()].comp_graph.num_h << endl;
                lower_bound_tmp *= set_edge_component[i - set_tree_component.size()].comp_graph.num_h;
            } else if (i < core_height_ind){
                if (set_edge_component[i - set_tree_component.size()].no_feasible_pair){
                    st_max[i] = 1;
                } else {
                    st_max[i] = set_edge_component[i - set_tree_component.size()].comp_graph.all_seq_h.size() + set_edge_component[i - set_tree_component.size()].comp_graph.all_seq.size();
                }
                // cout << "num_h + num = " << set_edge_component[i - set_tree_component.size()].comp_graph.num_h + set_edge_component[i - set_tree_component.size()].comp_graph.num << endl;
                lower_bound_tmp *= (set_edge_component[i - set_tree_component.size()].comp_graph.num_h + set_edge_component[i - set_tree_component.size()].comp_graph.num);
            } else {
                if (set_edge_component[i - set_tree_component.size()].no_feasible_pair && !set_edge_component[i - set_tree_component.size()].with_core_height){
                    st_max[i] = 1;
                } else {
                    st_max[i] = set_edge_component[i - set_tree_component.size()].comp_graph.all_seq.size();
                }
                // cout << "num = " << set_edge_component[i - set_tree_component.size()].comp_graph.num << endl;
                lower_bound_tmp *= set_edge_component[i - set_tree_component.size()].comp_graph.num;
            }
        }
        if (st_max[i] == 0) return;
        possible_num_tmp *= st_max[i];
    }    
    possible_num += possible_num_tmp;
    lower_bound += lower_bound_tmp;
    // cout << "lower_bound_tmp = " << lower_bound_tmp << endl;

    int h = 0;

    while (h >= 0){
        if (_stop_gen) return;
        // if (_DEBUG) {
        //     cout << "h = "  << h<< endl;
        //     for (size_t i = 0; i < st.size(); ++i){
        //         cout << st[i] << " " ;
        //     }
        //     cout << endl;
        // }
        
        if (h == set_tree_component.size() + set_edge_component.size()){    
            prepare_SDF(g, coreset, set_tree_component, set_edge_component, g_n, g_m, core_height_ind, st, outputfilename);
            --h;
            // if (_DEBUG) cout << "h = " << h<< " end!!!!!"<< endl;
        } else {
            if (st[h] + 1 < st_max[h]){
                ++st[h];
                ++h;
                st[h] = -1;
            } else {
                --h;
            }
        }
    }
}

// A function to generate SDF format files
void Gen_Graph(    
    const Graph& g,
    const vector <Vertex>& coreset,

    const vector <component>& set_tree_component,
    const vector <component>& set_edge_component,
    size_t& g_n,
    size_t& g_m,

    size_t& possible_num,
    size_t& lower_bound,

    string& outputfilename
){
    ofstream output(outputfilename, ios::out);
    output.close();

    size_t n_comp = set_tree_component.size() + set_edge_component.size();
    possible_num = 0;
    lower_bound = 0;

    for (size_t i = 0; i < n_comp; ++i){
        generate(g, coreset, set_tree_component, set_edge_component, g_n, g_m, i, possible_num, lower_bound, outputfilename);
    }
}

#endif
