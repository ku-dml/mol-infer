// This file contains the code of the necessary functions about chemical graph.

#ifndef __CHEMICAL_GRAPH_HPP__INCLUDED
#define __CHEMICAL_GRAPH_HPP__INCLUDED

#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>
#include <queue>

#include "cross_timer.h"
#include "data_structures.hpp"
// #include "mem_osx.h"

using namespace std;

bool _stop_gen = false;
size_t num_limit = 100;  //  limit of the number of graph to output
size_t seq_limit = 2;    // limit of sequences generated per component
size_t global_time_limit = 3600;
// double memory_limit = 20000000;
size_t global_path_limit = 0;

size_t total_num = 0;

double time_count_path;
double _start_time;
double _start_memory;

bool _MEMORY_OUT = false;
bool _GLOBAL_TIME_OUT = false;

mt19937 engine(1);

typedef size_t Vertex;
class Graph {
public:
  string CID;                   // CID
  size_t numAtom, numBond;      // number of vertices and edges
  vector<vector<Vertex>> adj;   // adjacent vertices
  vector<string> alpha;         // atom types (C, H, O, ...)
  vector<vector<size_t>> beta;  // bond types (1, 2, ...)
  vector<bool> status;          // used in calculating branch height
  // new, 0430
  vector <int> charge;

  Graph(){}

  Graph(const Graph& graph, vector <Vertex>& efflist){
      CID = graph.CID;
      numAtom = efflist.size();
      initialize(numAtom);
      for (Vertex u = 0; u < numAtom; u++){
        alpha[u] = graph.alpha[efflist[u]];
        charge[u] = graph.charge[efflist[u]];
        status[u] = true;
        for (auto v : graph.adj[efflist[u]]){
            auto itr = find(efflist.begin(), efflist.end(), v);
            adj[u].push_back(distance(efflist.begin(), itr));
        }
        for (Vertex v = 0; v < numAtom; v++){
            beta[u][v] = graph.beta[efflist[u]][efflist[v]];
        }
      }
  }

  vector <Vertex> calcEffectiveAtom(){
        vector <Vertex> efflist = {};
        for (Vertex u = 0; u < numAtom; u++){
            if (status[u]){
                efflist.push_back(u);
            }
        }
        return efflist;
    }

  void initialize(int num){
      numAtom = num;
        alpha.resize(num);
        status.resize(num);
    charge.resize(num);
        
        adj = {};
        vector <Vertex> tmp;
        for (Vertex u = 0; u < num; u++){
            adj.push_back(tmp);
            alpha[u] = "";
      charge[u] = 0;
            status[u] = false;
        }
        beta.resize(num);
        for (Vertex u = 0; u < num; u++){
            beta[u].resize(num);
            for (Vertex v = 0; v < num; v++){
                beta[u][v] = 0;
            }
        }   
    }
};


int calcEffectiveVertexNum(const Graph& g) {
  int numv = 0;
  for (Vertex u = 0; u < g.numAtom; ++u) {
    if (g.status[u] == true) {
      numv++;
    }
  }

  return numv;
}

int calcNumOfEdges(const Graph& g) {
  int nume = 0;
  for (Vertex u = 0; u < g.numAtom; u++) {
    nume += g.adj[u].size();
  }
  return nume / 2;
}

void cutLeaf_oneLayer(Graph& g) {
  if (calcEffectiveVertexNum(g) <= 2) {
    return;
  }
  auto _adj = g.adj;
  for (Vertex u = 0; u < g.numAtom; ++u) {
    if (_adj[u].size() == 1) {
      for (auto v : g.adj[u]) {
        g.adj[u].erase(remove(g.adj[u].begin(), g.adj[u].end(), v),
                       g.adj[u].end());
        g.adj[v].erase(remove(g.adj[v].begin(), g.adj[v].end(), u),
                       g.adj[v].end());
      }
      g.status[u] = false;
    }
  }
}

void calcInternalVertexSet(const Graph& g, vector<Vertex>& internal_set) {
  internal_set.clear();
  Graph _g = g;
  cutLeaf_oneLayer(_g);
  cutLeaf_oneLayer(_g);

  for (Vertex u = 0; u < g.numAtom; ++u) {
    if (_g.status[u] == true) {
      internal_set.push_back(u);
    }
  }
}

void calcCoreVertexSet(const Graph& g, vector<Vertex>& coreset) {
  coreset.clear();
  Graph _g = g;
  bool loop = true;
  while (loop) {
    size_t numv = calcEffectiveVertexNum(_g);
    if (numv > 2) {
      cutLeaf_oneLayer(_g);
      size_t numv2 = calcEffectiveVertexNum(_g);
      if (numv == numv2) {
        loop = false;
      }
    } else {
      loop = false;
    }
  }

  for (Vertex u = 0; u < g.numAtom; ++u) {
    if (_g.status[u] == true) {
      coreset.push_back(u);
    }
  }
}

size_t calcDegree(const Graph& g, Vertex u) {
  size_t ans = 0;
  for (auto& v : g.adj[u]) {
    if (g.status[v]) {
      ++ans;
    }
  }
  return ans;
}

size_t calc_cs(const Graph& g) {
  vector<Vertex> coreset;
  coreset.clear();
  calcCoreVertexSet(g, coreset);
  return coreset.size();
}

size_t calc_delta(const Graph& g, const vector<Vertex>& coreset, Vertex u) {
  size_t ans = 0;
  for (auto& v : g.adj[u]) {
    auto itr = find(coreset.begin(), coreset.end(), v);
    if (itr != coreset.end()) {
      ++ans;
    }
  }
  return ans;
}

/////////////////
// classes and functions used to compute TS

void searchBranch(Vertex u, vector <Vertex>& forbiddenset, const Graph& graph, Graph& branch){
    for (auto v : graph.adj[u]){
        if (!graph.status[v]) continue;
        auto itr = find(forbiddenset.begin(), forbiddenset.end(), v);
        if (itr == forbiddenset.end()){
            branch.status[v] = true;
            branch.adj[u].push_back(v);
            branch.adj[v].push_back(u);
            branch.beta[u][v] = graph.beta[u][v];
            branch.beta[v][u] = graph.beta[v][u];
            branch.alpha[v] = graph.alpha[v];
            branch.charge[v] = graph.charge[v];
            auto forbiddenset_new = forbiddenset;
            forbiddenset_new.push_back(v);
            searchBranch(v, forbiddenset_new, graph, branch);
        }
    }
}

class RootedTree{

public:
    Vertex root;
    Graph graph;
    vector <int> depth;

    RootedTree(){
      root = 0;
      depth.clear();
    }

    RootedTree(Vertex _root, const Graph& _graph){
        root = _root;
        graph = _graph;
        depth = vector <int>(graph.numAtom, -1);
        if (graph.numAtom > 0){
            calcdepth();
        }
    }

    RootedTree(const RootedTree& _rt){
        root = _rt.root;
        graph = _rt.graph;
        depth = _rt.depth;
    }

    void calcdepth(){
        depth[root] = 0;
        queue <Vertex> q;
        q.push(root);
        while (!q.empty()){
            Vertex u = q.front();
            q.pop();
            for (auto v : graph.adj[u]){
                if (graph.status[v] && depth[v] == -1){
                    depth[v] = depth[u] + 1;
                    q.push(v);
                }
            }
        }
    }

    void buildFromCoreVertex(Vertex _root, vector <Vertex> CoreSet, const Graph& _graph){
        root = _root;
        graph = _graph;
        depth = vector <int>(graph.numAtom, -1);
        graph.initialize(_graph.numAtom);
        graph.status[root] = true;
        graph.alpha[root] = _graph.alpha[root];
        graph.charge[root] = _graph.charge[root];
        searchBranch(root, CoreSet, _graph, graph);
        calcdepth();
    }
};

RootedTree RemoveNullVertex(RootedTree& RT){
    int EAN = calcEffectiveVertexNum(RT.graph);
    vector <Vertex> efflist = RT.graph.calcEffectiveAtom();
    Graph new_graph(RT.graph, efflist);
    int new_root = -1;
    for (int i = 0; i < efflist.size(); ++i){
        if (efflist[i] == RT.root){
            new_root = i;
            break;
        }
    }
    RootedTree new_RT(new_root, new_graph); 
    return new_RT;
}

class TopologyGraph{
    
public:
    string name;
    int numVertex;

    vector <vector <Vertex>> adj;
    vector <bool> status;
    vector <int> charge;

    TopologyGraph(){
        name = "";
        numVertex = 0;
        adj = {};
        status = {};
        charge = {};
    }

    TopologyGraph(const TopologyGraph& Tgraph){
        name = Tgraph.name;
        numVertex = Tgraph.numVertex;
        adj = Tgraph.adj;
        status = Tgraph.status;
        charge = Tgraph.charge;
    }

    TopologyGraph(const Graph& graph){
        name = graph.CID;
        numVertex = graph.numAtom;
        adj = graph.adj;
        status = graph.status;
        charge = graph.charge;
    }

    TopologyGraph(const TopologyGraph& Tgraph, vector <Vertex>& efflist){
        numVertex = efflist.size();
        initialize(numVertex);
        name = Tgraph.name;
        for (Vertex u = 0; u < numVertex; u++){
            status[u] = true;
            charge[u] = Tgraph.charge[efflist[u]];
            for (auto v : Tgraph.adj[efflist[u]]){
                auto itr = find(efflist.begin(), efflist.end(), v);
                adj[u].push_back(distance(efflist.begin(), itr));
            }
        }
    }

    void initialize(int num){
        name = "";
        numVertex = num;
        status.resize(num); 
        adj = {};
        charge.resize(num);
        vector <Vertex> tmp;
        for (Vertex u = 0; u < num; u++){
            adj.push_back(tmp);
            status[u] = false;
            charge[u] = 0;
        }   
    }

    vector <int> calcAdjofVertex(Vertex u){
        vector <int> adjofvertex(numVertex, 0);
        for (auto v : adj[u]){
            if (status[v]) {
                adjofvertex[v]++;
            }
        }
        return adjofvertex;
    }


    int degree(Vertex u){
        vector <int> adju = calcAdjofVertex(u);
        int deg = 0;
        for (Vertex v = 0; v < numVertex; v++){
            deg += adju[v];
            if (v == u){
                deg += adju[v];
            }
        }
        return deg;
    }
};

class RootedGraph{
public:
    TopologyGraph noColorGraph;

    vector <int> depth;

    Vertex root;
    int numVertexAll;
    int numVertexNow;

    RootedGraph(){}

    RootedGraph(const RootedGraph& _RG){
        noColorGraph = _RG.noColorGraph;
        root = _RG.root;
        numVertexAll = _RG.numVertexAll;
        numVertexNow = _RG.numVertexNow;
        depth = _RG.depth;

    }

    RootedGraph(const RootedTree& _RT){
        noColorGraph = TopologyGraph(_RT.graph);
        root = _RT.root;
        numVertexAll = _RT.graph.numAtom;
        numVertexNow = _RT.graph.numAtom;
        depth = _RT.depth;
    }

    RootedGraph(int n, Vertex u = 0){
        //noColorGraph = TopologyGraph();
        numVertexAll = n;
        numVertexNow = 1;
        noColorGraph.initialize(n);
        depth = vector <int> (n, -1);
        root = u;
        noColorGraph.status[u] = true;
        depth[u] = 0;
    }
};


class RootedMultiTree{
private:

    Color _edge_color(Vertex u){   //edge u to parent(u)
        // return max_edge_color * mul[u] + edge_color[u];
        return mul[u];
    }

    void _signature_rmt(TreeSignature& TS, Vertex root, int depth){
        int child_num = RG.noColorGraph.degree(root) - 1;
        if (root == r) ++child_num;

        if (child_num <= 0){
            TS.delta.clear();
            TS.delta.push_back(label[root]);
            TS.delta.push_back(depth);
            TS.mu.clear();
            TS.chg.clear();
            TS.chg.push_back(RG.noColorGraph.charge[root]);
            return;
        }

        vector <TreeSignature> tmp_TS(child_num);

        int ind = 0;
        for (auto& u : RG.noColorGraph.adj[root]){
            if (RG.noColorGraph.status[u] && u == parent[root]) continue;

            tmp_TS[ind].delta.clear();
            tmp_TS[ind].mu.clear();
            tmp_TS[ind].chg.clear();

            _signature_rmt(tmp_TS[ind], u, depth + 1);
            tmp_TS[ind].mu.insert(tmp_TS[ind].mu.begin(), _edge_color(u));
            ++ind;
        }

        sort(tmp_TS.begin(), tmp_TS.end(), TS_cmp);

        TS.delta.push_back(label[root]);
        TS.delta.push_back(depth);
        TS.chg.push_back(RG.noColorGraph.charge[root]);

        for (int i = 0; i < child_num; ++i){
            TS.extend(tmp_TS[i]);
        }

        return;
    }

public:
    RootedGraph RG;
    int n;
    int r;
    int max_edge_color;
    vector <Color> label;
    vector <Vertex> parent;
    vector <int> mul;   // the mul of the edge i to parent(i)
    vector <Color> edge_color;

    RootedMultiTree(RootedGraph& _RG){
        RG = _RG;
        n = _RG.numVertexNow;
        r = _RG.root;
        max_edge_color = 1;
        mul = vector <int>(n, 1);
        label = vector <Color>(n, 1);
        edge_color = vector <Color>(n, 1);
        parent = vector <Vertex>(n, -1);
        
        for (Vertex u = 0; u < n; u++){
            
            if (u == _RG.root){
                parent[u] = -1;
                mul[u] = 0;
            } else {
                parent[u] = _RG.noColorGraph.adj[u][0];
            }
        }
    }

    RootedMultiTree(RootedGraph& _RG, RootedTree& _RT){
        RG = _RG;
        n = _RG.numVertexNow;
        r = _RG.root;
        max_edge_color = 1;
        mul = vector <int>(n, 1);
        label = vector <Color>(n, 1);
        edge_color = vector <Color>(n, 1);
        parent = vector <Vertex>(n, -1);
        
        for (Vertex u = 0; u < n; u++){
            
            if (u == _RG.root){
                parent[u] = -1;
                mul[u] = 0;
            } else {
                parent[u] = _RG.noColorGraph.adj[u][0];
                mul[u] = _RT.graph.beta[u][parent[u]];
            }
        }
    }

    TreeSignature getSignature(){
        TreeSignature TS;
        _signature_rmt(TS, r, 0);

        return TS;
    }
};
/////////////////////////////////////////

vector <RootedTree> getFringeTrees(Graph& graph){
    Graph _graph = graph;
    vector <Vertex> BranchTree;
    calcInternalVertexSet(_graph, BranchTree);
    vector <RootedTree> FringeTrees = {};
    for (auto& u : BranchTree){
        RootedTree ct;
        ct.buildFromCoreVertex(u, BranchTree, graph);
        // if (calcEffectiveVertexNum(ct.graph) > 1){
        FringeTrees.push_back(ct);
        // }
    }
    return FringeTrees;
}

Graph read_graph_sdf(const string& inFileName) {
  Graph read_graphs;  // return value
  // open the file
  ifstream infile;
  try {
    infile.open(inFileName);
  } catch (const exception& e) {
    cerr << "Couldn't open '" << inFileName << "' for reading!" << endl;
    throw(e);
  }

  if (!infile.is_open()) {
    cerr << "Error, infile not initialized!" << endl;
    throw(-1);
  }
  // read the file
  string line;
  int flag = 0;

  while (getline(infile, line)) {
    while (flag == 0) {
      read_graphs.CID = line;
      // Read and discard the next two lines of the file
      getline(infile, line);
      getline(infile, line);

      int n, m;
      stringstream st;
      char* threeChars = new char[3];
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

      // new, 0430
      read_graphs.charge.resize(n);

      // a local variable to store the atom symbol
      string alpha;
      // read atom data
      for (int i = 0; i < n; ++i) {
        getline(infile, line);
        stringstream st1;
        // temporary storage for the 3D data from a line in the sdf file
        char* tdData = new char[31];
        st1 << line;
        // read 3D coordinate data in a temporary string storage
        st1.read(tdData, 31);  // read the 3D information from the line
        st1 >> alpha;
        read_graphs.alpha[i] = alpha;
        read_graphs.status[i] = true;

        // new, 0430
        read_graphs.charge[i] = 0;

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
        vector<Vertex> tmp;
        read_graphs.adj.push_back(tmp);
      }

      // read bond data
      for (int i = 0; i < m; ++i) {
        // local variable for vertex ids
        int v1, v2, mul;
        getline(infile, line);
        stringstream st2;
        char* vChars = new char[3];
        st2 << line;
        st2.read(vChars, 3);
        v1 = atoi(vChars);  // id, 1-n of vertex
        st2.read(vChars, 3);
        v2 = atoi(vChars);  // id, 1-n of vertex
        st2.read(vChars, 3);
        mul = atoi(vChars);  // v1-v2 multiplicity
        delete[] vChars;

        // account for off-by-one in the indexing
        read_graphs.adj[v1 - 1].push_back(v2 - 1);
        read_graphs.adj[v2 - 1].push_back(v1 - 1);
        read_graphs.beta[v1 - 1][v2 - 1] = mul;
        read_graphs.beta[v2 - 1][v1 - 1] = mul;
      }
      flag = 1;
    }

    // new, 0430
    if (line.substr(0, 6) == "M  CHG"){
// read the charge and store it
      string charge_value = line.substr(6);
      string stt = charge_value.substr(0, 3);
      stt.erase(remove(stt.begin(), stt.end(), ' '), stt.end());
      int n_chg = stoi(stt);
      for (int i = 0; i < n_chg; i++) {
        int ind_atom, atom_chg;
        stt = charge_value.substr(3 + i*8, 4);
        stt.erase(remove(stt.begin(), stt.end(), ' '), stt.end());
        ind_atom = stoi(stt);
        stt = charge_value.substr(3 + 4 + i * 8, 4);
        stt.erase(remove(stt.begin(), stt.end(), ' '), stt.end());
        atom_chg = stoi(stt);
        read_graphs.charge[ind_atom - 1] = atom_chg;
      }
    }

    // The end of the graph information in the file is marked with a "$$$$"
    if (line == "$$$$" && flag == 1) {
      break;
    }
  }
  return read_graphs;
}

// new, 0430
// rename the atoms in the graph w.r.t charge and multi-valence
void change_atom_name(Graph& h){
  for (Vertex u = 0; u < h.numAtom; ++u){
    int beta_u = 0;
    for (Vertex v : h.adj[u]){
      beta_u += h.beta[u][v];
    }
    auto alpha_org = h.alpha[u];
    h.alpha[u] += to_string(beta_u - h.charge[u]);
    valence_map.emplace(h.alpha[u], beta_u - h.charge[u]); // valence of the new atom is beta_u
    mass_map.emplace(h.alpha[u], mass_map.at(alpha_org));
    // cout << h.alpha[u] << " " << beta_u - h.charge[u] << endl;
  }
}

/**
 * Take a chemical graph H and convert it to its H-suppressed model
 * @param h: a graph possibly containing H-atoms
 * @param Graph g: the H-suppressed model of h
 * @return the number of H atoms in h
 * H_suppressed_convert_old is used when there is no H in the input
 * H_suppressed_convert is used when there is H in the input
 */

//int H_suppressed_convert_old(const Graph& h, Graph& g) {
//  int numH = 0;
//  // get H-vertices from the input graph
//  // unordered_set <Vertex> Hvertex, nonHvertex;
//  unordered_map<Vertex, Vertex> g_h;
//  int count = 0;  // count the number on non-H atoms
//  // itereate over the atoms in graph h
//  for (size_t i = 0; i < h.alpha.size(); ++i) {
//    if (h.alpha[i] != "H") {
//      g.alpha.push_back(h.alpha[i]);
//      g.status.push_back(true);  // by Zhu
//      g_h[count] = i;
//      count++;
//    } else {
//      // increase the count of H-atoms
//      ++numH;
//    }
//  }
//
//  g.numAtom = g.alpha.size();
//  g.numBond = 0;
//  g.beta.resize((g.numAtom), vector<size_t>(g.numAtom));
//  for (size_t i = 0; i < g.numAtom; ++i) {
//    for (size_t j = 0; j < g.numAtom; ++j) {
//      g.beta[i][j] = h.beta[g_h[i]][g_h[j]];
//      if (h.beta[g_h[i]][g_h[j]] != 0) {
//        g.numBond++;
//      }
//    }
//  }
//  g.numBond = g.numBond * 0.5;
//  g.CID = h.CID;
//
//  g.adj.resize(g.numAtom);
//  for (size_t i = 0; i < g.beta.size(); ++i) {
//    for (size_t j = 0; j < g.beta.size(); ++j) {
//      if (g.beta[i][j] != 0) {
//        g.adj[i].push_back(Vertex(j));
//      }
//    }
//  }
//  return numH;
//}

void H_suppressed_convert(const Graph& h, Graph& g) {
	g = h;
	for (Vertex u = 0; u < h.numAtom; ++u) {
    // new, 0430
		if (h.alpha[u] == "H1") {
			for (auto& v : h.adj[u]) {
				g.adj[u].erase(remove(g.adj[u].begin(), g.adj[u].end(), v), g.adj[u].end());
				g.adj[v].erase(remove(g.adj[v].begin(), g.adj[v].end(), u), g.adj[v].end());

			}
			g.status[u] = false;
		}
	}
}


string get_string(input_info& IN_INFO, Graph& g, Vertex u, Vertex parent,
                  const vector<Vertex>& internal_set, size_t height) {
  string ans = "";
  size_t mul;
  if (parent == u) {
    mul = 0;
  } else {
    mul = g.beta[parent][u];
  }
  size_t root_index = IN_INFO.find_index(g.alpha[u]);
  ans += to_string(mul) + to_string(root_index);

  vector<string> subtree_string;
  subtree_string.clear();

  for (auto& v : g.adj[u]) {
    auto itr = find(internal_set.begin(), internal_set.end(), v);
    if (v != parent && itr == internal_set.end()) {
      string v_string = get_string(IN_INFO, g, v, u, internal_set, height + 1);
      subtree_string.push_back(v_string);
    }
  }
  if (height == 1) {  // u is not an internal vertex or a leaf
    size_t subtree_string_size = subtree_string.size();
    for (size_t i = 0; i < IN_INFO.dmax - 1 - subtree_string_size; ++i) {
      string tmp_str = "00";
      subtree_string.push_back(tmp_str);
    }
  }
  if (height == 0 && subtree_string.size() == 0) {
    for (size_t i = 0; i < IN_INFO.dmax; ++i) {
      string tmp_str = "00";
      subtree_string.push_back(tmp_str);
    }
  }
  sort(subtree_string.rbegin(), subtree_string.rend());
  for (auto& tmp : subtree_string) {
    ans += tmp;
  }
  return ans;
}

void get_F1(input_info& IN_INFO, const string& inFileName,
            unordered_set<string>& T_W1_F1, unordered_set<string>& T_W2_F1,
            unordered_set<string>& T_W3_F1) {
  Graph h = read_graph_sdf(inFileName);
  Graph g;
  g.adj.clear();
  g.alpha.clear();
  g.beta.clear();
  g.status.clear();
  H_suppressed_convert(h, g);
  vector<Vertex> internal_set;
  calcInternalVertexSet(g, internal_set);

  for (auto& u : internal_set) {
    size_t neightbour_num_int = 0;
    for (auto& v : g.adj[u]) {
      auto itr = find(internal_set.begin(), internal_set.end(), v);
      if (itr != internal_set.end()) {
        ++neightbour_num_int;
      }
    }
    string fringrtree_string = get_string(IN_INFO, g, u, u, internal_set, 0);
    // cout << fringrtree_string << endl;
    if (neightbour_num_int == 1) {
      T_W1_F1.insert(fringrtree_string);
    } else if (neightbour_num_int == 2) {
      T_W2_F1.insert(fringrtree_string);
    } else {
      T_W3_F1.insert(fringrtree_string);
    }
  }
}

void get_partition(const Graph& g, const vector<Vertex>& coreset,
                   vector<Vertex>& base_vertices,
                   vector_2D<Vertex>& base_edges) {
  base_vertices.clear();
  base_edges.clear();

  for (auto& u : coreset) {
    if (calc_delta(g, coreset, u) > 2) {
      base_vertices.push_back(u);
    }
  }

  if (base_vertices.size() == 0) {
    for (auto& u : coreset) {
      auto itr = find(base_vertices.begin(), base_vertices.end(), u);
      if (itr == base_vertices.end()) {
        base_vertices.push_back(u);
        if (base_vertices.size() == 1) break;
      }
    }
  }

  unordered_set<Vertex> visited;
  vector_2D<Vertex> cycle;
  for (auto& u : base_vertices) {
    visited.insert(u);
    for (auto& v : g.adj[u]) {
      auto itr = find(coreset.begin(), coreset.end(), v);
      if (itr != coreset.end()) {
        vector<Vertex> tmp;
        tmp.push_back(u);
        tmp.push_back(v);
        auto uu = u;
        auto vv = v;
        while (visited.find(vv) == visited.end() &&
               find(base_vertices.begin(), base_vertices.end(), vv) ==
                   base_vertices.end()) {
          visited.insert(vv);
          for (auto& ww : g.adj[vv]) {
            if (ww != uu &&
                find(coreset.begin(), coreset.end(), ww) != coreset.end()) {
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
          if (tmp[0] == tmp[tmp.size() - 1]) {  // avoid self-loop
            cycle.push_back(tmp);
          } else {
            base_edges.push_back(tmp);
          }
        }
      }
    }
  }

  for (auto& tmp : cycle) {
    base_vertices.push_back(tmp[1]);
    vector<Vertex> edge_1 = {tmp[0], tmp[1]};
    vector<Vertex> edge_2;
    for (size_t i = 1; i < tmp.size(); ++i) {
      edge_2.push_back(tmp[i]);
    }
    base_edges.push_back(edge_1);
    base_edges.push_back(edge_2);
  }
}

void get_partition(const Graph& g, const vector<Vertex>& coreset,
                   const vector<Vertex>& internal_set,
                   vector<Vertex>& base_vertices,
                   vector_2D<Vertex>& base_edges) {
  base_vertices.clear();
  base_edges.clear();

  for (auto& u : coreset) {
    if (calc_delta(g, internal_set, u) > 2) {
      base_vertices.push_back(u);
    }
  }

  if (base_vertices.size() == 0) {
    for (auto& u : coreset) {
      auto itr = find(base_vertices.begin(), base_vertices.end(), u);
      if (itr == base_vertices.end()) {
        base_vertices.push_back(u);
        if (base_vertices.size() == 1) break;
      }
    }
  }

  unordered_set<Vertex> visited;
  vector_2D<Vertex> cycle;
  for (auto& u : base_vertices) {
    visited.insert(u);
    for (auto& v : g.adj[u]) {
      auto itr = find(coreset.begin(), coreset.end(), v);
      if (itr != coreset.end()) {
        vector<Vertex> tmp;
        tmp.push_back(u);
        tmp.push_back(v);
        auto uu = u;
        auto vv = v;
        while (visited.find(vv) == visited.end() &&
               find(base_vertices.begin(), base_vertices.end(), vv) ==
                   base_vertices.end()) {
          visited.insert(vv);
          for (auto& ww : g.adj[vv]) {
            if (ww != uu &&
                find(coreset.begin(), coreset.end(), ww) != coreset.end()) {
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
          if (tmp[0] == tmp[tmp.size() - 1]) {  // avoid self-loop
            cycle.push_back(tmp);
          } else {
            base_edges.push_back(tmp);
          }
        }
      }
    }
  }

  for (auto& tmp : cycle) {
    base_vertices.push_back(tmp[1]);
    vector<Vertex> edge_1 = {tmp[0], tmp[1]};
    vector<Vertex> edge_2;
    for (size_t i = 1; i < tmp.size(); ++i) {
      edge_2.push_back(tmp[i]);
    }
    base_edges.push_back(edge_1);
    base_edges.push_back(edge_2);
  }
}

void get_partition(vector<Vertex>& base_vertices, vector_2D<Vertex>& base_edges,
                   vector<size_t>& chLB_v, vector<size_t>& chUB_v,
                   vector<size_t>& chLB_e, vector<size_t>& chUB_e,
                   vector<bool>& fixed_v, vector<bool>& fixed_e,

                   vector_2D <size_t>& fringe_tree_indices_v,
                   vector_2D <size_t>& fringe_tree_indices_e,  

                   vector <Vertex>& core_set,

                   string inFileName) {
  ifstream infile;
  try {
    infile.open(inFileName);
  } catch (const exception& e) {
    cerr << "Couldn't open '" << inFileName << "' for reading!" << endl;
    throw(e);
  }

  if (!infile.is_open()) {
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

  for (size_t i = 0; i < vc_size; ++i) {
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
    if (st2 >> fixed_tmp) {
      if (fixed_tmp != 0) {
        fixed_v.push_back(true);
      } else {
        fixed_v.push_back(false);
      }
    } else {
      fixed_v.push_back(false);
    }

    st2.clear();

    stringstream st3;
    getline(infile, line);
    st3 << line;
    size_t ind_tmp;
    vector <size_t> ft_ind_tmp;
    ft_ind_tmp.clear();

    while (st3 >> ind_tmp){
      ft_ind_tmp.push_back(ind_tmp);
    }
    fringe_tree_indices_v.push_back(ft_ind_tmp);

    st3.clear();
  }

  getline(infile, line);
  st << line;
  size_t ec_size;
  st >> ec_size;
  st.clear();

  for (size_t i = 0; i < ec_size; ++i) {
    stringstream st1;

    getline(infile, line);
    st1 << line;

    Vertex u;
    vector<Vertex> tmp;
    while (st1 >> u) {
      tmp.push_back(u - 1);
      if (find(core_set.begin(), core_set.end(), u - 1) == core_set.end()){
        core_set.push_back(u - 1);
      }
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
    if (st2 >> fixed_tmp) {
      if (fixed_tmp != 0) {
        fixed_e.push_back(true);
      } else {
        fixed_e.push_back(false);
      }
    } else {
      fixed_e.push_back(false);
    }

    st2.clear();

    stringstream st3;
    getline(infile, line);
    st3 << line;
    size_t ind_tmp;
    vector <size_t> ft_ind_tmp;
    ft_ind_tmp.clear();

    while (st3 >> ind_tmp){
      ft_ind_tmp.push_back(ind_tmp);
    }
    fringe_tree_indices_e.push_back(ft_ind_tmp);

    st3.clear();
  }

  infile.close();
}

atom_degree get_atom_degree(const input_info& IN_INFO, const Graph& g,
                            Vertex u) {
  atom_degree a(IN_INFO.find_index(g.alpha[u]), calcDegree(g, u));
  return a;
}

void get_atom_set_from_TS(vector <Atom>& atom_set, const map <size_t, TreeSignature>& TS_map) {

  for (auto& tmp : TS_map){
    for (size_t i = 0; i < tmp.second.delta_2.size(); i += 2){
      bool flag = true;

      // new, 0430
      if (tmp.second.delta_2[i] == "H1") continue;
      string alpha_org = tmp.second.delta_2[i];

      // alpha_org += "_" + to_string(tmp.second.chg[i]);

      for (auto& aa : atom_set){
        if (aa.NAME == tmp.second.delta_2[i]){
          flag = false;
          break;
        }
      }

      if (flag){
        Atom a;
        a.NAME = tmp.second.delta_2[i];

        if (valence_map.find(a.NAME) == valence_map.end()){
          stringstream st_tmp;
          st_tmp << a.NAME[a.NAME.size() - 1];
          int val_tmp;
          st_tmp >> val_tmp;
          string atom_symbol = "";
          size_t j = 0;
          while ((a.NAME[j] < '0') || (a.NAME[j] > '9')){
            atom_symbol += a.NAME[j];
            ++j;
          }
          valence_map.emplace(a.NAME, val_tmp);
          mass_map.emplace(a.NAME, mass_map.at(atom_symbol));
          map_atomic_number.emplace(a.NAME, val_tmp * 200 + map_atomic_number.at(atom_symbol));
        }

          
        a.VALENCE = valence_map.at(a.NAME);
        a.MASS = mass_map.at(a.NAME);
        atom_set.push_back(a);
      }
    }
  }

}

void get_atom_set_tree(input_info& IN_INFO, const Graph& g,
                       const vector<Vertex>& coreset, Vertex u, Vertex parent,
                       bool root_only = false) {
  if (IN_INFO.find_index(g.alpha[u]) == -1) {
    Atom a;
    if (_DEBUG) cout << "u = " << u <<   " g.alpha[u] = " << g.alpha[u] << endl;
    a.NAME = g.alpha[u];
    a.VALENCE = valence_map.at(g.alpha[u]);
    a.MASS = mass_map.at(g.alpha[u]);
    ++IN_INFO.num_kind_atoms;
    IN_INFO.atoms.push_back(a);
  }

  if (root_only) return;

  for (auto& v : g.adj[u]) {
    auto itr = find(coreset.begin(), coreset.end(), v);
    if (!g.status[v] || v == parent || itr != coreset.end()) continue;
    get_atom_set_tree(IN_INFO, g, coreset, v, u);
  }
}

void get_descriptors_tree_dfs(
    input_info& IN_INFO, const Graph& g, const Graph& g_H, const vector<Vertex>& coreset,
    const vector<Vertex>& internal_set,  // set of branch vertices
    // const vector <Vertex>& BranchTree,
    Vertex u, Vertex parent, size_t height, size_t& core_height,
    const map <size_t, TreeSignature>& TS_map
    // the  descriptors will store in IN_INFO.rv
    // need to first get the atom set
) {
  ++IN_INFO.n;
  if (height > core_height) core_height = height;
  if (u != parent) {
    atom_degree a = get_atom_degree(IN_INFO, g, parent);
    atom_degree b = get_atom_degree(IN_INFO, g, u);
    size_t m = g.beta[parent][u];
    auto itr1 = find(internal_set.begin(), internal_set.end(), parent);
    auto itr2 = find(internal_set.begin(), internal_set.end(), u);
    if (itr1 != internal_set.end() && itr2 != internal_set.end()) {
      ++IN_INFO.rv
            .ec_int[IN_INFO.atom_deg_map.at(a)][IN_INFO.atom_deg_map.at(b)][m];
      ++IN_INFO.rv
            .ec_int[IN_INFO.atom_deg_map.at(b)][IN_INFO.atom_deg_map.at(a)][m];
    } 
    // else {
    //   ++IN_INFO.rv
    //         .ec_ex[IN_INFO.atom_deg_map.at(a)][IN_INFO.atom_deg_map.at(b)][m];
    // }

  }

  if (find(internal_set.begin(), internal_set.end(), u) != internal_set.end()){
    // find the fringe tree at the vertex u
    RootedTree RT;
    RT.buildFromCoreVertex(u, internal_set, g_H);
    // get the sequence of fringe tree
    auto RT_tmp = RemoveNullVertex(RT);
    RootedGraph RG(RT_tmp);
    RootedMultiTree RMT(RG, RT_tmp);
    for (Vertex u = 0; u < RMT.n; ++u){
    //RMT.label[u] = atom_map.at(RT_tmp.graph.alpha[u]);
      int len = RT_tmp.graph.alpha[u].size();
      int ind = 0;
      for (int i = 0; i < len; ++i) {
        if (RT_tmp.graph.alpha[u][i] >= '0' && RT_tmp.graph.alpha[u][i] <= '9') {
          ind = i;
          break;
        }
      }
      string str_tmp = RT_tmp.graph.alpha[u].substr(0, ind);
      int val_tmp = stoi(RT_tmp.graph.alpha[u].substr(ind));

      RMT.label[u] = val_tmp * 200 + map_atomic_number.at(str_tmp); 
        // we are assigning colors as val(a)*200+atomic_number
    }
    TreeSignature TS = RMT.getSignature(); 

    // check TS_map, to get the id
    size_t id = 0;
    for (auto& _TS : TS_map){
      if (TS.is_equal(_TS.second)){
        id = _TS.first;
        break;
      }
    }

    // update num_fc and fc_map
    // if (IN_INFO.fc_map.find(id) == IN_INFO.fc_map.end()){
    //   ++IN_INFO.num_fc;
    //   ++IN_INFO.rv1D_size;
    //   IN_INFO.fc_map.emplace(id, IN_INFO.num_fc - 1);
    //   IN_INFO.rv.fc.push_back(1);
    // } else {
    //   ++IN_INFO.rv.fc[IN_INFO.fc_map.at(id)];
    // }
    if (id != 0) ++IN_INFO.rv.fc[IN_INFO.fc_map.at(id)];
  }
  

  size_t deg = calcDegree(g, u);
  // if (deg > 3) IN_INFO.dmax = 4;

  for (auto& v : g.adj[u]) {
    auto itr = find(coreset.begin(), coreset.end(), v);
    if (!g.status[v] || v == parent || itr != coreset.end()) continue;
    get_descriptors_tree_dfs(IN_INFO, g, g_H, coreset, internal_set, v, u,
                             height + 1, core_height, TS_map);
  }
}
// To find resource vector and other related information of a v-component
void get_descriptors_tree(
    input_info& IN_INFO, const Graph& g, const Graph& g_H, const vector<Vertex>& coreset,
    const vector<Vertex>& internal_set,  // set of branch vertices
    Vertex u,
    const vector <Atom>& atom_set,
    const map <size_t, TreeSignature>& TS_map
    // the  descriptors will store in IN_INFO.rv
    // need to first get the atom set
) {
  IN_INFO.dmax = 4;
  get_atom_set_tree(IN_INFO, g, coreset, u, u);
  IN_INFO.atoms = atom_set;
  IN_INFO.num_kind_atoms = atom_set.size();
  sort(IN_INFO.atoms.begin(), IN_INFO.atoms.end());
  IN_INFO.initialize_atom_deg();
  IN_INFO.rv.initialize_ec(IN_INFO.num_kind_atom_degree);
  IN_INFO.initialize_rv_edge_map();

  IN_INFO.rv.initialize_fc(IN_INFO.num_fc);
  IN_INFO.rv1D_size += IN_INFO.num_fc;

  size_t height = 0;

  get_descriptors_tree_dfs(IN_INFO, g, g_H, coreset, internal_set, u, u, 0, height, TS_map);
  if (height > IN_INFO.k_star) {
    IN_INFO.rv.bl += 1;
  }

  if (height > IN_INFO.height) {
    IN_INFO.height = height;
  }

  IN_INFO.rv1D = rv2rv1D(IN_INFO.rv, IN_INFO);

  IN_INFO.delta_1 = calc_delta(g, coreset, u);
  IN_INFO.base_vertices = vector<size_t>({u});
  IN_INFO.rs1_index = u;

  size_t d = calcDegree(g, u);
  size_t a = IN_INFO.find_index(g.alpha[u]);
  int m = IN_INFO.atoms[a].VALENCE; // valence of root
  for (auto& v : g.adj[u]) {
    auto itr = find(coreset.begin(), coreset.end(), v);
    if (itr != coreset.end()) {
      m -= g.beta[u][v];
    }
  }
  // m -= g.charge[u]; // 0509, add charge info
  root_status rs(a, m, d);
  IN_INFO.target_rs1 = rs;
}

void get_atom_set_edge(input_info& IN_INFO, const Graph& g,
                       const vector<Vertex>& coreset,
                       vector<Vertex>& edge_vertex) {
  for (size_t i = 0; i < edge_vertex.size(); ++i) {
    auto u = edge_vertex[i];

    if (i != 0 && i != edge_vertex.size() - 1) {
      get_atom_set_tree(IN_INFO, g, coreset, u, u);
    } else {
      get_atom_set_tree(IN_INFO, g, coreset, u, u, true);
    }
  }
}

void get_descriptors_edge(
    input_info& IN_INFO, const Graph& g, const Graph& g_H, const vector<Vertex>& coreset,
    const vector<Vertex>& internal_set,  // set of branch vertices
    vector<Vertex>& edge_vertex,
    const vector <Atom>& atom_set,
    const map <size_t, TreeSignature>& TS_map
    // the  descriptors will store in IN_INFO.rv
    // need to first get the atom set
) {
  IN_INFO.dmax = 4;
  get_atom_set_edge(IN_INFO, g, coreset, edge_vertex);
  IN_INFO.atoms = atom_set;
  IN_INFO.num_kind_atoms = atom_set.size();
  sort(IN_INFO.atoms.begin(), IN_INFO.atoms.end());
  IN_INFO.initialize_atom_deg();
  IN_INFO.rv.initialize_ec(IN_INFO.num_kind_atom_degree);
  IN_INFO.initialize_rv_edge_map();

  IN_INFO.rv.initialize_fc(IN_INFO.num_fc);
  IN_INFO.rv1D_size += IN_INFO.num_fc;

  for (size_t i = 0; i < edge_vertex.size(); ++i) {
    auto u = edge_vertex[i];
    atom_degree a_u = get_atom_degree(IN_INFO, g, u);
    size_t num_neighbour = 0;

    if (i != 0) {
      auto v = edge_vertex[i - 1];
      atom_degree a_v = get_atom_degree(IN_INFO, g, v);
      size_t m = g.beta[u][v];
      ++IN_INFO.rv.ec_int[IN_INFO.atom_deg_map.at(a_u)]
                         [IN_INFO.atom_deg_map.at(a_v)][m];
      ++IN_INFO.rv.ec_int[IN_INFO.atom_deg_map.at(a_v)]
                         [IN_INFO.atom_deg_map.at(a_u)][m];
    }

    if (i != 0 && i != edge_vertex.size() - 1) {
      size_t height = 0;

      get_descriptors_tree_dfs(IN_INFO, g, g_H, coreset, internal_set, u, u, 0,
                               height, TS_map);
      if (height > 2) {
        IN_INFO.rv.bl += 1;
      }
      if (height > IN_INFO.height) {
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

void output_partition(const Graph& g, const string& outputfilename,
                      const vector<Vertex>& base_vertices,
                      const vector_2D<Vertex>& base_edges,
                      const vector<component>& set_tree_component,
                      const vector<component>& set_edge_component) {
  ofstream output(outputfilename, ios::out);

  output << base_vertices.size() << endl;
  for (size_t i = 0; i < base_vertices.size(); ++i) {
    string sample_str = g.alpha[base_vertices[i]];
    output << base_vertices[i] + 1 << " # " << sample_str << "\n";
    output << "0 " << set_tree_component[i].IN_INFO.height << " 0\n";
  }

  output << base_edges.size() << endl;
  for (size_t i = 0; i < base_edges.size(); ++i) {
    string sample_str = g.alpha[base_edges[i][0]];
    for (size_t j = 1; j < base_edges[i].size(); ++j) {
      sample_str += to_string(g.beta[base_edges[i][j - 1]][base_edges[i][j]]) +
                    g.alpha[base_edges[i][j]];
    }
    for (size_t j = 0; j < base_edges[i].size(); ++j) {
      output << base_edges[i][j] + 1 << " ";
    }
    output << "# " << sample_str << "\n";
    output << "0 " << set_edge_component[i].IN_INFO.height << " 0\n";
  }

  output.close();
}

// add H in a H-suppressed graph
void add_H(size_t & graph_ind,
	vector_2D<size_t>& graph_adj,
	vector<string>& graph_col,
  vector <int>& graph_chg) {
	size_t graph_ind_p = graph_ind;
	for (size_t i = 0; i < graph_ind_p; ++i) {
		int val_i = -graph_chg[i]; // valence of a in G*
		for (size_t j = 0; j < graph_ind_p; ++j)
			if (i != j) val_i += graph_adj[i][j];
    // cout << "i = " << i << " graph_ind = " << graph_ind << " val_i = " << val_i << " graph_col[i] = " << graph_col[i] << endl;
		if (valence_map.at(graph_col[i]) != val_i) { // not all used
			for (int H_counter = 0; H_counter <
				valence_map.at(graph_col[i]) - val_i; ++H_counter) {
				
				graph_col[graph_ind] = "H1";
        graph_chg[graph_ind] = 0;
				graph_adj[graph_ind][i] = 1;
				graph_adj[i][graph_ind] = 1;
				++graph_ind;
			}
		}
	}
}

void output_SDF(size_t& n, size_t& m, vector_2D<size_t>& graph_adj,
                vector<string>& graph_col, vector <int>& graph_chg, string& outputfilename) {
  ofstream output(outputfilename, ios::app);

    size_t real_n = 0;
    size_t real_m = 0;

    for (size_t i = 0; i < n; ++i){
        if (graph_col[i] == ""){
            real_n = i;
            break;
        }
    }
    for (size_t i = 0; i < real_n; ++i){
        for (size_t j = i + 1; j < real_n; ++j){
            if (graph_adj[i][j] > 0){
                ++real_m;
            }
        }
    }

  output << total_num << "\n";
  output << "BH-2LM-M"
         << "\n";
  output << "BH-2LM-M"
         << "\n";
  output << std::setw(3) << real_n << std::setw(3) << real_m
         << "  0  0  0  0  0  0  0  0999 V2000 "
         << "\n";

  for (size_t i = 0; i < real_n; ++i) {
    // new, 0430
    string atom_symbol = "";
    size_t j = 0;
    while ((graph_col[i][j] < '0') || (graph_col[i][j] > '9')){
      atom_symbol += graph_col[i][j];
      ++j;
    }
    int chg_ind = 0;
    if (graph_chg[i] != 0) chg_ind = 4 - graph_chg[i];
    output << "    0.0000    0.0000    0.0000" << std::setw(3) << atom_symbol
           << "  0" << std::setw(3) << chg_ind << "  0  0  0  0  0  0  0  0  0  0"
           << "\n";
  }

  for (size_t i = 0; i < real_n; ++i) {
    for (size_t j = i + 1; j < real_n; ++j) {
      if (graph_adj[i][j] > 0) {
        output << std::setw(3) << i + 1 << std::setw(3) << j + 1 << std::setw(3)
               << graph_adj[i][j] << "  0  0  0  0"
               << "\n";
      }
    }
  }

  size_t n_chg = 0;
  vector <int> vector_chg;
  vector_chg.clear();

  for (size_t i = 0; i < real_n; ++i){

      if (graph_chg[i] != 0){
        ++n_chg;
        vector_chg.push_back(i + 1);
        vector_chg.push_back(graph_chg[i]);
      }
    
  }

  // new, 0430
  if (n_chg == 0){
    output << "M  END"
           << "\n";
  } else {
    output << "M  CHG" << std::setw(3) << n_chg;
    for (auto& tmp : vector_chg){
      output << std::setw(4) << tmp;
    }
    output << "\n";
    output << "M  END"
           << "\n";
  }

  output << "$$$$"
         << "\n";
  output.close();
}

// A function used to generate SDF, especially adding a fringe tree to a vertex
void add_fringe_tree(const input_info& IN_INFO, const FringeTree& T,
                     size_t& graph_ind, vector_2D<size_t>& graph_adj,
                     vector<string>& graph_col, vector <int>& graph_chg) {
  auto& prt = _prt[IN_INFO.dmax];

  vector<size_t> ind_map(T.seq.size());
  size_t root_ind = graph_ind;
  ind_map[0] = root_ind;
  graph_col[root_ind] = IN_INFO.atoms[T.seq[0].second].NAME;
  graph_chg[graph_ind] = T.charge[0]; //add charge

  for (size_t i = 1; i < T.seq.size(); ++i) {
    // if (_DEBUG) cout << "i = " << i << " T.seq.size() = " << T.seq.size() <<
    // endl;
    if (T.seq[i].first != 0) {
      ++graph_ind;
      graph_col[graph_ind] = IN_INFO.atoms[T.seq[i].second].NAME;
      
      graph_chg[graph_ind] = T.charge[i]; //add charge
      // if (_DEBUG) cout <<  "i = " << i << " graph_ind = " << graph_ind << "
      // ind_map[prt[i]] = " << ind_map[prt[i]] << " graph_col[graph_ind] = " <<
      // graph_col[graph_ind] << endl;
      ind_map[i] = graph_ind;
      graph_adj[graph_ind][ind_map[prt[i]]] = T.seq[i].first;
      graph_adj[ind_map[prt[i]]][graph_ind] = T.seq[i].first;
    }
  }
  // if (_DEBUG) cout << "add_fringe_tree end !!!" << endl;
}

void prepare_SDF_branch_tree(const input_info& IN_INFO,
                             const component_result& comp_graph,
                             // vector <int>& st,
                             const map_rv& seq,

                             size_t& graph_ind, vector_2D<size_t>& graph_adj,
                             vector<string>& graph_col,
                             vector <int>& graph_chg,

                             map<Vertex, size_t>& base_vertex_index,

                             size_t root_type,  // 1: W3, 2: W2_core
                             bool tree_component) {
  // cout << "graph_ind = " << graph_ind << endl;
  size_t last_n = 0;
  for (size_t h = 0; h < seq.w_ind.size(); ++h) {
    // if (_DEBUG) cout << "bt h = " <<  h  << " seq.w_ind.size() = " <<
    // seq.w_ind.size() <<  " graph_ind = " << graph_ind << " last_n = " <<
    // last_n << endl; cout << "seq.mul.size() = " << seq.mul.size() << endl;
    if (h != 0) {
      graph_adj[graph_ind][last_n] = seq.mul[h - 1];
      graph_adj[last_n][graph_ind] = seq.mul[h - 1];
    }

    last_n = graph_ind;

    if (h == 0) {
      if (tree_component) {
        base_vertex_index[IN_INFO.rs1_index] = graph_ind;
      }
      if (root_type == 1) {  // always first tree-component then edge-component
        // auto& FT_ind = comp_graph.map_FT_W3.at(make_pair(seq.rs[h],
        // seq.w_ind[h]))[st[h]];
        auto& FT_ind = seq.w_ind[h];
        // cout  << " seq.rs[h] = " << seq.rs[h]  << " seq.w_ind[h] = " <<
        // seq.w_ind[h] << " FT_ind =  " << FT_ind << endl;
        add_fringe_tree(IN_INFO, comp_graph.All_FT_W3[seq.rs[h].color][FT_ind],
                        graph_ind, graph_adj, graph_col, graph_chg);
      } else {
        // auto& FT_ind = comp_graph.map_FT_W2_core.at(make_pair(seq.rs[h],
        // seq.w_ind[h]))[st[h]];
        auto& FT_ind = seq.w_ind[h];
        // cout  << " seq.rs[h] = " << seq.rs[h]  << " seq.w_ind[h] = " <<
        // seq.w_ind[h] << " FT_ind =  " << FT_ind << endl;
        add_fringe_tree(IN_INFO, comp_graph.All_FT_W2_core[seq.rs[h].color][FT_ind],
                        graph_ind, graph_adj, graph_col, graph_chg);
      }
    } else if (h == seq.w_ind.size() - 1) {
      // auto& FT_ind = comp_graph.map_FT_W1.at(make_pair(seq.rs[h],
      // seq.w_ind[h]))[st[h]];
      auto& FT_ind = seq.w_ind[h];
      // cout  << " seq.rs[h] = " << seq.rs[h]  << " seq.w_ind[h] = " <<
      // seq.w_ind[h] << " FT_ind =  " << FT_ind << endl;
      add_fringe_tree(IN_INFO, comp_graph.All_FT_W1[seq.rs[h].color][FT_ind],
                      graph_ind, graph_adj, graph_col, graph_chg);
    } else {
      // auto& FT_ind = comp_graph.map_FT_W2.at(make_pair(seq.rs[h],
      // seq.w_ind[h]))[st[h]];
      auto& FT_ind = seq.w_ind[h];
      // cout  << " seq.rs[h] = " << seq.rs[h]  << " seq.w_ind[h] = " <<
      // seq.w_ind[h] << " FT_ind =  " << FT_ind << endl;
      add_fringe_tree(IN_INFO, comp_graph.All_FT_W2[seq.rs[h].color][FT_ind],
                      graph_ind, graph_adj, graph_col, graph_chg);
    }

    ++graph_ind;
    // if (_DEBUG) cout << "bt h = " <<  h  << " seq.w_ind.size() = " <<
    // seq.w_ind.size() <<  " graph_ind = " << graph_ind << " last_n = " <<
    // last_n << " end!!!" << endl;
  }
}

void prepare_SDF_component(const input_info& IN_INFO,
                           const component_result& comp_graph,
                           const map_tree& seq,

                           size_t& graph_ind, vector_2D<size_t>& graph_adj,
                           vector<string>& graph_col,
                           vector <int>& graph_chg,

                           map<Vertex, size_t>& base_vertex_index,

                           bool tree_component) { 
  if (tree_component) {
    auto& seq_rv = seq.trees[0];
    // vector <int> st(seq_rv.w_ind.size() + 1, 0);
    size_t root_type;
    if (IN_INFO.height > IN_INFO.k_star) {
      root_type = 1;
    } else {
      root_type = 2;
    }
    // cout << "seq_rv.w_ind.size() = " << seq_rv.w_ind.size() << endl;
    prepare_SDF_branch_tree(IN_INFO, comp_graph, seq_rv, graph_ind, graph_adj,
                            graph_col, graph_chg, base_vertex_index, root_type,
                            tree_component);
  } else {
    size_t last_n = graph_ind;
    for (size_t h = 0; h < seq.trees.size(); ++h) {
      // if (_DEBUG) cout << "co h = " << h << " seq.trees.size() = " <<
      // seq.trees.size() << " graph_ind = " << graph_ind << " last_n = " <<
      // last_n << endl; cout << "seq.mul.size() = " << seq.mul.size() << endl;
      // cout << "graph_col.size() = " << graph_col.size() << endl;
      if (h != 0) {
        size_t ind_tmp_1, ind_tmp_2;
        if (h == 1 && h == seq.trees.size() - 1) {
          ind_tmp_1 = base_vertex_index.at(IN_INFO.rs1_index);
          ind_tmp_2 = base_vertex_index.at(IN_INFO.rs2_index);
        } else if (h == 1) {
          ind_tmp_1 = base_vertex_index.at(IN_INFO.rs1_index);
          ind_tmp_2 = graph_ind;
        } else if (h == seq.trees.size() - 1) {
          ind_tmp_1 = last_n;
          ind_tmp_2 = base_vertex_index.at(IN_INFO.rs2_index);
        } else {
          ind_tmp_1 = last_n;
          ind_tmp_2 = graph_ind;
        }
        // if (_DEBUG) cout << "ind_tmp_1 = " << ind_tmp_1 << " ind_tmp_2 = " <<
        // ind_tmp_2 << endl;
        graph_adj[ind_tmp_1][ind_tmp_2] = seq.mul[h - 1];
        graph_adj[ind_tmp_2][ind_tmp_1] = seq.mul[h - 1];
      }

      last_n = graph_ind;

      if (h != 0 && h != seq.trees.size() - 1) {
        auto& seq_rv = seq.trees[h];
        // vector <int> st(seq_rv.w_ind.size() + 1, 0);
        size_t root_type;
        // cout << "seq_rv.w_ind.size() = " << seq_rv.w_ind.size() << endl;
        if (seq_rv.w_ind.size() > 1) {
          root_type = 1;
        } else {
          root_type = 2;
        }

        prepare_SDF_branch_tree(IN_INFO, comp_graph, seq_rv, graph_ind,
                                graph_adj, graph_col, graph_chg, base_vertex_index,
                                root_type, tree_component);
      }
      // if (_DEBUG) cout << "co h = " << h << " seq.trees.size() = " <<
      // seq.trees.size() << " graph_ind = " << graph_ind << " last_n = " <<
      // last_n << " end!!!!!!" << endl;
    }
  }
}

void prepare_SDF_branch_tree_from_g(const Graph& g, vector<bool>& visited,
                                    Vertex u,

                                    size_t& graph_ind,
                                    vector_2D<size_t>& graph_adj,
                                    vector<string>& graph_col,
                                    vector <int>& graph_chg,

                                    map<Vertex, size_t>& base_vertex_index,

                                    bool tree_component) {
  stack<Vertex> _stack;
  _stack.push(u);

  map<Vertex, size_t> g_prt_map_ind;
  map<Vertex, Vertex> g_prt_map;

  base_vertex_index[u] = graph_ind;

  while (!_stack.empty()) {
    Vertex _u = _stack.top();
    _stack.pop();

    graph_col[graph_ind] = g.alpha[_u];

    // add charge
    graph_chg[graph_ind] = g.charge[_u];

    if (_u != u) {
      Vertex prt_ind = g_prt_map_ind.at(_u);
      Vertex _v = g_prt_map.at(_u);
      graph_adj[graph_ind][prt_ind] = g.beta[_u][_v];
      graph_adj[prt_ind][graph_ind] = g.beta[_v][_u];
    }

    for (auto& _v : g.adj[_u]) {
      if (visited[_v]) {
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

void prepare_SDF_component_from_g(const Graph& g, const vector<Vertex>& coreset,
                                  const vector<Vertex>& base_vertices,

                                  size_t& graph_ind,
                                  vector_2D<size_t>& graph_adj,
                                  vector<string>& graph_col,
                                  vector <int>& graph_chg,

                                  map<Vertex, size_t>& base_vertex_index,

                                  bool tree_component) {
  vector<bool> visited(g.alpha.size(), false);
  for (auto& v : coreset) {
    visited[v] = true;
  }

  if (tree_component) {
    Vertex u = base_vertices[0];
    prepare_SDF_branch_tree_from_g(g, visited, u, graph_ind, graph_adj,
                                   graph_col, graph_chg, base_vertex_index,
                                   tree_component);
  } else {
    size_t size_tmp = base_vertices.size();
    size_t last_n = graph_ind;

    for (size_t h = 0; h < size_tmp; ++h) {
      auto& u = base_vertices[h];

      if (h != 0) {
        size_t ind_tmp_1, ind_tmp_2;
        if (h == 1 && h == size_tmp - 1) {
          ind_tmp_1 = base_vertex_index.at(base_vertices[0]);
          ind_tmp_2 = base_vertex_index.at(base_vertices[size_tmp - 1]);
        } else if (h == 1) {
          ind_tmp_1 = base_vertex_index.at(base_vertices[0]);
          ind_tmp_2 = graph_ind;
        } else if (h == size_tmp - 1) {
          ind_tmp_1 = last_n;
          ind_tmp_2 = base_vertex_index.at(base_vertices[size_tmp - 1]);
        } else {
          ind_tmp_1 = last_n;
          ind_tmp_2 = graph_ind;
        }
        // if (_DEBUG) cout << "ind_tmp_1 = " << ind_tmp_1 << " ind_tmp_2 = " <<
        // ind_tmp_2 << endl;
        graph_adj[ind_tmp_1][ind_tmp_2] =
            g.beta[base_vertices[h - 1]][base_vertices[h]];
        graph_adj[ind_tmp_2][ind_tmp_1] =
            g.beta[base_vertices[h - 1]][base_vertices[h]];
      }

      last_n = graph_ind;
      if (h != 0 && h != size_tmp - 1) {
        prepare_SDF_branch_tree_from_g(g, visited, u, graph_ind, graph_adj,
                                       graph_col, graph_chg, base_vertex_index,
                                       tree_component);
      }
    }
  }
}

void prepare_SDF(const Graph& g, const vector<Vertex>& coreset,

                 const vector<component>& set_tree_component,
                 const vector<component>& set_edge_component, size_t& g_n,
                 size_t& g_m, const size_t& core_height_ind,

                 vector<int>& st, string& outputfilename) {
  vector_2D<size_t> graph_adj(g_n, vector<size_t>(g_n, 0));
  vector<string> graph_col(g_n, "");

  vector <int> graph_chg(g_n, 0); // charge info

  size_t graph_ind = 0;

  map<Vertex, size_t> base_vertex_index;
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

  for (size_t i = 0; i < set_tree_component.size(); ++i) {
    auto& IN_INFO = set_tree_component[i].IN_INFO;
    auto& comp_graph = set_tree_component[i].comp_graph;

    if (set_tree_component[i].no_feasible_pair) {
      prepare_SDF_component_from_g(g, coreset, IN_INFO.base_vertices, graph_ind,
                                   graph_adj, graph_col, graph_chg, base_vertex_index,
                                   true);
      continue;
    }

    if (i == core_height_ind) {
      prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_h[st[i]],
                            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index,
                            true);
    } else if (i < core_height_ind) {
      if (st[i] < comp_graph.all_seq_h.size()) {
        prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_h[st[i]],
                              graph_ind, graph_adj, graph_col, graph_chg,
                              base_vertex_index, true);
      } else {
        prepare_SDF_component(
            IN_INFO, comp_graph,
            comp_graph.all_seq[st[i] - comp_graph.all_seq_h.size()], graph_ind,
            graph_adj, graph_col, graph_chg, base_vertex_index, true);
      }
    } else {
      prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq[st[i]],
                            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index,
                            true);
    }
  }

  size_t tree_comp_num = set_tree_component.size();
  for (size_t i = 0; i < set_edge_component.size(); ++i) {
    auto& IN_INFO = set_edge_component[i].IN_INFO;
    auto& comp_graph = set_edge_component[i].comp_graph;

    if (set_edge_component[i].no_feasible_pair) {
      prepare_SDF_component_from_g(g, coreset, IN_INFO.base_vertices, graph_ind,
                                   graph_adj, graph_col, graph_chg, base_vertex_index,
                                   false);
      continue;
    }

    if (i + tree_comp_num == core_height_ind) {
      prepare_SDF_component(
          IN_INFO, comp_graph, comp_graph.all_seq_h[st[i + tree_comp_num]],
          graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
    } else if (i + tree_comp_num < core_height_ind) {
      if (st[i + tree_comp_num] < comp_graph.all_seq_h.size()) {
        prepare_SDF_component(
            IN_INFO, comp_graph, comp_graph.all_seq_h[st[i + tree_comp_num]],
            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
      } else {
        prepare_SDF_component(
            IN_INFO, comp_graph,
            comp_graph
                .all_seq[st[i + tree_comp_num] - comp_graph.all_seq_h.size()],
            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
      }
    } else {
      prepare_SDF_component(
          IN_INFO, comp_graph, comp_graph.all_seq[st[i + tree_comp_num]],
          graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
    }
  }
  ++total_num;
  if (total_num >= num_limit) _stop_gen = true;
  add_H(graph_ind, graph_adj, graph_col, graph_chg);
  output_SDF(g_n, g_m, graph_adj, graph_col, graph_chg, outputfilename);
  // if (_DEBUG) cout << "total_num = " << total_num << endl;
}

void prepare_SDF_unfixed_ch(const Graph& g, const vector<Vertex>& coreset,

                            const vector<component>& set_tree_component,
                            const vector<component>& set_edge_component,
                            size_t& g_n, size_t& g_m,
                            // const size_t& core_height_ind,

                            vector<int>& st, string& outputfilename) {
  vector_2D<size_t> graph_adj(g_n, vector<size_t>(g_n, 0));
  vector<string> graph_col(g_n, "");
  size_t graph_ind = 0;

  vector <int> graph_chg(g_n, 0); // charge info

  map<Vertex, size_t> base_vertex_index;
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

  for (size_t i = 0; i < set_tree_component.size(); ++i) {
    auto& IN_INFO = set_tree_component[i].IN_INFO;
    auto& comp_graph = set_tree_component[i].comp_graph;

    if (set_tree_component[i].no_feasible_pair) {
      prepare_SDF_component_from_g(g, coreset, IN_INFO.base_vertices, graph_ind,
                                   graph_adj, graph_col, graph_chg, base_vertex_index,
                                   true);
      continue;
    }

    // if (i == core_height_ind){
    //     prepare_SDF_component(IN_INFO, comp_graph,
    //     comp_graph.all_seq_h[st[i]], graph_ind, graph_adj, graph_col,
    //     base_vertex_index, true);
    // } else if (i < core_height_ind){
    if (st[i] < comp_graph.all_seq_h.size()) {
      prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_h[st[i]],
                            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index,
                            true);
    } else {
      prepare_SDF_component(
          IN_INFO, comp_graph,
          comp_graph.all_seq[st[i] - comp_graph.all_seq_h.size()], graph_ind,
          graph_adj, graph_col, graph_chg, base_vertex_index, true);
    }
    // } else {
    //     prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq[st[i]],
    //     graph_ind, graph_adj, graph_col, base_vertex_index, true);
    // }
  }

  size_t tree_comp_num = set_tree_component.size();
  for (size_t i = 0; i < set_edge_component.size(); ++i) {
    auto& IN_INFO = set_edge_component[i].IN_INFO;
    auto& comp_graph = set_edge_component[i].comp_graph;

    if (set_edge_component[i].no_feasible_pair) {
      prepare_SDF_component_from_g(g, coreset, IN_INFO.base_vertices, graph_ind,
                                   graph_adj, graph_col, graph_chg, base_vertex_index,
                                   false);
      continue;
    }

    // if (i + tree_comp_num == core_height_ind){
    //     prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_h[st[i
    //     + tree_comp_num]], graph_ind, graph_adj, graph_col,
    //     base_vertex_index, false);
    // } else if (i + tree_comp_num < core_height_ind) {
    if (st[i + tree_comp_num] < comp_graph.all_seq_h.size()) {
      prepare_SDF_component(
          IN_INFO, comp_graph, comp_graph.all_seq_h[st[i + tree_comp_num]],
          graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
    } else {
      prepare_SDF_component(
          IN_INFO, comp_graph,
          comp_graph
              .all_seq[st[i + tree_comp_num] - comp_graph.all_seq_h.size()],
          graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
    }
    // } else {
    //     prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq[st[i +
    //     tree_comp_num]], graph_ind, graph_adj, graph_col, base_vertex_index,
    //     false);
    // }
  }
  ++total_num;
  if (total_num >= num_limit) _stop_gen = true;
  add_H(graph_ind, graph_adj, graph_col, graph_chg);
  output_SDF(g_n, g_m, graph_adj, graph_col, graph_chg, outputfilename);
  // if (_DEBUG) cout << "total_num = " << total_num << endl;
}

void prepare_SDF_loose(const Graph& g, const vector<Vertex>& coreset,

                 const vector<component>& set_tree_component,
                 const vector<component>& set_edge_component, size_t& g_n,
                 size_t& g_m, const size_t& core_height_ind,

                 vector<int>& st, string& outputfilename,
                 const map <size_t, size_t>& fc_map, const vector <unsigned short>& fc_all,
                 const map <size_t, TreeSignature>& TS_map) {
  vector_2D<size_t> graph_adj(g_n, vector<size_t>(g_n, 0));
  vector<string> graph_col(g_n, "");
  size_t graph_ind = 0;

  vector <int> graph_chg(g_n, 0); // charge info

  map<Vertex, size_t> base_vertex_index;
  base_vertex_index.clear();

  int num_loose = 0;

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

  for (size_t i = 0; i < set_tree_component.size(); ++i) {
    auto& IN_INFO = set_tree_component[i].IN_INFO;
    auto& comp_graph = set_tree_component[i].comp_graph;

    if (set_tree_component[i].no_feasible_pair_loose) {
      prepare_SDF_component_from_g(g, coreset, IN_INFO.base_vertices, graph_ind,
                                   graph_adj, graph_col, graph_chg, base_vertex_index,
                                   true);
      continue;
    }

    if (i == core_height_ind) {
      if (st[i] < comp_graph.all_seq_h.size()){
        prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_h[st[i]],
                              graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index,
                              true);
      } else {
        prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_h_loose[st[i] - comp_graph.all_seq_h.size()],
                              graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index,
                              true);
        ++num_loose;
      }
      
    } else if (i < core_height_ind) {
      if (st[i] < comp_graph.all_seq_h.size()) {
        prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_h[st[i]],
                              graph_ind, graph_adj, graph_col, graph_chg,
                              base_vertex_index, true);
      } else if (st[i] < comp_graph.all_seq_h.size() + comp_graph.all_seq.size()) {
        prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq[st[i] - comp_graph.all_seq_h.size()],
                              graph_ind, graph_adj, graph_col, graph_chg,
                              base_vertex_index, true);
      } else if (st[i] < comp_graph.all_seq_h.size() + comp_graph.all_seq.size() + comp_graph.all_seq_h_loose.size()) {
        prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_h_loose[st[i] - comp_graph.all_seq_h.size() - comp_graph.all_seq.size()],
                              graph_ind, graph_adj, graph_col, graph_chg,
                              base_vertex_index, true);
        ++num_loose;
      } else {
        prepare_SDF_component(
            IN_INFO, comp_graph,
            comp_graph.all_seq_loose[st[i] - comp_graph.all_seq_h.size() - comp_graph.all_seq.size() - comp_graph.all_seq_h_loose.size()], graph_ind,
            graph_adj, graph_col, graph_chg, base_vertex_index, true);
        ++num_loose;
      }
    } else {
      if (st[i] < comp_graph.all_seq.size()){
        prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq[st[i]],
                              graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index,
                              true);
      } else {
        prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_loose[st[i] - comp_graph.all_seq.size()],
                              graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index,
                              true);
        ++num_loose;
      }
      
    }
  }

  size_t tree_comp_num = set_tree_component.size();
  for (size_t i = 0; i < set_edge_component.size(); ++i) {
    auto& IN_INFO = set_edge_component[i].IN_INFO;
    auto& comp_graph = set_edge_component[i].comp_graph;

    if (set_edge_component[i].no_feasible_pair_loose) {
      prepare_SDF_component_from_g(g, coreset, IN_INFO.base_vertices, graph_ind,
                                   graph_adj, graph_col, graph_chg, base_vertex_index,
                                   false);
      continue;
    }

    if (i + tree_comp_num == core_height_ind) {
      if (st[i + tree_comp_num] < comp_graph.all_seq_h.size()){
        prepare_SDF_component(
            IN_INFO, comp_graph, comp_graph.all_seq_h[st[i + tree_comp_num]],
            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
      } else {
        prepare_SDF_component(
            IN_INFO, comp_graph, comp_graph.all_seq_h_loose[st[i + tree_comp_num] - comp_graph.all_seq_h.size()],
            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
        ++num_loose;
      }
    } else if (i + tree_comp_num < core_height_ind) {
      if (st[i + tree_comp_num] < comp_graph.all_seq_h.size()) {
        prepare_SDF_component(
            IN_INFO, comp_graph, comp_graph.all_seq_h[st[i + tree_comp_num]],
            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
      } else if (st[i + tree_comp_num] < comp_graph.all_seq_h.size() + comp_graph.all_seq.size()) {
        prepare_SDF_component(
            IN_INFO, comp_graph, comp_graph.all_seq[st[i + tree_comp_num] - comp_graph.all_seq_h.size()],
            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
      } else if (st[i + tree_comp_num] < comp_graph.all_seq_h.size() + comp_graph.all_seq.size() + comp_graph.all_seq_h_loose.size()) {
        prepare_SDF_component(
            IN_INFO, comp_graph, comp_graph.all_seq_h_loose[st[i + tree_comp_num] - comp_graph.all_seq_h.size() - comp_graph.all_seq.size()],
            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
        ++num_loose;
      } else {
        prepare_SDF_component(
            IN_INFO, comp_graph,
            comp_graph
                .all_seq_loose[st[i + tree_comp_num] - comp_graph.all_seq_h.size() - comp_graph.all_seq.size() - comp_graph.all_seq_h_loose.size()],
            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
        ++num_loose;
      }
    } else {
      if (st[i + tree_comp_num] < comp_graph.all_seq.size()){
        prepare_SDF_component(
            IN_INFO, comp_graph, comp_graph.all_seq[st[i + tree_comp_num]],
            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
      } else {
        prepare_SDF_component(
            IN_INFO, comp_graph, comp_graph.all_seq_loose[st[i + tree_comp_num] - comp_graph.all_seq.size()],
            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
        ++num_loose;
      }
      
    }
  }

  if (num_loose < 2) return;

  // TO DO add check fc

  Graph _g;
  size_t real_n = 0;

  for (size_t i = 0; i < graph_col.size(); ++i){
      if (graph_col[i] == ""){
          real_n = i;
          break;
      }
  }
  _g.initialize(real_n);
  for (Vertex u = 0; u < real_n; ++u){
    _g.alpha[u] = graph_col[u];
    _g.status[u] = true;
    _g.charge[u] = graph_chg[u];
  }
  for  (Vertex u = 0; u < real_n; ++u){
    for  (Vertex v = 0; v < real_n; ++v){
      if (graph_adj[u][v] != 0){
        _g.adj[u].push_back(v);
        _g.beta[u][v] = graph_adj[u][v];
      }
    }
  }

  vector <unsigned short> fc(fc_all.size(), 0);

  vector <RootedTree> RTs = getFringeTrees(_g);
  for (auto& RT : RTs){
    auto RT_tmp = RemoveNullVertex(RT);
    RootedGraph RG(RT_tmp);
    RootedMultiTree RMT(RG, RT_tmp);

    for (Vertex u = 0; u < RMT.n; ++u){
      // RMT.label[u] = map_atomic_number.at(RT_tmp.graph.alpha[u]);
      int len = RT_tmp.graph.alpha[u].size();
      int ind = 0;
      for (int i = 0; i < len; ++i) {
        if (RT_tmp.graph.alpha[u][i] >= '0' && RT_tmp.graph.alpha[u][i] <= '9') {
          ind = i;
          break;
        }
      }
      string str_tmp = RT_tmp.graph.alpha[u].substr(0, ind);
      int val_tmp = stoi(RT_tmp.graph.alpha[u].substr(ind));

      RMT.label[u] = val_tmp * 200 + map_atomic_number.at(str_tmp); 
        // we are assigning colors as val(a)*200+atomic_number
    }
    TreeSignature TS = RMT.getSignature();

    for (auto& _TS : TS_map){
      if (TS.is_equal(_TS.second)){
        size_t id = _TS.first;
        ++fc[fc_map.at(id)];
        break;
      }
    }
  }

  if (lexicographical_comparison(fc, fc_all) != 0) return;

  ++total_num;
  if (total_num >= num_limit) _stop_gen = true;
  add_H(graph_ind, graph_adj, graph_col, graph_chg);
  output_SDF(g_n, g_m, graph_adj, graph_col, graph_chg, outputfilename);
  // if (_DEBUG) cout << "total_num = " << total_num << endl;
}

void prepare_SDF_unfixed_ch_loose(const Graph& g, const vector<Vertex>& coreset,

                            const vector<component>& set_tree_component,
                            const vector<component>& set_edge_component,
                            size_t& g_n, size_t& g_m,
                            // const size_t& core_height_ind,

                            vector<int>& st, string& outputfilename,
                            const map <size_t, size_t>& fc_map, const vector <unsigned short>& fc_all,
                            const map <size_t, TreeSignature>& TS_map) {
  vector_2D<size_t> graph_adj(g_n, vector<size_t>(g_n, 0));
  vector<string> graph_col(g_n, "");
  size_t graph_ind = 0;

  vector <int> graph_chg(g_n, 0); // charge info

  map<Vertex, size_t> base_vertex_index;
  base_vertex_index.clear();

  int num_loose = 0;

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

  for (size_t i = 0; i < set_tree_component.size(); ++i) {
    auto& IN_INFO = set_tree_component[i].IN_INFO;
    auto& comp_graph = set_tree_component[i].comp_graph;

    if (set_tree_component[i].no_feasible_pair_loose) {
      prepare_SDF_component_from_g(g, coreset, IN_INFO.base_vertices, graph_ind,
                                   graph_adj, graph_col, graph_chg, base_vertex_index,
                                   true);
      continue;
    }

    // if (i == core_height_ind){
    //     prepare_SDF_component(IN_INFO, comp_graph,
    //     comp_graph.all_seq_h[st[i]], graph_ind, graph_adj, graph_col, graph_chg,
    //     base_vertex_index, true);
    // } else if (i < core_height_ind){
    if (st[i] < comp_graph.all_seq_h.size()) {
        prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_h[st[i]],
                              graph_ind, graph_adj, graph_col, graph_chg,
                              base_vertex_index, true);
      } else if (st[i] < comp_graph.all_seq_h.size() + comp_graph.all_seq.size()) {
        prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq[st[i] - comp_graph.all_seq_h.size()],
                              graph_ind, graph_adj, graph_col, graph_chg,
                              base_vertex_index, true);
      } else if (st[i] < comp_graph.all_seq_h.size() + comp_graph.all_seq.size() + comp_graph.all_seq_h_loose.size()) {
        prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_h_loose[st[i] - comp_graph.all_seq_h.size() - comp_graph.all_seq.size()],
                              graph_ind, graph_adj, graph_col, graph_chg,
                              base_vertex_index, true);
        ++num_loose;
      } else {
        prepare_SDF_component(
            IN_INFO, comp_graph,
            comp_graph.all_seq_loose[st[i] - comp_graph.all_seq_h.size() - comp_graph.all_seq.size() - comp_graph.all_seq_h_loose.size()], graph_ind,
            graph_adj, graph_col, graph_chg, base_vertex_index, true);
        ++num_loose;
      }
    // } else {
    //     prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq[st[i]],
    //     graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, true);
    // }
  }

  size_t tree_comp_num = set_tree_component.size();
  for (size_t i = 0; i < set_edge_component.size(); ++i) {
    auto& IN_INFO = set_edge_component[i].IN_INFO;
    auto& comp_graph = set_edge_component[i].comp_graph;

    if (set_edge_component[i].no_feasible_pair_loose) {
      prepare_SDF_component_from_g(g, coreset, IN_INFO.base_vertices, graph_ind,
                                   graph_adj, graph_col, graph_chg, base_vertex_index,
                                   false);
      continue;
    }

    // if (i + tree_comp_num == core_height_ind){
    //     prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq_h[st[i
    //     + tree_comp_num]], graph_ind, graph_adj, graph_col, graph_chg,
    //     base_vertex_index, false);
    // } else if (i + tree_comp_num < core_height_ind) {
    if (st[i + tree_comp_num] < comp_graph.all_seq_h.size()) {
        prepare_SDF_component(
            IN_INFO, comp_graph, comp_graph.all_seq_h[st[i + tree_comp_num]],
            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
      } else if (st[i + tree_comp_num] < comp_graph.all_seq_h.size() + comp_graph.all_seq.size()) {
        prepare_SDF_component(
            IN_INFO, comp_graph, comp_graph.all_seq[st[i + tree_comp_num] - comp_graph.all_seq_h.size()],
            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
      } else if (st[i + tree_comp_num] < comp_graph.all_seq_h.size() + comp_graph.all_seq.size() + comp_graph.all_seq_h_loose.size()) {
        prepare_SDF_component(
            IN_INFO, comp_graph, comp_graph.all_seq_h_loose[st[i + tree_comp_num] - comp_graph.all_seq_h.size() - comp_graph.all_seq.size()],
            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
        ++num_loose;
      } else {
        prepare_SDF_component(
            IN_INFO, comp_graph,
            comp_graph
                .all_seq_loose[st[i + tree_comp_num] - comp_graph.all_seq_h.size() - comp_graph.all_seq.size() - comp_graph.all_seq_h_loose.size()],
            graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index, false);
        ++num_loose;
      }
    // } else {
    //     prepare_SDF_component(IN_INFO, comp_graph, comp_graph.all_seq[st[i +
    //     tree_comp_num]], graph_ind, graph_adj, graph_col, graph_chg, base_vertex_index,
    //     false);
    // }
  }

  if (num_loose < 2) return;

  // TO DO add check fc
  Graph _g;
  size_t real_n = 0;

  for (size_t i = 0; i < graph_col.size(); ++i){
      if (graph_col[i] == ""){
          real_n = i;
          break;
      }
  }
  _g.initialize(real_n);
  for (Vertex u = 0; u < real_n; ++u){
    _g.alpha[u] = graph_col[u];
    _g.status[u] = true;
    _g.charge[u] = graph_chg[u];
  }
  for  (Vertex u = 0; u < real_n; ++u){
    for  (Vertex v = 0; v < real_n; ++v){
      if (graph_adj[u][v] != 0){
        _g.adj[u].push_back(v);
        _g.beta[u][v] = graph_adj[u][v];
      }
    }
  }

  vector <unsigned short> fc(fc_all.size(), 0);

  vector <RootedTree> RTs = getFringeTrees(_g);
  for (auto& RT : RTs){
    auto RT_tmp = RemoveNullVertex(RT);
    RootedGraph RG(RT_tmp);
    RootedMultiTree RMT(RG, RT_tmp);

    for (Vertex u = 0; u < RMT.n; ++u){
      // RMT.label[u] = map_atomic_number.at(RT_tmp.graph.alpha[u]);
      int len = RT_tmp.graph.alpha[u].size();
      int ind = 0;
      for (int i = 0; i < len; ++i) {
        if (RT_tmp.graph.alpha[u][i] >= '0' && RT_tmp.graph.alpha[u][i] <= '9') {
          ind = i;
          break;
        }
      }
      string str_tmp = RT_tmp.graph.alpha[u].substr(0, ind);
      int val_tmp = stoi(RT_tmp.graph.alpha[u].substr(ind));

      RMT.label[u] = val_tmp * 200 + map_atomic_number.at(str_tmp); 
        // we are assigning colors as val(a)*200+atomic_number
    }
    TreeSignature TS = RMT.getSignature();

    for (auto& _TS : TS_map){
      if (TS.is_equal(_TS.second)){
        size_t id = _TS.first;
        ++fc[fc_map.at(id)];
        break;
      }
    }
  }

  if (lexicographical_comparison(fc, fc_all) != 0) return;

  ++total_num;
  if (total_num >= num_limit) _stop_gen = true;
  add_H(graph_ind, graph_adj, graph_col, graph_chg);
  output_SDF(g_n, g_m, graph_adj, graph_col, graph_chg, outputfilename);
  // if (_DEBUG) cout << "total_num = " << total_num << endl;
}

void generate(const Graph& g, const vector<Vertex>& coreset,

              const vector<component>& set_tree_component,
              const vector<component>& set_edge_component, size_t& g_n,
              size_t& g_m, const size_t& core_height_ind,

              size_t& possible_num, size_t& lower_bound,

              string& outputfilename) {
  size_t n_comp = set_tree_component.size() + set_edge_component.size();
  vector<int> st(n_comp + 1, -1);
  vector<size_t> st_max(n_comp + 1, 0);
  size_t possible_num_tmp = 1;
  size_t lower_bound_tmp = 1;

  for (size_t i = 0; i < n_comp; ++i) {
    // cout << "i = " << i << endl;
    if (i < set_tree_component.size()) {
      // cout << set_tree_component[i].comp_graph.num_h << " " <<
      // set_tree_component[i].comp_graph.num << endl;
      if (i == core_height_ind) {
        if (set_tree_component[i].no_feasible_pair &&
            set_tree_component[i].with_core_height) {
          st_max[i] = 1;
        } else {
          st_max[i] = set_tree_component[i].comp_graph.all_seq_h.size();
        }
        lower_bound_tmp *= set_tree_component[i].comp_graph.num_h;
      } else if (i < core_height_ind) {
        if (set_tree_component[i].no_feasible_pair) {
          st_max[i] = 1;
        } else {
          st_max[i] = set_tree_component[i].comp_graph.all_seq_h.size() +
                      set_tree_component[i].comp_graph.all_seq.size();
        }
        lower_bound_tmp *= (set_tree_component[i].comp_graph.num_h +
                            set_tree_component[i].comp_graph.num);
      } else {
        if (set_tree_component[i].no_feasible_pair &&
            !set_tree_component[i].with_core_height) {
          st_max[i] = 1;
        } else {
          st_max[i] = set_tree_component[i].comp_graph.all_seq.size();
        }
        lower_bound_tmp *= set_tree_component[i].comp_graph.num;
      }
    } else {
      // cout << set_edge_component[i -
      // set_tree_component.size()].comp_graph.num_h << " " <<
      // set_edge_component[i - set_tree_component.size()].comp_graph.num <<
      // endl;
      if (i == core_height_ind) {
        if (set_edge_component[i - set_tree_component.size()]
                .no_feasible_pair &&
            set_edge_component[i - set_tree_component.size()]
                .with_core_height) {
          st_max[i] = 1;
        } else {
          st_max[i] = set_edge_component[i - set_tree_component.size()]
                          .comp_graph.all_seq_h.size();
        }
        lower_bound_tmp *=
            set_edge_component[i - set_tree_component.size()].comp_graph.num_h;
      } else if (i < core_height_ind) {
        if (set_edge_component[i - set_tree_component.size()]
                .no_feasible_pair) {
          st_max[i] = 1;
        } else {
          st_max[i] = set_edge_component[i - set_tree_component.size()]
                          .comp_graph.all_seq_h.size() +
                      set_edge_component[i - set_tree_component.size()]
                          .comp_graph.all_seq.size();
        }
        lower_bound_tmp *=
            (set_edge_component[i - set_tree_component.size()]
                 .comp_graph.num_h +
             set_edge_component[i - set_tree_component.size()].comp_graph.num);
      } else {
        if (set_edge_component[i - set_tree_component.size()]
                .no_feasible_pair &&
            !set_edge_component[i - set_tree_component.size()]
                 .with_core_height) {
          st_max[i] = 1;
        } else {
          st_max[i] = set_edge_component[i - set_tree_component.size()]
                          .comp_graph.all_seq.size();
        }
        lower_bound_tmp *=
            set_edge_component[i - set_tree_component.size()].comp_graph.num;
      }
    }
    if (st_max[i] == 0) return;
    possible_num_tmp *= st_max[i];
  }
  possible_num += possible_num_tmp;
  lower_bound += lower_bound_tmp;
  // cout << "lower_bound_tmp = " << lower_bound_tmp << endl;

  int h = 0;

  while (h >= 0) {
    if (_stop_gen) return;
    // if (_DEBUG) {
    //     cout << "h = "  << h<< endl;
    //     for (size_t i = 0; i < st.size(); ++i){
    //         cout << st[i] << " " ;
    //     }
    //     cout << endl;
    // }

    if (h == set_tree_component.size() + set_edge_component.size()) {
      prepare_SDF(g, coreset, set_tree_component, set_edge_component, g_n, g_m,
                  core_height_ind, st, outputfilename);
      --h;
      // if (_DEBUG) cout << "h = " << h<< " end!!!!!"<< endl;
    } else {
      if (st[h] + 1 < st_max[h]) {
        ++st[h];
        ++h;
        st[h] = -1;
      } else {
        --h;
      }
    }
  }
}

void generate_unfixed_ch(const Graph& g, const vector<Vertex>& coreset,

                         const vector<component>& set_tree_component,
                         const vector<component>& set_edge_component,
                         size_t& g_n, size_t& g_m,
                         // const size_t& core_height_ind,

                         size_t& possible_num, size_t& lower_bound,

                         string& outputfilename) {
  size_t n_comp = set_tree_component.size() + set_edge_component.size();
  vector<int> st(n_comp + 1, -1);
  vector<size_t> st_max(n_comp + 1, 0);
  size_t possible_num_tmp = 1;
  size_t lower_bound_tmp = 1;

  for (size_t i = 0; i < n_comp; ++i) {
    // cout << "i = " << i << endl;
    if (i < set_tree_component.size()) {
      // cout << set_tree_component[i].comp_graph.num_h << " " <<
      // set_tree_component[i].comp_graph.num << endl; if (i ==
      // core_height_ind){
      //     if (set_tree_component[i].no_feasible_pair &&
      //     set_tree_component[i].with_core_height){
      //         st_max[i] = 1;
      //     } else {
      //         st_max[i] = set_tree_component[i].comp_graph.all_seq_h.size();
      //     }
      //     lower_bound_tmp *= set_tree_component[i].comp_graph.num_h;
      // } else if (i < core_height_ind){
      if (set_tree_component[i].no_feasible_pair) {
        st_max[i] = 1;
      } else {
        st_max[i] = set_tree_component[i].comp_graph.all_seq_h.size() +
                    set_tree_component[i].comp_graph.all_seq.size();
      }
      lower_bound_tmp *= (set_tree_component[i].comp_graph.num_h +
                          set_tree_component[i].comp_graph.num);
      // } else {
      //     if (set_tree_component[i].no_feasible_pair &&
      //     !set_tree_component[i].with_core_height){
      //         st_max[i] = 1;
      //     } else {
      //         st_max[i] = set_tree_component[i].comp_graph.all_seq.size();
      //     }
      //     lower_bound_tmp *= set_tree_component[i].comp_graph.num;
      // }
    } else {
      // cout << set_edge_component[i -
      // set_tree_component.size()].comp_graph.num_h << " " <<
      // set_edge_component[i - set_tree_component.size()].comp_graph.num <<
      // endl; if (i == core_height_ind){
      //     if (set_edge_component[i -
      //     set_tree_component.size()].no_feasible_pair && set_edge_component[i
      //     - set_tree_component.size()].with_core_height){
      //         st_max[i] = 1;
      //     } else {
      //         st_max[i] = set_edge_component[i -
      //         set_tree_component.size()].comp_graph.all_seq_h.size();
      //     }
      //     lower_bound_tmp *= set_edge_component[i -
      //     set_tree_component.size()].comp_graph.num_h;
      // } else if (i < core_height_ind){
      if (set_edge_component[i - set_tree_component.size()].no_feasible_pair) {
        st_max[i] = 1;
      } else {
        st_max[i] = set_edge_component[i - set_tree_component.size()]
                        .comp_graph.all_seq_h.size() +
                    set_edge_component[i - set_tree_component.size()]
                        .comp_graph.all_seq.size();
      }
      lower_bound_tmp *=
          (set_edge_component[i - set_tree_component.size()].comp_graph.num_h +
           set_edge_component[i - set_tree_component.size()].comp_graph.num);
      // } else {
      //     if (set_edge_component[i -
      //     set_tree_component.size()].no_feasible_pair &&
      //     !set_edge_component[i -
      //     set_tree_component.size()].with_core_height){
      //         st_max[i] = 1;
      //     } else {
      //         st_max[i] = set_edge_component[i -
      //         set_tree_component.size()].comp_graph.all_seq.size();
      //     }
      //     lower_bound_tmp *= set_edge_component[i -
      //     set_tree_component.size()].comp_graph.num;
      // }
    }
    if (st_max[i] == 0) return;
    possible_num_tmp *= st_max[i];
  }
  possible_num += possible_num_tmp;
  lower_bound += lower_bound_tmp;
  // cout << "lower_bound_tmp = " << lower_bound_tmp << endl;

  int h = 0;

  while (h >= 0) {
    if (_stop_gen) return;
    // if (_DEBUG) {
    //     cout << "h = "  << h<< endl;
    //     for (size_t i = 0; i < st.size(); ++i){
    //         cout << st[i] << " " ;
    //     }
    //     cout << endl;
    // }

    if (h == set_tree_component.size() + set_edge_component.size()) {
      prepare_SDF_unfixed_ch(g, coreset, set_tree_component, set_edge_component,
                             g_n, g_m, st, outputfilename);
      --h;
      // if (_DEBUG) cout << "h = " << h<< " end!!!!!"<< endl;
    } else {
      if (st[h] + 1 < st_max[h]) {
        ++st[h];
        ++h;
        st[h] = -1;
      } else {
        --h;
      }
    }
  }
}

// change 0610
void generate_loose(const Graph& g, const vector<Vertex>& coreset,

              const vector<component>& set_tree_component,
              const vector<component>& set_edge_component, size_t& g_n,
              size_t& g_m, const size_t& core_height_ind,

              size_t& possible_num, size_t& lower_bound,

              string& outputfilename,
              const map <size_t, size_t>& fc_map, const vector <unsigned short>& fc_all,
              const map <size_t, TreeSignature>& TS_map) {
  size_t n_comp = set_tree_component.size() + set_edge_component.size();
  vector<int> st(n_comp + 1, -1);
  vector<size_t> st_max(n_comp + 1, 0);
  size_t possible_num_tmp = 1;
  size_t lower_bound_tmp = 1;
  
  // cout << "core_height_ind = " << core_height_ind << endl;
  for (size_t i = 0; i < n_comp; ++i) {
    // cout << "i = " << i << endl;
    if (i < set_tree_component.size()) {
      // cout << set_tree_component[i].comp_graph.num_h << " " <<
      // set_tree_component[i].comp_graph.num << endl;
      if (i == core_height_ind) {
        if (set_tree_component[i].no_feasible_pair_loose &&
            set_tree_component[i].with_core_height) {
          st_max[i] = 1;
        } else {
          st_max[i] = set_tree_component[i].comp_graph.all_seq_h.size() + set_tree_component[i].comp_graph.all_seq_h_loose.size();
        }
      } else if (i < core_height_ind) {
        if (set_tree_component[i].no_feasible_pair_loose) {
          st_max[i] = 1;
        } else {
          st_max[i] = set_tree_component[i].comp_graph.all_seq_h.size() + set_tree_component[i].comp_graph.all_seq.size() +
                      set_tree_component[i].comp_graph.all_seq_h_loose.size() + set_tree_component[i].comp_graph.all_seq_loose.size();
        }
      } else {
        if (set_tree_component[i].no_feasible_pair_loose &&
            !set_tree_component[i].with_core_height) {
          st_max[i] = 1;
        } else {
          st_max[i] = set_tree_component[i].comp_graph.all_seq.size() + set_tree_component[i].comp_graph.all_seq_loose.size();
        }
      }
    } else {
      // cout << set_edge_component[i -
      // set_tree_component.size()].comp_graph.num_h << " " <<
      // set_edge_component[i - set_tree_component.size()].comp_graph.num <<
      // endl;
      if (i == core_height_ind) {
        if (set_edge_component[i - set_tree_component.size()]
                .no_feasible_pair_loose &&
            set_edge_component[i - set_tree_component.size()]
                .with_core_height) {
          st_max[i] = 1;
        } else {
          st_max[i] = set_edge_component[i - set_tree_component.size()]
                          .comp_graph.all_seq_h.size() + 
                      set_edge_component[i - set_tree_component.size()]
                          .comp_graph.all_seq_h_loose.size();
        }
      } else if (i < core_height_ind) {
        if (set_edge_component[i - set_tree_component.size()]
                .no_feasible_pair_loose) {
          st_max[i] = 1;
        } else {
          st_max[i] = set_edge_component[i - set_tree_component.size()]
                          .comp_graph.all_seq_h.size() +
                      set_edge_component[i - set_tree_component.size()]
                          .comp_graph.all_seq.size() + 
                      set_edge_component[i - set_tree_component.size()]
                          .comp_graph.all_seq_h_loose.size() +
                      set_edge_component[i - set_tree_component.size()]
                          .comp_graph.all_seq_loose.size();
        }
      } else {
        if (set_edge_component[i - set_tree_component.size()]
                .no_feasible_pair_loose &&
            !set_edge_component[i - set_tree_component.size()]
                 .with_core_height) {
          st_max[i] = 1;
        } else {
          st_max[i] = set_edge_component[i - set_tree_component.size()]
                          .comp_graph.all_seq.size() + 
                      set_edge_component[i - set_tree_component.size()]
                          .comp_graph.all_seq_loose.size();
        }
      }
    }
    // for (auto& s : st_max) cout << s << " "; cout << endl; 
    if (st_max[i] == 0) return;
  }
  // cout << endl;
  // cout << "lower_bound_tmp = " << lower_bound_tmp << endl;

  int h = 0;

  while (h >= 0) {
    if (_stop_gen) return;
    // if (_DEBUG) {
    //     cout << "h = "  << h<< endl;
    //     for (size_t i = 0; i < st.size(); ++i){
    //         cout << st[i] << " " ;
    //     }
    //     cout << endl;
    // }

    if (h == set_tree_component.size() + set_edge_component.size()) {
      prepare_SDF_loose(g, coreset, set_tree_component, set_edge_component, g_n, g_m,
                  core_height_ind, st, outputfilename, fc_map, fc_all, TS_map);
      --h;
      // if (_DEBUG) cout << "h = " << h<< " end!!!!!"<< endl;
    } else {
      if (st[h] + 1 < st_max[h]) {
        ++st[h];
        ++h;
        st[h] = -1;
      } else {
        --h;
      }
    }
  }
}

void generate_unfixed_ch_loose(const Graph& g, const vector<Vertex>& coreset,

                         const vector<component>& set_tree_component,
                         const vector<component>& set_edge_component,
                         size_t& g_n, size_t& g_m,
                         // const size_t& core_height_ind,

                         size_t& possible_num, size_t& lower_bound,

                         string& outputfilename,
                         const map <size_t, size_t>& fc_map, const vector <unsigned short>& fc_all,
                         const map <size_t, TreeSignature>& TS_map) {
  size_t n_comp = set_tree_component.size() + set_edge_component.size();
  vector<int> st(n_comp + 1, -1);
  vector<size_t> st_max(n_comp + 1, 0);
  size_t possible_num_tmp = 1;
  size_t lower_bound_tmp = 1;

  for (size_t i = 0; i < n_comp; ++i) {
    // cout << "i = " << i << endl;
    if (i < set_tree_component.size()) {
      // cout << set_tree_component[i].comp_graph.num_h << " " <<
      // set_tree_component[i].comp_graph.num << endl; if (i ==
      // core_height_ind){
      //     if (set_tree_component[i].no_feasible_pair &&
      //     set_tree_component[i].with_core_height){
      //         st_max[i] = 1;
      //     } else {
      //         st_max[i] = set_tree_component[i].comp_graph.all_seq_h.size();
      //     }
      //     lower_bound_tmp *= set_tree_component[i].comp_graph.num_h;
      // } else if (i < core_height_ind){
      if (set_tree_component[i].no_feasible_pair_loose) {
        st_max[i] = 1;
      } else {
        st_max[i] = set_tree_component[i].comp_graph.all_seq_h.size() + set_tree_component[i].comp_graph.all_seq.size() +
                    set_tree_component[i].comp_graph.all_seq_h_loose.size() + set_tree_component[i].comp_graph.all_seq_loose.size();
      }
      // } else {
      //     if (set_tree_component[i].no_feasible_pair &&
      //     !set_tree_component[i].with_core_height){
      //         st_max[i] = 1;
      //     } else {
      //         st_max[i] = set_tree_component[i].comp_graph.all_seq.size();
      //     }
      //     lower_bound_tmp *= set_tree_component[i].comp_graph.num;
      // }
    } else {
      // cout << set_edge_component[i -
      // set_tree_component.size()].comp_graph.num_h << " " <<
      // set_edge_component[i - set_tree_component.size()].comp_graph.num <<
      // endl; if (i == core_height_ind){
      //     if (set_edge_component[i -
      //     set_tree_component.size()].no_feasible_pair && set_edge_component[i
      //     - set_tree_component.size()].with_core_height){
      //         st_max[i] = 1;
      //     } else {
      //         st_max[i] = set_edge_component[i -
      //         set_tree_component.size()].comp_graph.all_seq_h.size();
      //     }
      //     lower_bound_tmp *= set_edge_component[i -
      //     set_tree_component.size()].comp_graph.num_h;
      // } else if (i < core_height_ind){
      if (set_edge_component[i - set_tree_component.size()].no_feasible_pair_loose) {
        st_max[i] = 1;
      } else {
        st_max[i] = set_edge_component[i - set_tree_component.size()]
                        .comp_graph.all_seq_h.size() +
                    set_edge_component[i - set_tree_component.size()]
                        .comp_graph.all_seq.size() + 
                    set_edge_component[i - set_tree_component.size()]
                        .comp_graph.all_seq_h_loose.size() +
                    set_edge_component[i - set_tree_component.size()]
                        .comp_graph.all_seq_loose.size();
      }
      // } else {
      //     if (set_edge_component[i -
      //     set_tree_component.size()].no_feasible_pair &&
      //     !set_edge_component[i -
      //     set_tree_component.size()].with_core_height){
      //         st_max[i] = 1;
      //     } else {
      //         st_max[i] = set_edge_component[i -
      //         set_tree_component.size()].comp_graph.all_seq.size();
      //     }
      //     lower_bound_tmp *= set_edge_component[i -
      //     set_tree_component.size()].comp_graph.num;
      // }
    }
    if (st_max[i] == 0) return;
  }
  // cout << "lower_bound_tmp = " << lower_bound_tmp << endl;

  int h = 0;

  while (h >= 0) {
    if (_stop_gen) return;
    // if (_DEBUG) {
    //     cout << "h = "  << h<< endl;
    //     for (size_t i = 0; i < st.size(); ++i){
    //         cout << st[i] << " " ;
    //     }
    //     cout << endl;
    // }

    if (h == set_tree_component.size() + set_edge_component.size()) {
      prepare_SDF_unfixed_ch_loose(g, coreset, set_tree_component, set_edge_component,
                             g_n, g_m, st, outputfilename, fc_map, fc_all, TS_map);
      --h;
      // if (_DEBUG) cout << "h = " << h<< " end!!!!!"<< endl;
    } else {
      if (st[h] + 1 < st_max[h]) {
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
void Gen_Graph(const Graph& g, const vector<Vertex>& coreset,

               const vector<component>& set_tree_component,
               const vector<component>& set_edge_component, size_t& g_n,
               size_t& g_m,

               size_t& possible_num, size_t& lower_bound,

               string& outputfilename) {
  ofstream output(outputfilename, ios::out);
  output.close();

  size_t n_comp = set_tree_component.size() + set_edge_component.size();
  possible_num = 0;
  lower_bound = 0;

  for (size_t i = 0; i < n_comp; ++i) {
    generate(g, coreset, set_tree_component, set_edge_component, g_n, g_m, i,
             possible_num, lower_bound, outputfilename);
  }
}

// A function to generate SDF format files
void Gen_Graph_unfixed_ch(const Graph& g, const vector<Vertex>& coreset,

                          const vector<component>& set_tree_component,
                          const vector<component>& set_edge_component,
                          size_t& g_n, size_t& g_m,

                          size_t& possible_num, size_t& lower_bound,

                          string& outputfilename) {
  ofstream output(outputfilename, ios::out);
  output.close();

  size_t n_comp = set_tree_component.size() + set_edge_component.size();
  possible_num = 0;
  lower_bound = 0;

  // for (size_t i = 0; i < n_comp; ++i){
  generate_unfixed_ch(g, coreset, set_tree_component, set_edge_component, g_n,
                      g_m, possible_num, lower_bound, outputfilename);
  // }
}

// A function to generate SDF format files
void Gen_Graph_loose(const Graph& g, const vector<Vertex>& coreset,

               const vector<component>& set_tree_component,
               const vector<component>& set_edge_component, size_t& g_n,
               size_t& g_m,

               size_t& possible_num, size_t& lower_bound,

               string& outputfilename,
               const map <size_t, size_t>& fc_map, const vector <unsigned short>& fc_all,
               const map <size_t, TreeSignature>& TS_map) {
  // ofstream output(outputfilename, ios::out);
  // output.close();

  size_t n_comp = set_tree_component.size() + set_edge_component.size();
  // possible_num = 0;
  // lower_bound = 0;

  for (size_t i = 0; i < n_comp; ++i) {
    generate_loose(g, coreset, set_tree_component, set_edge_component, g_n, g_m, i,
             possible_num, lower_bound, outputfilename, fc_map, fc_all, TS_map);
  }
}

// A function to generate SDF format files
void Gen_Graph_unfixed_ch_loose(const Graph& g, const vector<Vertex>& coreset,

                          const vector<component>& set_tree_component,
                          const vector<component>& set_edge_component,
                          size_t& g_n, size_t& g_m,

                          size_t& possible_num, size_t& lower_bound,

                          string& outputfilename,
                          const map <size_t, size_t>& fc_map, const vector <unsigned short>& fc_all,
                          const map <size_t, TreeSignature>& TS_map) {
  ofstream output(outputfilename, ios::out);
  output.close();

  size_t n_comp = set_tree_component.size() + set_edge_component.size();
  // possible_num = 0;
  // lower_bound = 0;

  // for (size_t i = 0; i < n_comp; ++i){
  generate_unfixed_ch_loose(g, coreset, set_tree_component, set_edge_component, g_n,
                      g_m, possible_num, lower_bound, outputfilename, fc_map, fc_all, TS_map);
  // }
}

// A function to get paths from DAG
void Gen_Path_W(const unordered_map<DAG_node, vector<DAG_edge>>& map_DAG,
                const unordered_map<pair<root_status, size_t>,
                                    vector<size_t>>& map_FT_inl,
                const unordered_map<pair<root_status, size_t>,
                                    vector<size_t>>& map_FT_end,
                const DAG_node& node, map_rv& seq, vector<map_rv>& all_seq,
                size_t& global_path_num) {
  if (map_DAG.find(node) == map_DAG.end() || map_DAG.at(node).size() == 0) {
    // cout << "a, ind = " << node.rs.color << ", " << node.ind << endl;
    auto& tmp = map_FT_end.at(make_pair(node.rs, node.ind));
    for (auto& i : tmp) {
      map_rv map_rv_tmp(node.rs, i);
      all_seq.emplace_back(seq, map_rv_tmp, true);

      /* IDO --start-- */

      global_path_num += 1;
      if (global_path_limit != 0 && global_path_num >= global_path_limit) {
        _MEMORY_OUT = true;
        return;
      }

      /* IDO --end-- */
      if (seq_limit != 0 && all_seq.size() >= seq_limit) {
        break;
      }
    }
  } else {
    for (auto& e : map_DAG.at(node)) {
      auto& tmp = map_FT_inl.at(make_pair(e.rs, e.w_ind));
      for (auto& i : tmp) {
        map_rv map_rv_tmp(e.rs, i);
        map_rv new_seq(seq, map_rv_tmp, e.k, true);
        Gen_Path_W(map_DAG, map_FT_inl, map_FT_end, e.node, new_seq, all_seq,
                   global_path_num);

        if (global_path_limit != 0 && global_path_num >= global_path_limit) {
          _MEMORY_OUT = true;
          break;
        }
        if (seq_limit != 0 && all_seq.size() >= seq_limit) {
          break;
        }
      }
      if (global_path_limit != 0 && global_path_num >= global_path_limit) {
        _MEMORY_OUT = true;
        break;
      }
      if (seq_limit != 0 && all_seq.size() >= seq_limit) {
        break;
      }
    }
  }
  return;
}

void Gen_Path_W_path(const unordered_map<DAG_node, vector<DAG_edge>>& map_DAG,
                     const unordered_map<pair<root_status, size_t>,
                                         vector<size_t>>& map_FT_inl,
                     const unordered_map<pair<root_status, size_t>,
                                         vector<size_t>>& map_FT_end,
                     const DAG_node& node, map_rv& seq,
                     vector<map_rv>& all_seq) {
  if (globalTimeKeeper::tt.diff() - _start_time > global_time_limit) {
    _GLOBAL_TIME_OUT = true;
    return;
  }

  // if (getPhysicalMem() - _start_memory > memory_limit) {
  //   _MEMORY_OUT = true;
  //   return;
  // }

  if (map_DAG.find(node) == map_DAG.end() || map_DAG.at(node).size() == 0) {
    // cout << "a, ind = " << node.rs.color << ", " << node.ind << endl;
    map_rv map_rv_tmp(node.rs, node.ind);
    all_seq.emplace_back(seq, map_rv_tmp, true);
  } else {
    for (auto& e : map_DAG.at(node)) {
      auto& tmp = map_FT_inl.at(make_pair(e.rs, e.w_ind));
      map_rv map_rv_tmp(e.rs, e.w_ind);
      map_rv new_seq(seq, map_rv_tmp, e.k, true);
      Gen_Path_W_path(map_DAG, map_FT_inl, map_FT_end, e.node, new_seq,
                      all_seq);
      if (seq_limit != 0 && all_seq.size() >= seq_limit) {
        break;
      }
    }
  }
  return;
}

void Gen_Path_W_path_count(
    const unordered_map<DAG_node, vector<DAG_edge>>& map_DAG,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_inl,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_end,
    const DAG_node& node, map_rv& seq, size_t& all_seq_size) {
  if (globalTimeKeeper::tt.diff() - _start_time > global_time_limit) {
    _GLOBAL_TIME_OUT = true;
    return;
  }

  // if (getPhysicalMem() - _start_memory > memory_limit) {
  //   _MEMORY_OUT = true;
  //   return;
  // }

  if (map_DAG.find(node) == map_DAG.end() || map_DAG.at(node).size() == 0) {
    // cout << "a, ind = " << node.rs.color << ", " << node.ind << endl;
    map_rv map_rv_tmp(node.rs, node.ind);
    all_seq_size += 1;
  } else {
    for (auto& e : map_DAG.at(node)) {
      auto& tmp = map_FT_inl.at(make_pair(e.rs, e.w_ind));
      map_rv map_rv_tmp(e.rs, e.w_ind);
      map_rv new_seq(seq, map_rv_tmp, e.k, true);
      Gen_Path_W_path_count(map_DAG, map_FT_inl, map_FT_end, e.node, new_seq,
                            all_seq_size);
      if (seq_limit != 0 && all_seq_size >= seq_limit) {
        break;
      }
    }
  }
  return;
}

void Gen_Path_W_num(const unordered_map<DAG_node, vector<DAG_edge>>& map_DAG,
                    const unordered_map<pair<root_status, size_t>,
                                        vector<size_t>>& map_FT_inl,
                    const unordered_map<pair<root_status, size_t>,
                                        vector<size_t>>& map_FT_end,
                    const DAG_node& node, size_t num, size_t& total_num) {
  if (globalTimeKeeper::tt.diff() - _start_time > global_time_limit) {
    _GLOBAL_TIME_OUT = true;
    return;
  }

  // if (getPhysicalMem() - _start_memory > memory_limit) {
  //   _MEMORY_OUT = true;
  //   return;
  // }
  if (map_DAG.find(node) == map_DAG.end() || map_DAG.at(node).size() == 0) {
    // cout << "a, ind = " << node.rs.color << ", " << node.ind << endl;
    auto& tmp = map_FT_end.at(make_pair(node.rs, node.ind));
    total_num += num * tmp.size();
  } else {
    for (auto& e : map_DAG.at(node)) {
      auto& tmp = map_FT_inl.at(make_pair(e.rs, e.w_ind));
      Gen_Path_W_num(map_DAG, map_FT_inl, map_FT_end, e.node, num * tmp.size(),
                     total_num);
    }
  }
  return;
}

void Gen_Path_W_in_T(
    const unordered_map<DAG_node, vector<DAG_edge>>& map_DAG,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W1,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W2,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W3,
    const DAG_node& node, map_rv& seq, vector<map_rv>& all_seq,
    size_t& global_path_num) {
  if (map_DAG.find(node) == map_DAG.end() || map_DAG.at(node).size() == 0) {
    // This case should never happen.
    cout << "Error!" << endl;
  } else {
    for (auto& e : map_DAG.at(node)) {
      // cout << "W in T e" << endl;
      auto& tmp = map_FT_W3.at(make_pair(e.rs, e.w_ind));
      for (auto& i : tmp) {
        // cout << "W in T i" << endl;
        map_rv map_rv_tmp(e.rs, i);
        map_rv new_seq(seq, map_rv_tmp, e.k, true);
        // cout << "start Gen_Path_W" << endl;
        Gen_Path_W(map_DAG, map_FT_W2, map_FT_W1, e.node, new_seq, all_seq,
                   global_path_num);
        // cout << "end Gen_Path_W" << endl;
        if (seq_limit != 0 && all_seq.size() >= seq_limit) {
          break;
        }
      }
      if (seq_limit != 0 && all_seq.size() >= seq_limit) {
        break;
      }
    }
  }
  return;
}

void Gen_Path_W_in_T_path(
    const unordered_map<DAG_node, vector<DAG_edge>>& map_DAG,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W1,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W2,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W3,
    const DAG_node& node, map_rv& seq, vector<map_rv>& all_seq) {
  if (globalTimeKeeper::tt.diff() - _start_time > global_time_limit) {
    _GLOBAL_TIME_OUT = true;
    return;
  }

  // if (getPhysicalMem() - _start_memory > memory_limit) {
  //   _MEMORY_OUT = true;
  //   return;
  // }
  if (map_DAG.find(node) == map_DAG.end() || map_DAG.at(node).size() == 0) {
    // This case should never happen.
    cout << "Error!" << endl;
  } else {
    for (auto& e : map_DAG.at(node)) {
      // cout << "W in T e" << endl;

      map_rv map_rv_tmp(e.rs, e.w_ind);
      map_rv new_seq(seq, map_rv_tmp, e.k, true);
      // cout << "start Gen_Path_W" << endl;
      Gen_Path_W_path(map_DAG, map_FT_W2, map_FT_W1, e.node, new_seq, all_seq);

      if (seq_limit != 0 && all_seq.size() >= seq_limit) {
        break;
      }
    }
  }
  return;
}

void Gen_Path_W_in_T_path_count(
    const unordered_map<DAG_node, vector<DAG_edge>>& map_DAG,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W1,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W2,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W3,
    const DAG_node& node, map_rv& seq, size_t& all_seq_size) {
  if (globalTimeKeeper::tt.diff() - _start_time > global_time_limit) {
    _GLOBAL_TIME_OUT = true;
    return;
  }

  // if (getPhysicalMem() - _start_memory > memory_limit) {
  //   _MEMORY_OUT = true;
  //   return;
  // }

  if (map_DAG.find(node) == map_DAG.end() || map_DAG.at(node).size() == 0) {
    // This case should never happen.
    cout << "Error!" << endl;
  } else {
    for (auto& e : map_DAG.at(node)) {
      // cout << "W in T e" << endl;

      map_rv map_rv_tmp(e.rs, e.w_ind);
      map_rv new_seq(seq, map_rv_tmp, e.k, true);
      // cout << "start Gen_Path_W" << endl;
      Gen_Path_W_path_count(map_DAG, map_FT_W2, map_FT_W1, e.node, new_seq,
                            all_seq_size);

      if (seq_limit != 0 && all_seq_size >= seq_limit) {
        break;
      }
    }
  }
  return;
}

void Gen_Path_W_in_T_num(
    const unordered_map<DAG_node, vector<DAG_edge>>& map_DAG,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W1,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W2,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W3,
    const DAG_node& node, size_t num, size_t& total_num) {
  if (globalTimeKeeper::tt.diff() - _start_time > global_time_limit) {
    _GLOBAL_TIME_OUT = true;
    return;
  }

  // if (getPhysicalMem() - _start_memory > memory_limit) {
  //   _MEMORY_OUT = true;
  //   return;
  // }
  if (map_DAG.find(node) == map_DAG.end() || map_DAG.at(node).size() == 0) {
    // This case should never happen.
    cout << "Error!" << endl;
  } else {
    for (auto& e : map_DAG.at(node)) {
      // cout << "W in T e" << endl;
      auto& tmp = map_FT_W3.at(make_pair(e.rs, e.w_ind));
      Gen_Path_W_num(map_DAG, map_FT_W2, map_FT_W1, e.node, num * tmp.size(),
                     total_num);
    }
  }
  return;
}

void Gen_Path_T(
    const unordered_map<DAG_node_T, vector<DAG_edge_T>>& map_DAG_T,
    const unordered_map<DAG_node, vector<DAG_edge>>& map_DAG_W,

    const unordered_map<root_status_height,
                        map<size_t, vector<DAG_node_height>>>& map_T2_W3,

    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W1,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W2,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W3,
    const unordered_map<pair<root_status_height, size_t>, vector<size_t>>&
        map_FT_W2_core,
    const DAG_node_T& node_T, map_tree& seq, vector<map_tree>& all_seq,
    size_t& global_path_num) {
  if (map_DAG_T.find(node_T) == map_DAG_T.end() ||
      map_DAG_T.at(node_T).size() == 0) {
    map_rv tmp(node_T.rsT.first, 0);
    map_tree mt_tmp(seq, map_tree(tmp), seq.height, true);
    all_seq.emplace_back(mt_tmp, false);

    global_path_num += 1;

    if (global_path_limit != 0 && global_path_num >= global_path_limit) {
      _MEMORY_OUT = true;
      return;
    }

  } else {
    vector<int> e_ind(map_DAG_T.at(node_T).size());
    iota(e_ind.begin(), e_ind.end(), 0);
    shuffle(e_ind.begin(), e_ind.end(), engine);

    for (auto& _e_ind : e_ind) {
      // for (auto& e : map_DAG_T.at(node_T)){
      auto& e = map_DAG_T.at(node_T)[_e_ind];
      // cout << "e" << endl;
      auto& node_vector = map_T2_W3.at(e.rsh).at(e.T2_ind);
      for (auto& node_tmp : node_vector) {
        // cout << "node_tmp.l = " << node_tmp.l << endl;
        vector<map_rv> all_seq_W1;
        all_seq_W1.clear();
        map_rv seq_tmp;
        if (node_tmp.l == 0) {
          auto& tmp =
              map_FT_W2_core.at(make_pair(node_tmp.rsh, node_tmp.ind));
          for (auto& i : tmp) {
            // cout << "i" << endl;
            map_rv map_rv_tmp(node_tmp.rsh, i);
            all_seq_W1.push_back(map_rv_tmp);

            if (seq_limit != 0 && all_seq_W1.size() >= seq_limit) {
              break;
            }
          }
          for (auto& seq_i : all_seq_W1) {
            // cout << "seq_i" << endl;
            map_tree new_seq(seq, map_tree(seq_i), e.k,
                             max(seq.height, node_tmp.l), true);
            Gen_Path_T(map_DAG_T, map_DAG_W, map_T2_W3, map_FT_W1, map_FT_W2,
                       map_FT_W3, map_FT_W2_core, e.node, new_seq, all_seq,
                       global_path_num);
            if (seq_limit != 0 && all_seq.size() >= seq_limit) {
              break;
            }
          }
        } else {
          root_status rs_tmp(node_tmp.rsh.color, node_tmp.rsh.val,
                             node_tmp.rsh.deg);
          DAG_node DAG_node_tmp(rs_tmp, node_tmp.l, node_tmp.ind);
          // cout << "start Gen_Path_W_in_T" << endl;
          Gen_Path_W_in_T(map_DAG_W, map_FT_W1, map_FT_W2, map_FT_W3,
                          DAG_node_tmp, seq_tmp, all_seq_W1, global_path_num);
          // cout << "end Gen_Path_W_in_T" << endl;
          for (auto& seq_i : all_seq_W1) {
            // cout << "seq_i" << endl;
            map_tree new_seq(seq, map_tree(seq_i), e.k,
                             max(seq.height, node_tmp.l), true);
            Gen_Path_T(map_DAG_T, map_DAG_W, map_T2_W3, map_FT_W1, map_FT_W2,
                       map_FT_W3, map_FT_W2_core, e.node, new_seq, all_seq,
                       global_path_num);
            if (global_path_limit != 0 &&
                global_path_num >= global_path_limit) {
              _MEMORY_OUT = true;
              break;
            }
            if (seq_limit != 0 && all_seq.size() >= seq_limit) {
              break;
            }
          }
        }
        if (global_path_limit != 0 && global_path_num >= global_path_limit) {
          _MEMORY_OUT = true;
          break;
        }
        if (seq_limit != 0 && all_seq.size() >= seq_limit) {
          break;
        }
      }
      if (global_path_limit != 0 && global_path_num >= global_path_limit) {
        _MEMORY_OUT = true;
        break;
      }
      if (seq_limit != 0 && all_seq.size() >= seq_limit) {
        break;
      }
    }
  }
  return;
}

void Gen_Path_T_only_path(
    const unordered_map<DAG_node_T, vector<DAG_edge_T>>& map_DAG_T,
    const unordered_map<DAG_node, vector<DAG_edge>>& map_DAG_W,

    const unordered_map<root_status_height,
                        map<size_t, vector<DAG_node_height>>>& map_T2_W3,

    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W1,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W2,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W3,
    const unordered_map<pair<root_status_height, size_t>, vector<size_t>>&
        map_FT_W2_core,
    const DAG_node_T& node_T, map_tree& seq, vector<map_tree>& all_seq) {
  if (globalTimeKeeper::tt.diff() - _start_time > global_time_limit) {
    _GLOBAL_TIME_OUT = true;
    return;
  }

  // if (getPhysicalMem() - _start_memory > memory_limit) {
  //   _MEMORY_OUT = true;
  //   return;
  // }
  if (map_DAG_T.find(node_T) == map_DAG_T.end() ||
      map_DAG_T.at(node_T).size() == 0) {
    map_rv tmp(node_T.rsT.first, 0);
    map_tree mt_tmp(seq, map_tree(tmp), seq.height, true);
    all_seq.emplace_back(mt_tmp, false);
  } else {
    for (auto& e : map_DAG_T.at(node_T)) {
      // cout << "e" << endl;
      auto& node_vector = map_T2_W3.at(e.rsh).at(e.T2_ind);
      // for (auto& node_tmp : node_vector){
      //     // cout << "node_tmp.l = " << node_tmp.l << endl;
      //     vector <map_rv> all_seq_W1;
      //     all_seq_W1.clear();
      //     map_rv seq_tmp;
      //     if (node_tmp.l == 0){
      //         map_rv map_rv_tmp(node_tmp.rsh.color, node_tmp.ind);
      //         all_seq_W1.push_back(map_rv_tmp);

      //         for (auto& seq_i : all_seq_W1){
      //             // cout << "seq_i" << endl;
      //             map_tree new_seq(seq, map_tree(seq_i), e.k, max(seq.height,
      //             node_tmp.l), true); Gen_Path_T_path(map_DAG_T, map_DAG_W,
      //             map_T2_W3, map_FT_W1, map_FT_W2, map_FT_W3, map_FT_W2_core,
      //             e.node, new_seq, all_seq); if (seq_limit != 0 &&
      //             all_seq.size() >= seq_limit){
      //                 break;
      //             }
      //         }
      //     } else {
      //         root_status rs_tmp(node_tmp.rsh.color, node_tmp.rsh.val,
      //         node_tmp.rsh.deg); DAG_node DAG_node_tmp(rs_tmp, node_tmp.l,
      //         node_tmp.ind);
      //         // cout << "start Gen_Path_W_in_T" << endl;
      //         Gen_Path_W_in_T_path(map_DAG_W, map_FT_W1, map_FT_W2,
      //         map_FT_W3, DAG_node_tmp, seq_tmp, all_seq_W1);
      //         // cout << "end Gen_Path_W_in_T" << endl;
      //         for (auto& seq_i : all_seq_W1){
      //             // cout << "seq_i" << endl;
      //             map_tree new_seq(seq, map_tree(seq_i), e.k, max(seq.height,
      //             node_tmp.l), true); Gen_Path_T_path(map_DAG_T, map_DAG_W,
      //             map_T2_W3, map_FT_W1, map_FT_W2, map_FT_W3, map_FT_W2_core,
      //             e.node, new_seq, all_seq); if (seq_limit != 0 &&
      //             all_seq.size() >= seq_limit){
      //                 break;
      //             }
      //         }
      //     }
      //     if (seq_limit != 0 && all_seq.size() >= seq_limit){
      //         break;
      //     }
      // }

      if (seq_limit != 0 && all_seq.size() >= seq_limit) {
        break;
      }

      map_tree tmp;
      map_tree new_seq(seq, tmp, e.k, seq.height, true);
      Gen_Path_T_only_path(map_DAG_T, map_DAG_W, map_T2_W3, map_FT_W1,
                           map_FT_W2, map_FT_W3, map_FT_W2_core, e.node,
                           new_seq, all_seq);
    }
  }
  return;
}

void Gen_Path_T_path(
    const unordered_map<DAG_node_T, vector<DAG_edge_T>>& map_DAG_T,
    const unordered_map<DAG_node, vector<DAG_edge>>& map_DAG_W,

    const unordered_map<root_status_height,
                        map<size_t, vector<DAG_node_height>>>& map_T2_W3,

    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W1,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W2,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W3,
    const unordered_map<pair<root_status_height, size_t>, vector<size_t>>&
        map_FT_W2_core,
    const DAG_node_T& node_T, map_tree& seq, vector<map_tree>& all_seq) {
  if (globalTimeKeeper::tt.diff() - _start_time > global_time_limit) {
    _GLOBAL_TIME_OUT = true;
    return;
  }

  // if (getPhysicalMem() - _start_memory > memory_limit) {
  //   _MEMORY_OUT = true;
  //   return;
  // }
  if (map_DAG_T.find(node_T) == map_DAG_T.end() ||
      map_DAG_T.at(node_T).size() == 0) {
    map_rv tmp(node_T.rsT.first, 0);
    map_tree mt_tmp(seq, map_tree(tmp), seq.height, true);
    all_seq.emplace_back(mt_tmp, false);
  } else {
    for (auto& e : map_DAG_T.at(node_T)) {
      // cout << "e" << endl;
      auto& node_vector = map_T2_W3.at(e.rsh).at(e.T2_ind);
      for (auto& node_tmp : node_vector) {
        // cout << "node_tmp.l = " << node_tmp.l << endl;
        vector<map_rv> all_seq_W1;
        all_seq_W1.clear();
        map_rv seq_tmp;
        if (node_tmp.l == 0) {
          map_rv map_rv_tmp(node_tmp.rsh, node_tmp.ind);
          all_seq_W1.push_back(map_rv_tmp);

          for (auto& seq_i : all_seq_W1) {
            // cout << "seq_i" << endl;
            map_tree new_seq(seq, map_tree(seq_i), e.k,
                             max(seq.height, node_tmp.l), true);
            Gen_Path_T_path(map_DAG_T, map_DAG_W, map_T2_W3, map_FT_W1,
                            map_FT_W2, map_FT_W3, map_FT_W2_core, e.node,
                            new_seq, all_seq);
            if (seq_limit != 0 && all_seq.size() >= seq_limit) {
              break;
            }
          }
        } else {
          root_status rs_tmp(node_tmp.rsh.color, node_tmp.rsh.val,
                             node_tmp.rsh.deg);
          DAG_node DAG_node_tmp(rs_tmp, node_tmp.l, node_tmp.ind);
          // cout << "start Gen_Path_W_in_T" << endl;
          Gen_Path_W_in_T_path(map_DAG_W, map_FT_W1, map_FT_W2, map_FT_W3,
                               DAG_node_tmp, seq_tmp, all_seq_W1);
          // cout << "end Gen_Path_W_in_T" << endl;
          for (auto& seq_i : all_seq_W1) {
            // cout << "seq_i" << endl;
            map_tree new_seq(seq, map_tree(seq_i), e.k,
                             max(seq.height, node_tmp.l), true);
            Gen_Path_T_path(map_DAG_T, map_DAG_W, map_T2_W3, map_FT_W1,
                            map_FT_W2, map_FT_W3, map_FT_W2_core, e.node,
                            new_seq, all_seq);
            if (seq_limit != 0 && all_seq.size() >= seq_limit) {
              break;
            }
          }
        }
        if (seq_limit != 0 && all_seq.size() >= seq_limit) {
          break;
        }
      }
      if (seq_limit != 0 && all_seq.size() >= seq_limit) {
        break;
      }
    }
  }
  return;
}

void Gen_Path_T_num(
    const unordered_map<DAG_node_T, vector<DAG_edge_T>>& map_DAG_T,
    const unordered_map<DAG_node, vector<DAG_edge>>& map_DAG_W,

    const unordered_map<root_status_height,
                        map<size_t, vector<DAG_node_height>>>& map_T2_W3,

    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W1,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W2,
    const unordered_map<pair<root_status, size_t>, vector<size_t>>&
        map_FT_W3,
    const unordered_map<pair<root_status_height, size_t>, vector<size_t>>&
        map_FT_W2_core,
    const DAG_node_T& node_T, size_t num, size_t& total_num) {
  if (globalTimeKeeper::tt.diff() - _start_time > global_time_limit) {
    _GLOBAL_TIME_OUT = true;
    return;
  }

  // if (getPhysicalMem() - _start_memory > memory_limit) {
  //   _MEMORY_OUT = true;
  //   return;
  // }
  if (map_DAG_T.find(node_T) == map_DAG_T.end() ||
      map_DAG_T.at(node_T).size() == 0) {
    total_num += num;
  } else {
    for (auto& e : map_DAG_T.at(node_T)) {
      // cout << "e" << endl;
      auto& node_vector = map_T2_W3.at(e.rsh).at(e.T2_ind);
      for (auto& node_tmp : node_vector) {
        // cout << "node_tmp.l = " << node_tmp.l << endl;
        vector<map_rv> all_seq_W1;
        all_seq_W1.clear();
        map_rv seq_tmp;
        if (node_tmp.l == 0) {
          auto& tmp =
              map_FT_W2_core.at(make_pair(node_tmp.rsh, node_tmp.ind));
          Gen_Path_T_num(map_DAG_T, map_DAG_W, map_T2_W3, map_FT_W1, map_FT_W2,
                         map_FT_W3, map_FT_W2_core, e.node, num * tmp.size(),
                         total_num);
        } else {
          root_status rs_tmp(node_tmp.rsh.color, node_tmp.rsh.val,
                             node_tmp.rsh.deg);
          DAG_node DAG_node_tmp(rs_tmp, node_tmp.l, node_tmp.ind);
          // cout << "start Gen_Path_W_in_T" << endl;
          size_t num_tmp = 0;
          Gen_Path_W_in_T_num(map_DAG_W, map_FT_W1, map_FT_W2, map_FT_W3,
                              DAG_node_tmp, 1, num_tmp);
          // cout << "end Gen_Path_W_in_T" << endl;
          Gen_Path_T_num(map_DAG_T, map_DAG_W, map_T2_W3, map_FT_W1, map_FT_W2,
                         map_FT_W3, map_FT_W2_core, e.node, num * num_tmp,
                         total_num);
        }
      }
    }
  }
  return;
}


#endif
