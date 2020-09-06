/*===========================================================================
Feature vector generation
Author : Discrete Mathematics Lab

	This file implements functions that given an sdf file
	of acyclic molecular graphs, calculates the feature vector
	for each one, and stores the values as a comma separated text file,
	whose first column gives the descriptor name.

    LICENSE: MIT license, see License.txt

	* Version 2020/04/08: uses bond configuration as descriptor
	* Version 2020/07/14: separates descriptor between internal/external
	* Version 2020/08/24: includes branch-leaf-number descriptor
  ===========================================================================*/

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <stack>
#include <queue>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <algorithm>

bool debug = false; // Used for debugging during development

using namespace std;
typedef size_t Vertex;

/**
 * Write hash functions for vectors and pairs
 * and include them into the std namespace for convenience
 */
namespace std{
	/**
	 * A function that given a seed hash value
	 * and an object v, combines the seed value
	 * and the hash of the object.
	 * Follows Boost's hash_combine function
	 */
	template <class T>
	size_t hash_combine(std::size_t& seed, const T& v){
		std::hash<T> hasher;
		size_t ans = seed;
		ans ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
		return ans;
	}

    /**
     * A hash function for std::vector
     * Note, the elements of the vector must be hashable in std.
     */
    template <class T>
    class hash <vector <T>>{
    public:
        size_t operator()(const vector <T>& aa) const {
            size_t res = 0;
            for (const auto &t : aa){
                res = hash_combine<T>(res, t);
            }
            return res;
        }
    };

	/**
	 * A hash function for std::pair
	 * Note, each member of the pair must be hashable in std.
	 */
    template <class T1, class T2>
    class hash <pair<T1, T2>>{
    public:
        size_t operator()(const pair<T1, T2>& pp) const {
            size_t res = std::hash<T1>{} (pp.first);
            res = hash_combine<T2>(res, pp.second);
            return res;
        }
    };
} // Namespace std

/**
 * A simple structure to store a chemical graph object,
 * together with vertex colors (element names) and bond multiplicities.
 * Graphs are stored in adjacency lists, and bond multiplicities in a 2D matrix
 */
struct Graph{
	string CID;							// CID
	size_t numAtom, numBond;			// number of vertices and edges
	vector < vector < Vertex > > adj;	// adjacent vertices
	vector < string > alpha;			// atom types (C, H, O, ...)
	vector < vector < size_t > > beta;	// bond types (1, 2, ...)
	vector < bool > status;				// used in calculating branch height
};

/**
 * A helper encapsulating structure to handle rooted trees
 */
struct RootedTree{
	Graph graph;
	Vertex root;

};

// Typedef a map to store the mass of each type of atom
using MassMap = unordered_map<string, size_t>;

/**
 * A function to initialize a map tying an element's symbol
 * to its atomic mass.
 * @return
 */
MassMap
init_MassMap() {
	MassMap M;
	// Initialize entries in the map manually
	M["B"]  = 108;
	M["C"]  = 120;
	M["O"]  = 160;
	M["N"]  = 140;
	M["F"]  = 190;
	M["Si"] = 280;
	M["P"]  = 310;
	M["S"]  = 320;
	M["Cl"] = 355;
	M["V"]  = 510;
	M["Br"] = 800;
	M["Cd"] = 1124;
	M["I"]  = 1270;
	M["Hg"] = 2006;
	M["Pb"] = 2072;
    M["Al"] = 269;
	// Finally return the the initialized map
	return M;
}


// Typedef a map to store the valence of each type of atom
using ValenceMap = unordered_map<string, size_t>;


/**
 * A function to initialize a map tying an element's symbol
 * to its most common valence.
 * @return
 */
ValenceMap
init_ValenceMap() {
	ValenceMap M;
	// Initialize entries in the map manually
	// Check https://sciencenotes.org/valences-of-the-elements/
	M["B"]  = 3;
	M["C"]  = 4;
	M["O"]  = 2;
	M["N"]  = 3;
	M["F"]  = 1;
	M["Si"] = 4;
	M["P"]  = 3;
	M["S"]  = 2;
	M["Cl"] = 1;
	M["V"]  = 3;
	M["Br"] = 1;
	M["Cd"] = 2;
	M["I"]  = 1;
	M["Hg"] = 2;
	M["Pb"] = 2;
    M["Al"] = 3;
	// Finally return the the initialized map
	return M;
}


// The set Bc, used to help calculating the number of bond-configuration
vector <vector <size_t>>
set_Bc = {	{1, 2, 1},
			{1, 2, 2},
			{1, 2, 3},
			{1, 3, 1},
			{1, 3, 2},
			{1, 4, 1},
			{2, 2, 1},
			{2, 2, 2},
			{2, 2, 3},
			{2, 3, 1},
			{2, 3, 2},
			{2, 4, 1},
			{3, 3, 1},
			{3, 3, 2},
			{3, 4, 1},
			{4, 4, 1}
		};

/**
 * Read a collection of chemical graphs from an sdf file
 * @param string inFileName: the filename of the input sd file
 * @return a vector storing the chemical graphs
 */
vector < Graph >
read_graph_sdf(const string& inFileName) {

	vector < Graph > read_graphs; // return value
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

	if (debug) {
		cout << endl << "----- Input file " << inFileName << " -----" << endl;
	}

	//read the file
	int count = -1;
	int flag = 0;
	string line;

	while (getline(infile, line)) {

		if (flag == 0) {
			Graph g;
			read_graphs.push_back(g);
			count++;
			flag = 1;
		}

		while (flag == 1) {
			read_graphs[count].CID = line;
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

			read_graphs[count].numAtom = n;
			read_graphs[count].numBond = m;

			read_graphs[count].alpha.resize(n);
			read_graphs[count].status.resize(n);
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
				if (debug)
					cout << "3D Data, alpha: " << tdData
						 << "; " << alpha << endl;
				read_graphs[count].alpha[i] = alpha;
				read_graphs[count].status[i] = true;
				delete[] tdData;
			}

			// initialize a matrix for storing bond multiplicities
			read_graphs[count].beta.resize(n);
			for (int i = 0; i < n; ++i) {
				read_graphs[count].beta[i].resize(n);
				for (int j = 0; j < n; ++j) {
					read_graphs[count].beta[i][j] = 0;
				}
			}

			for (int i = 0; i < n; ++i) {
				vector < Vertex > tmp;
				read_graphs[count].adj.push_back(tmp);
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
				read_graphs[count].adj[v1-1].push_back(v2-1);
				read_graphs[count].adj[v2-1].push_back(v1-1);
				read_graphs[count].beta[v1-1][v2-1] = mul;
				read_graphs[count].beta[v2-1][v1-1] = mul;
			}
			flag = 2;
		}
		// The end of the graph information in the file is marked with a "$$$$"
		if (line == "$$$$" && flag == 2) {
			flag = 0;
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


/**
 * Get the number of each element type in the molecular graph
 * @param Graph g: chemical graph
 * @param element: map of element symbols and integer ids
 * @param Graph globalElement: store all element symbols in a vector
 */
void
getElementCount( const Graph &g,
		    	 map <string, size_t> &element_in,
		    	 map <string, size_t> &element_ex,
				 vector < string > &globalElement,
				 vector < bool > &_status ) {
	// Assign carbon to the first order
	element_in["C"] = 0;
	element_ex["C"] = 0;
	for (size_t i = 0; i < g.alpha.size(); ++i) {
		if (_status[i]){
			if (element_in.find(g.alpha[i]) == element_in.end()) {
				element_in[g.alpha[i]] = 1;
				auto it = find( globalElement.begin(),
								globalElement.end(),
								g.alpha[i] );
				if (it == globalElement.end()) {
					globalElement.push_back(g.alpha[i]);
				}
			} else {
				++element_in[g.alpha[i]];
			}
		} else{
			if (element_ex.find(g.alpha[i]) == element_ex.end()) {
				element_ex[g.alpha[i]] = 1;
				auto it = find( globalElement.begin(),
								globalElement.end(),
								g.alpha[i] );
				if (it == globalElement.end()) {
					globalElement.push_back(g.alpha[i]);
				}
			} else {
				++element_ex[g.alpha[i]];
			}
		}
	}
	// print element map
	if (debug) {
		for (const auto &a : element_in) {
			cout << "Element (in): "<< a.first << " " << a.second << endl;
		}
		for (const auto &a : element_ex) {
			cout << "Element (ex): "<< a.first << " " << a.second << endl;
		}
	}
}


/**
 * Assign a total order on the elements
 * @param Graph globalElement
 * @param elementOrder
 */
void assignElemOrder( const vector < string > &globalElement,
					  map<string, size_t> &elementOrder) {
	int order = 0;
	for (size_t i = 0; i < globalElement.size(); ++i) {
		if (elementOrder.find(globalElement[i]) == elementOrder.end()) {
			elementOrder[globalElement[i]] = order;
			order++;
		}
	}
}


/**
 * Get the edges in the graph as element-bond_mutliplicity-element string
 * @param Graph g: input graph
 * @param elementOrder: a pre-assigned total order on the elements
 * @param onepath_map: store the number of occurences of each edge
 * @param Graph globalOnePath: store all edges in a vector
 */
void
getOnePath( const Graph &g,
			map <string, size_t> &elementOrder,
			map <string, size_t> &onepath_map_in,
			map <string, size_t> &onepath_map_ex,
			vector <string> &globalOnePath,
			vector <bool>& _status) {

	vector < vector < int > > onePath;
	vector < string > onePath_code;

	for (size_t i = 0; i < g.numAtom; ++i) {
		for (size_t j = i; j < g.numAtom; ++j) {
			if (g.beta[i][j] != 0) {
				vector < int > temp;
				temp.push_back(i);
				temp.push_back(j);
				temp.push_back(g.beta[i][j]);
				onePath.push_back(temp);
			}
		}
	}

	for (size_t i = 0; i < onePath.size(); ++i) {
		string label1 = g.alpha[onePath[i][0]];
		string label2 = g.alpha[onePath[i][1]];
		string bondType = to_string(onePath[i][2]);

		if (elementOrder[label1] <= elementOrder[label2]) {
			string temp;
			temp.append(label1);
			temp.append(bondType);
			temp.append(label2);
			onePath_code.push_back(temp);
		} else {
			string temp;
			temp.append(label2);
			temp.append(bondType);
			temp.append(label1);
			onePath_code.push_back(temp);
		}
	}

	for (size_t i = 0; i < onePath_code.size(); ++i) {
		if (_status[onePath[i][0]] && _status[onePath[i][1]]){
			if (onepath_map_in.find(onePath_code[i]) == onepath_map_in.end()) {
				onepath_map_in[onePath_code[i]] = 1;
				auto it = find( globalOnePath.begin(),
								globalOnePath.end(),
								onePath_code[i]);
				if (it == globalOnePath.end()) {
					globalOnePath.push_back(onePath_code[i]);
				}
			} else {
				onepath_map_in[onePath_code[i]]++;
			}
		} else {
			if (onepath_map_ex.find(onePath_code[i]) == onepath_map_ex.end()) {
				onepath_map_ex[onePath_code[i]] = 1;
				auto it = find( globalOnePath.begin(),
								globalOnePath.end(),
								onePath_code[i]);
				if (it == globalOnePath.end()) {
					globalOnePath.push_back(onePath_code[i]);
				}
			} else {
				onepath_map_ex[onePath_code[i]]++;
			}
		}
	}
	// print the map of the frequency of each edge
	if (debug) {
		cout << "Path map (in): " << endl;
		for (auto a : onepath_map_in) {
			cout << "Path: " << a.first << " " << a.second << endl;
		}
		cout << "Path map (ex): " << endl;
		for (auto a : onepath_map_ex) {
			cout << "Path: " << a.first << " " << a.second << endl;
		}
	}
}


/**
 * Calculate the number of vertices of degree d, d = 1, 2, ..., 6
 * @param Graph g: input graph
 * @param degree_vertices: store the numbers in a vector
 */
void calcDegree (const Graph &g, 
	vector <size_t > &degree_vertices_in,
	vector <size_t > &degree_vertices_ex,
	vector <bool > &_status) {
	degree_vertices_in = vector < size_t > (6, 0);
	degree_vertices_ex = vector < size_t > (6, 0);
	for (size_t i = 0; i < 6; ++i) {
		degree_vertices_in[i] = 0;
		degree_vertices_ex[i] = 0;
		for (size_t j = 0; j < g.numAtom; ++j) {
			if (g.adj[j].size() == i+1){
				if (_status[j]){
					degree_vertices_in[i] += 1;
				} else {
					degree_vertices_ex[i] += 1;
				}

			}
		}
	}
}

/**
 * Calculate the number of edges of bond-configuration 
 * @param Graph g: input graph
 * @param degree_vertices: store the numbers in a map
 */
void 
calcBondConfiguration (const Graph &g, 
	map <vector <size_t>, size_t> &bond_conf_in,
	map <vector <size_t>, size_t> &bond_conf_ex,
	vector <bool> &_status
	) {
	bond_conf_in.clear();
	bond_conf_ex.clear();
	for (size_t u = 0; u < g.numAtom; ++u) {
        for (auto v : g.adj[u]){
            if (u < v){
                size_t d1 = g.adj[u].size();
                size_t d2 = g.adj[v].size();
                size_t k = g.beta[u][v];
                
                if (_status[u] && _status[v]){
                	if (d1 >= d2){
	                    bond_conf_in[{d2, d1, k}]++;
	                } else {
	                    bond_conf_in[{d1, d2, k}]++;
	                }
                } else {
                	if (d1 >= d2){
	                    bond_conf_ex[{d2, d1, k}]++;
	                } else {
	                    bond_conf_ex[{d1, d2, k}]++;
	                }
                }
            }
        }             
   }
}


/**
 * Calculate the numbers of double and triple bonds in a graph g
 * @param Graph g: input chemical graph
 * @param double_bond: stores the number of double bonds
 * @param triple_bond: stores the number of triple bonds
 */
void
calcDoubleTripleBond( const Graph &g,
					  size_t &double_bond_in,
					  size_t &double_bond_ex,
					  size_t &triple_bond_in,
					  size_t &triple_bond_ex,
					  vector < bool > &_status) {
	double_bond_in = 0;
	double_bond_ex = 0;
	triple_bond_in = 0;
	triple_bond_ex = 0;	

	for (size_t i = 0; i < g.numAtom; ++i) {
		for (size_t j = i; j < g.numAtom; ++j) {
			if (g.beta[i][j] == 2 && _status[i] && _status[j])
				double_bond_in += 1;
			if (g.beta[i][j] == 2 && (!_status[i] || !_status[j]))
				double_bond_ex += 1;
			if (g.beta[i][j] == 3 && _status[i] && _status[j])
				triple_bond_in += 1;
			if (g.beta[i][j] == 3 && (!_status[i] || !_status[j]))
				triple_bond_ex += 1;
		}
	}
}


/**
 * Calculate the total atomic mass of a molecule
 * @param Graph g: input chemical graph
 * @param MassMap mass_map: input map element symbol<>atomic mass
 * @return the sum of atomic masses over the graph
 */
double calcM( const Graph &g,
			  const MassMap &mass_map ) {
	int mass = 0;
	for (size_t i = 0; i < g.numAtom; ++i) {
		try {
			mass += mass_map.at(g.alpha[i]);
		} catch (const exception &e) {
			cerr << "Key error: " << g.alpha[i] << endl;
			throw(e);
		}
	}
	return mass;
}


int calcEffectiveVertexNum(const Graph& g){
	int numv = 0;
	for (Vertex u = 0; u < g.numAtom; ++u) {
		if (g.status[u] == true) {
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
	return nume;
}


/**
 * A helper function in calculating branch height/number
 * Modifies the input graph and removes all leaf edges.
 */
void cutLeaf_oneLayer(Graph &g){
    if (calcEffectiveVertexNum(g) <= 2) {
        return;
    }
    auto _adj = g.adj;
	for (Vertex u = 0; u < g.numAtom; ++u) {
		if (_adj[u].size() == 1){
			for (auto v : g.adj[u]) {
				if (debug) {
					cerr << "#### leaf :" << u << " ####" << endl;
					cerr << "g.adj[" << u << "]: ";
					for (auto nn :  g.adj[u])
						cout << nn << "; ";
					cout << endl;
				}
				g.adj[u].pop_back();
				if (debug) {
					cerr << "#### leaf popped :" << u << " ####" << endl;
					cerr << "g.adj[" << u << "]: ";
					for (auto nn :  g.adj[u])
						cout << nn << "; ";
					cout << endl;
				}
				if (debug) {
					cerr << "#### leaf neighbor :" << v << " ####" << endl;
					cerr << "g.adj[" << v << "]: ";
					for (auto nn :  g.adj[v])
						cout << nn << "; ";
					cout << endl;
				}
				g.adj[v].erase(
						remove(g.adj[v].begin(), g.adj[v].end(), u),
						g.adj[v].end()
						);
				if (debug) {
					cerr << "#### leaf neighbor remove :" <<
							v << " ####" << endl;
					cerr << "g.adj[" << v << "]: ";
					for (auto nn :  g.adj[v])
						cout << nn << "; ";
					cout << endl;
				}
			}
			g.status[u] = false;
		}
	}
}

void calcCoreVertexSet(const Graph& g, vector <Vertex> &coreset){

	Graph _g = g;
	bool loop = true;
	while (loop) {
		int numv = calcEffectiveVertexNum(_g);
		if (numv > 2){
			cutLeaf_oneLayer(_g);
			int numv2 = calcEffectiveVertexNum(_g);
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


/**
 * A helper function in calculating branch height/width
 * Modifies the supplied graph by pruning leafs up to height k
 * @param g
 * @param k
 */
void cutTree(Graph &g, const size_t& k){
	for (size_t i = 0; i < k; ++i) {
		cutLeaf_oneLayer(g);
	}
}


/**
 * A helper function in calculating branch height/number
 */
void searchBranch(
		const Graph& g,
		const Vertex& u,
		vector < Vertex >& forbiddenset,
		RootedTree &branch
		){

    for (auto v : g.adj[u]) {
    	auto itr = find(forbiddenset.begin(), forbiddenset.end(), v);
    	if (itr == forbiddenset.end()) {
    		branch.graph.adj[u].push_back(v);
            branch.graph.adj[v].push_back(u);
            branch.graph.status[v] = true;
    		auto forbiddenset_new = forbiddenset;
    		// push v at the end of the vector forbiddenset
    		forbiddenset.push_back(v);
    		searchBranch(g, v, forbiddenset, branch);
    		// restore vector forbiddenset after recursion
    		forbiddenset.pop_back();
    	}
    }
}

void getCoreTree(
		const Graph& g,
		vector <RootedTree> &coretree,
		const vector <Vertex> &coreset
		){
    for (auto u : coreset) {
	    for (auto v : g.adj[u]) {
    		auto itr = find(coreset.begin(), coreset.end(), v);
    		if (itr == coreset.end()) {
    			RootedTree subtree;
    			subtree.graph = g;
				for (Vertex uv = 0; uv < g.numAtom; ++uv) {
					subtree.graph.adj[uv].clear();
					subtree.graph.status[uv] = false;
	    		}
	    		subtree.root = u;
	    		subtree.graph.status[u] = true;
	    		subtree.graph.status[v] = true;
    			subtree.graph.adj[u].push_back(v);
    			subtree.graph.adj[v].push_back(u);
    			auto forbiddenset = coreset;
    			forbiddenset.push_back(v);
    			searchBranch(g, v, forbiddenset, subtree);
    			coretree.push_back(subtree);
    		}
    	}
    }
}


/**
 * A helper function in calculating branch-height
 * 	 and branch-leaf-number
 * @param branch - a rooted tree
 * @param u - last visited vertex in dfs
 * @param visited - vector of visited vertices
 * @param bh - store the result for branch-height
 * @param bl - store the result for branch-leaf-number
 */
void dfsBranch(
		const RootedTree& branch,
		const Vertex& u,
		vector < bool >& visited,
		int &bh, int &bl){

    if (branch.graph.adj[u].size() > 2){
    	++bh;
    }
    if (branch.graph.adj[u].size() == 1){
        ++bh;
        ++bl;
    }
    int bhmax = bh;
    if (branch.graph.adj[u].size() > 1) {
        for (auto v : branch.graph.adj[u]) {
        	if (visited[v] == false) {
                int bh1 = bh;
                // append vertex v at the end of the visited vector
                visited[v] = true;
                dfsBranch(branch, v, visited, bh1, bl);
                // restore the vector after recursive call
                visited[v] = false;
                if (bh1 > bhmax) {
                	bhmax = bh1;
                } // if (bh1 > bhmax)
        	} // if (itr == visited.end())
        } // for (auto v : branch.graph.adj[u])
    } // if (branch.graph.adj[u].size() > 1)
    bh = bhmax;
}


/**
 * A helper function in calculationg branch-height
 *   and branch-leaf-number
 * Initializes a dfs search starting at the root of a rooted tree
 */
void calcBranchHeightOne(const RootedTree& branch, int &bh, int &bl){
	Vertex u = branch.root;
	Vertex v = branch.graph.adj[u][0];
	// Initialize a boolean vector giving a boolean flag
	// for each vertex whether it has been visited or not
	vector <bool> visited(branch.graph.numAtom+1, false);
	visited[u] = visited[v] = true;
    dfsBranch(branch, v, visited, bh, bl);
}


/**
 * Calculate tyhe k-branch height/number of a graph g
 * @param g - input graph
 * @param k - integer >= 1
 * @param branchheight - store calculated answer
 * @param branchnumber - store calculated answer
 */
void calcBranchHeight(
		const Graph &g,
		const size_t& k,
		size_t &branch_height,
		size_t &branch_leaf_number
		){

    // this function calculates the k-branch-height and k-branch-number of graph g.
    vector <Vertex> coreset = {};
    calcCoreVertexSet(g, coreset);
    Graph gp(g);
    cutTree(gp, k);
    if (debug)
    	cerr << "Calculated cut tree" << endl;

    if (coreset.size() == 1 || coreset.size() > 2) {
    	if (debug)
    		cerr << "coreset.cize == 1" << endl;
        vector <RootedTree> coretree;
        getCoreTree(gp, coretree, coreset);
        if (debug)
        	cerr << "Got core height" << endl;
        int bhmax = 0;
        int bl = 0;
        // if (coreset.size() == 1) ++bn;
        for (auto branch : coretree) {
            int bh = 0;
            calcBranchHeightOne(branch, bh, bl);
            if (bh > bhmax){
                bhmax = bh;
            }
        }
        branch_height = bhmax;
        branch_leaf_number = bl;
    } else if (coreset.size() == 2) { // this situation only happens when graph g is a tree
    	if (debug)
    		cerr << "coreset.size == 2" << endl;
        vector <RootedTree> coretree = {};
        vector <Vertex> coreset1 = {coreset[0]};
        getCoreTree(gp, coretree, coreset1);
        if (debug)
        	cerr << "Got core tree" << endl;
        int bhmax1 = 0;
        int bl1 = 0;
        for (auto branch : coretree) {
            int bh1 = 0;
            calcBranchHeightOne(branch, bh1, bl1);
            if (bh1 > bhmax1){
                bhmax1 = bh1;
            }
        } // for (auto branch : coretree)
        coretree.clear();
        vector <Vertex> coreset2 = {coreset[0]};
        getCoreTree(gp, coretree, coreset2);
        int bhmax2 = 0;
        int bl2 = 0;
        for (auto branch : coretree) {
            int bh2 = 0;
            calcBranchHeightOne(branch, bh2, bl2);
            if (bh2 > bhmax2) {
                bhmax2 = bh2;
            }
        } // for (auto branch : coretree)

        if (bl1 <= bl2) {
            if (bl1 == bl2 && bhmax1 > bhmax2) {
                branch_height = bhmax2;
                branch_leaf_number = bl2;
            } else {
                branch_height = bhmax1;
                branch_leaf_number = bl1;
            } // (bl1 == bl2 && bhmax1 > bhmax2)
        } else {
            branch_height = bhmax2;
            branch_leaf_number = bl2;
        } // if (bl1 <= bl2)
    } // if (coreset.size() == 2)
}


/**
 * Calculate the total feature vector
 */
void getFeatureVec(
		const vector < Graph > &g,
		vector < string > &f_name,
		vector < double > f_value[],
		const vector < size_t > & diameter,
		const vector < map< string, size_t > > &element_map_in,
		const vector < map< string, size_t > > &element_map_ex,
		const vector < map< string, size_t > > &onePath_in,
		const vector < map< string, size_t > > &onePath_ex,
		const vector < map< string, size_t > > &mult_path_map,
		const vector < string > &globalElement,
		const vector < string > &globalOnePath,
		const vector < string > &globalMoreThanOnePath,
		const vector < size_t > &num_atom,
		const vector < double > &M,
		const vector < vector < size_t > > &degree_vertices_in,
		const vector < vector < size_t > > &degree_vertices_ex,
		const vector < size_t > &double_bond_in,
		const vector < size_t > &double_bond_ex,
		const vector < size_t > &triple_bond_in,
		const vector < size_t > &triple_bond_ex,
		const vector < size_t > &num_H_atoms,
		const vector < map < vector < size_t >, size_t > > &bond_conf_in,
		const vector < map < vector < size_t >, size_t > > &bond_conf_ex,
		const vector < size_t > &branch_height,
		const vector < size_t > &branch_leaf_number
		){

	// Add the descriptor names to a vector
	f_name.push_back("CID");
	f_name.push_back("n");
	f_name.push_back("M");

	// Write the title for element descriptors
	for (size_t i = 0; i < globalElement.size(); ++i) {
		f_name.push_back(globalElement[i] + "_in");
		f_name.push_back(globalElement[i] + "_ex");
	}
	f_name.push_back("H"); // add H after all other element symbols
	for (size_t i = 0; i < globalOnePath.size(); ++i) {
		f_name.push_back(globalOnePath[i] + "_in");
		f_name.push_back(globalOnePath[i] + "_ex");
	}
	for (size_t i = 0; i < globalMoreThanOnePath.size(); ++i) {
		f_name.push_back(globalMoreThanOnePath[i]);
	}

	f_name.push_back("#degree1_in");
	f_name.push_back("#degree1_ex");
	f_name.push_back("#degree2_in");
	f_name.push_back("#degree2_ex");
	f_name.push_back("#degree3_in");
	f_name.push_back("#degree3_ex");
	f_name.push_back("#degree4_in");
	f_name.push_back("#degree4_ex");
	f_name.push_back("#double_bond_in");
	f_name.push_back("#double_bond_ex");
	f_name.push_back("#triple_bond_in");
	f_name.push_back("#triple_bond_ex");
	f_name.push_back("Diameter");
	f_name.push_back("Bc_121_in");
	f_name.push_back("Bc_121_ex");
	f_name.push_back("Bc_122_in");
	f_name.push_back("Bc_122_ex");
	f_name.push_back("Bc_123_in");
	f_name.push_back("Bc_123_ex");
	f_name.push_back("Bc_131_in");
	f_name.push_back("Bc_131_ex");
	f_name.push_back("Bc_132_in");
	f_name.push_back("Bc_132_ex");
	f_name.push_back("Bc_141_in");
	f_name.push_back("Bc_141_ex");
	f_name.push_back("Bc_221_in");
	f_name.push_back("Bc_221_ex");
	f_name.push_back("Bc_222_in");
	f_name.push_back("Bc_222_ex");
	f_name.push_back("Bc_223_in");
	f_name.push_back("Bc_223_ex");
	f_name.push_back("Bc_231_in");
	f_name.push_back("Bc_231_ex");
	f_name.push_back("Bc_232_in");
	f_name.push_back("Bc_232_ex");
	f_name.push_back("Bc_241_in");
	f_name.push_back("Bc_241_ex");
	f_name.push_back("Bc_331_in");
	f_name.push_back("Bc_331_ex");
	f_name.push_back("Bc_332_in");
	f_name.push_back("Bc_332_ex");
	f_name.push_back("Bc_341_in");
	f_name.push_back("Bc_341_ex");
	f_name.push_back("Bc_441_in");
	f_name.push_back("Bc_441_ex");

	f_name.push_back("2-branch_height");
	f_name.push_back("2-branch_leaf_number");

	for (size_t i = 0; i < g.size(); ++i) {
		f_value[i].push_back(stoi(g[i].CID));
		f_value[i].push_back(num_atom[i]);
		f_value[i].push_back((double)M[i]/num_atom[i]);
		for (size_t j = 0; j < globalElement.size(); ++j) {
			auto loc = element_map_in[i].find(globalElement[j]);
			if (loc == element_map_in[i].end()) {
				f_value[i].push_back(0);
			} else {
				f_value[i].push_back(
						element_map_in[i].at(globalElement[j]));
			}
			loc = element_map_ex[i].find(globalElement[j]);
			if (loc == element_map_ex[i].end()) {
				f_value[i].push_back(0);
			} else {
				f_value[i].push_back(
						element_map_ex[i].at(globalElement[j]));
			}
		}
		// write number of H atoms value after writing all elements
		f_value[i].push_back(num_H_atoms[i]);
		for (size_t j = 0; j < globalOnePath.size(); ++j) {
			auto loc = onePath_in[i].find(globalOnePath[j]);
			if (loc == onePath_in[i].end()) {
				f_value[i].push_back(0);
			} else {
				f_value[i].push_back(loc->second);
			}
			loc = onePath_ex[i].find(globalOnePath[j]);
			if (loc == onePath_ex[i].end()) {
				f_value[i].push_back(0);
			} else {
				f_value[i].push_back(loc->second);
			}
		}
		for (size_t j = 0; j < globalMoreThanOnePath.size(); ++j) {
			auto loc = mult_path_map[i].find(globalMoreThanOnePath[j]);
			if (loc == mult_path_map[i].end()) {
				f_value[i].push_back(0);
			} else {
				f_value[i].push_back(loc->second);
			}
		}
		for (size_t j = 0; j < 4; ++j) {
			f_value[i].push_back(degree_vertices_in[i][j]);
			f_value[i].push_back(degree_vertices_ex[i][j]);
		}
		f_value[i].push_back(double_bond_in[i]);
		f_value[i].push_back(double_bond_ex[i]);
		f_value[i].push_back(triple_bond_in[i]);
		f_value[i].push_back(triple_bond_ex[i]);
		f_value[i].push_back((double)diameter[i] / num_atom[i]);
		auto deg_t_in = bond_conf_in[i];
		auto deg_t_ex = bond_conf_ex[i];
		for (auto bc_t : set_Bc) {
			f_value[i].push_back(deg_t_in[bc_t]);
			f_value[i].push_back(deg_t_ex[bc_t]);
		}
		// for (const auto &nc : globalNc)
		// 	f_value[i].push_back(ncMaps[i].at(nc));
		f_value[i].push_back(branch_height[i]);
		f_value[i].push_back(branch_leaf_number[i]);
	}
}


/**
 * Calculate sum of distances between each pair of two vertices
 * @param Graph g: input graph
 * @param sumDist: stores the answer value of sum of distances
 * @param diameter: stores the answer value of diamter
 */
void
calcDist(const Graph &g, int &sumDist, int &diameter) {

	// Create a 2D array to store the distances between all pairs of vertices
	vector < vector <int> > all_dist;

	for (size_t i = 0; i < g.numAtom; ++i) {
		vector <int> dist(g.numAtom, 0);
		vector <bool> visited(g.numAtom, false);

		queue <Vertex> Q;
		Q.push(Vertex(i));

		while (!Q.empty()) {
			Vertex u = Q.front();
			for (Vertex v : g.adj[u]) {
				if (!visited[v]) {
					dist[v] = dist[u] + 1;
					visited[v] = 1;
					Q.push(v);
				}
			}
			Q.pop();
		}
		all_dist.push_back(dist);
	}

	sumDist = 0;
	diameter = 0;
	// Calculate the sum of distances and diameter
	for (size_t i = 0; i < all_dist.size(); ++i) {
		for (size_t j = i+1; j < all_dist[0].size(); ++j) {
			sumDist += all_dist[i][j];
			diameter = max(diameter, all_dist[i][j]);
		}
	}
}


/**
 * Print a table of feature vector values to a file
 * @param num_instance - number of chemical graphs in the input sdf
 * @param f_name - names of the descriptors
 * @param f_value - values of the descriptors for each chemical graph
 * @param out_file_name - filename to write the output
 */
void printFeatureVec(
		size_t &num_instance,
		vector < string > &f_name,
		vector < double > f_value[],
		string out_file_name
	){

	// write out feature vector to the supplied file
	ofstream outputfile(out_file_name, ios::out);


	for (size_t i = 0; i < f_name.size()-1; ++i) {
		outputfile << f_name[i] << ",";
	}
	outputfile << f_name[f_name.size()-1] << "\n";

	for (size_t i = 0; i < num_instance; ++i) {
		for (size_t j = 0; j < f_value[i].size()-1; ++j) {
			if (j == 0) outputfile << (int) f_value[i][j] << ",";
			else outputfile << f_value[i][j] << ",";
		}
		outputfile << f_value[i][f_value[i].size()-1] << "\n";
	}

	outputfile.close();
}


int main(int argc, char **argv){
	// Check if there are input parameters supplied
	if (argc != 3) {
		cout << "Please supply as command line arguments " <<
				"an input file name (sdf) and an output file name." << endl;
		cout << "The program will now exit." << endl;
		return(1);
	}
	string input_filename = argv[1];
	string output_filename = argv[2];

	cerr << "Working on input file\n  " <<
			input_filename << "..." << endl;

	// Initialize a map with atomic mass for different elements
	MassMap massMap = init_MassMap();

	// Initialize a map for storing valence of different atoms
	ValenceMap valMap = init_ValenceMap();

	// Storing the chemical graphs in the input sdf
	vector < Graph > H;
	try {
		H = read_graph_sdf(input_filename);
	} catch (const exception &e) {
		cout << "Error in reading input file.\n" <<
				e.what() << endl <<
				"The program will now terminate." << endl;
		return(1);
	}

	size_t num_graphs = H.size();

	// Prepare storage for H-suppressed graphs
	vector < Graph > G(num_graphs);

	// Prepare storage for each descriptor
	vector < size_t > diameter(num_graphs);		 // Graph diameter
	vector < size_t > sumDist(num_graphs);		 // Stores sum of distances
	vector < size_t > num_atom(num_graphs); 	 // Number of non-hydrogen atoms
	vector < size_t > double_bond_in(num_graphs);	 // Number of internal double bonds
	vector < size_t > double_bond_ex(num_graphs);	 // Number of external double bonds
	vector < size_t > triple_bond_in(num_graphs); 	 // Number of internal triple bonds
	vector < size_t > triple_bond_ex(num_graphs); 	 // Number of external triple bonds
	vector < double > M(num_graphs);		 	 // Sum of atomic mass
	vector < size_t > num_H_atoms(num_graphs); 	 // Number of H-atoms
	vector < size_t > branch_height(num_graphs); // Branch height
	vector < size_t > branch_leaf_number(num_graphs);   // Branch-leaf number
	// Introduce integer IDs in mappings
	vector < map < string, size_t > > element_map_in(num_graphs);
	vector < map < string, size_t > > element_map_ex(num_graphs);
	vector < map < string, size_t > > onepath_map_in(num_graphs);
	vector < map < string, size_t > > onepath_map_ex(num_graphs);
	vector < map < string, size_t > > mult_path_map(num_graphs);

	vector < string > globalElement;
	globalElement.push_back("C");	// Assign carbon to the first order
	vector < string > globalOnePath;
	vector < string > globalMoreThanOnePath;
	vector < vector < size_t > > degree_vertices_in(num_graphs);
	vector < vector < size_t > > degree_vertices_ex(num_graphs);

	for (size_t i = 0; i < num_graphs; ++i) {
		num_H_atoms[i] = H_suppressed_convert(H[i], G[i]);

		Graph _g = G[i];
		cutTree(_g, 2);

		try {
			calcDegree(
					G[i],
					degree_vertices_in[i],
					degree_vertices_ex[i],
					_g.status );
		} catch (const exception& e) {
			cout << "Exception caught in calculating degree:" << endl;
			cout << e.what() << endl;
			cout << "The program will now exit" << endl;
			return (1);
		}
		num_atom[i] = G[i].numAtom;
		calcDoubleTripleBond(G[i], double_bond_in[i], double_bond_ex[i],
						triple_bond_in[i], triple_bond_ex[i], _g.status);
		try {
			M[i] = calcM(G[i], massMap);
		} catch (const exception &ex) {
			cerr << "An exception when calculating the mass "
					<< "of graph with id " << G[i].CID << endl;
			cerr << "Perhaps it contains an element not registered "
					<< "in the table of atomic masses?" << endl;
			cerr << ex.what() << endl;
			cerr << "The program will now terminate" << endl;
			return 1;
		}
		getElementCount(G[i],
						element_map_in[i],
						element_map_ex[i],
						globalElement,
						_g.status);
		calcBranchHeight(G[i], 2, branch_height[i], branch_leaf_number[i]);
	}

	map<string, size_t> elementOrder;
	assignElemOrder(globalElement, elementOrder);

	vector <map <vector <size_t>, size_t>> bond_conf_in;
	vector <map <vector <size_t>, size_t>> bond_conf_ex;

	for (size_t i = 0; i < num_graphs; ++i) {
		Graph _g = G[i];
		cutTree(_g, 2);
		getOnePath( G[i], elementOrder, onepath_map_in[i], onepath_map_ex[i],
					globalOnePath, _g.status );
		int sumD, dia;
		calcDist(G[i], sumD, dia);
		sumDist[i] = sumD;
		diameter[i] = dia;
		map <vector <size_t>, size_t> bond_conf_t_in;
		bond_conf_in.push_back(bond_conf_t_in);
		map <vector <size_t>, size_t> bond_conf_t_ex;
		bond_conf_ex.push_back(bond_conf_t_ex);
		calcBondConfiguration(G[i], bond_conf_in[i],
							  bond_conf_ex[i], _g.status);
	}

	// A vector that stores the descriptor names
	vector < string > f_name;
	// A 2D structure to store the descriptor values for each input graph
	vector < double > f_value[num_graphs];

	// Store the feature vector values in a 2D structure
	getFeatureVec( G, f_name, f_value,
				   diameter, element_map_in, element_map_ex,
				   onepath_map_in, onepath_map_ex, mult_path_map,
				   globalElement, globalOnePath, globalMoreThanOnePath,
				   num_atom, M, degree_vertices_in, degree_vertices_ex,
				   double_bond_in, double_bond_ex, triple_bond_in, triple_bond_ex,
				   num_H_atoms, bond_conf_in, bond_conf_ex,
				   branch_height, branch_leaf_number
				 );

	cerr << "Finished calculating FV4+Bc+(in/ex)" << endl;

	// Print the stored values to a given file
	printFeatureVec(num_graphs, f_name, f_value, output_filename);

	// Report that the program has finished and exit
	cerr << "The feature vectors have been written to the file\n  " <<
			output_filename << endl;
	return 0;
}
