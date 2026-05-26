#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <vector>

#include "commonused.hpp"

using namespace std;

const int maxAtomNum = 1000;
const int maxDegree = 6;

typedef size_t Vertex;
class ChemicalGraph;

void searchBranch(Vertex u, vector <Vertex>& forbiddenset, ChemicalGraph& graph, ChemicalGraph& branch);

class ChemicalGraph{

public:
    string CID;
    int numAtom;
    vector <vector <Vertex>> adj;
    vector <string> alpha;
    vector <vector <int>> beta;
    vector <bool> status;
	vector <int> charge;

    ChemicalGraph(){
        CID = "";
        numAtom = 0;
        adj = {};
        alpha = {};
        beta = {};
        status = {};
		charge = {};
    }

    ChemicalGraph(const ChemicalGraph& graph){
        CID = graph.CID;
        numAtom = graph.numAtom;
        adj = graph.adj;
        alpha = graph.alpha;
        beta = graph.beta;
        status = graph.status;
		charge = graph.charge;
    }

    ChemicalGraph(const ChemicalGraph& graph, vector <Vertex>& efflist){
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

    int EffectiveAtomNum(){
        int EAN = 0;
        for (Vertex u = 0; u < numAtom; u++){
            if (status[u]) {
                EAN++;
            }
        }
        return EAN;
    }

    int degree(Vertex u){
        if (status[u]) {
            return adj[u].size();
        } else {
            return 0;
        }
    }

	// for 2LM-M //
	int getVal(Vertex u) {
		int beta_u = 0;
		for (Vertex v : adj[u]) {
			beta_u = beta_u + beta[u][v];
		}
        beta_u -= charge[u]; // 0503 v013, remove the effect of charge
		return beta_u;
	}

	// for 2LM-M //
	void atom_name() {
		for (Vertex u = 0; u < numAtom; u++) { // traversing over all vertices
			alpha[u] += to_string(getVal(u));
		}
	}

    int maxd(){
        int ans = 0;
        for (Vertex u = 0; u < numAtom; u++){
            int dtemp = degree(u);
            if (dtemp > ans){
                ans = dtemp;
            }
        }
        return ans;
    }

    void HSuppress(){
        for (Vertex u = 0; u < numAtom; u++){
            if (alpha[u] == "H" || alpha[u] == "H1"){
                RemoveVertex_Simple(u);
            }
        }
    }

    void CutLeaf_OneLayer(){
        auto _adj = adj;
        if (EffectiveAtomNum() <= 2){
            return;
        }
        for (Vertex u = 0; u < numAtom; u++){
            if (_adj[u].size() == 1){
                RemoveVertex_Simple(u);
            }
        }
    }

    void CutLeaf(int k){
        for (int i = 0; i < k; i++){
            CutLeaf_OneLayer();
        }
    }

    void CutUntilCore(){
        bool loop = true;
        while (loop){
            int numv = EffectiveAtomNum();
            if (numv > 2){
                CutLeaf_OneLayer();
                int numv2 = EffectiveAtomNum();
                if (numv == numv2){
                    loop = false;
                }
            } else {
                loop = false;
            }
        }
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

    vector <Vertex> calcEffectiveAtom(){
        vector <Vertex> efflist = {};
        for (Vertex u = 0; u < numAtom; u++){
            if (status[u]){
                efflist.push_back(u);
            }
        }
        return efflist;
    }

    int readFromFile(ifstream& infile, bool noion = true, bool noradical = true){
        bool flag = true;
        string line;
        int ans = 1;
        while (true){
            getline(infile, line);
            if (infile.eof()){
                return 0;
            }
            if (flag) {
                CID = line;
                getline(infile, line);
                getline(infile, line);
                getline(infile, line);
                if (infile.eof()){
                    return 0;
                }
                int n, m;
                stringstream st;
                st << line;
                string strtemp = st.str();
                string stt = strtemp.substr(0, 3);
                stt.erase(remove(stt.begin(), stt.end(), ' '), stt.end());
                n = stoi(stt);
                stt = strtemp.substr(3, 3);
                stt.erase(remove(stt.begin(), stt.end(), ' '), stt.end());
                m = stoi(stt);
                initialize(n);
                double x, y, z;
                string al;
                for (Vertex u = 0; u < n; u++){
                    getline(infile, line);
                    status[u] = true;
                    stringstream st1;
                    st1 << line;
                    st1 >> x;
                    st1 >> y;
                    st1 >> z;
                    st1 >> al;
                    alpha[u] = al;
                }
                int v1, v2, mul;
                for (int i = 0; i < m; i++){
                    getline(infile, line);
                    stringstream st2;
                    st2 << line;
                    strtemp = st2.str();
                    stt = strtemp.substr(0, 3);
                    stt.erase(remove(stt.begin(), stt.end(), ' '), stt.end());
                    v1 = stoi(stt);
                    stt = strtemp.substr(3, 3);
                    stt.erase(remove(stt.begin(), stt.end(), ' '), stt.end());
                    v2 = stoi(stt);
                    stt = strtemp.substr(6, 3);
                    stt.erase(remove(stt.begin(), stt.end(), ' '), stt.end());
                    mul = stoi(stt);
        
                    adj[v1 - 1].push_back(v2 - 1);
                    adj[v2 - 1].push_back(v1 - 1);

                    beta[v1 - 1][v2 - 1] = mul;
                    beta[v2 - 1][v1 - 1] = mul;
                }
                flag = false;
            } else {            
                if (line.substr(0, 4) == "$$$$"){
                    break;
                }
				if (line.substr(0, 6) == "M  RAG") {
					if (noion) {
						ans = -2;
					}
				}
                if (line.substr(0, 6) == "M  CHG"){
                    if (noion){
                        ans = -1;
                    }
					else {// read the charge and store it
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
							charge[ind_atom - 1] = atom_chg;
						}
					
					}
                }
            }
        }
        return ans;
    }

private:
    void RemoveVertex_Simple(Vertex u){
        for (auto v : adj[u]){
            adj[u].erase(remove(adj[u].begin(), adj[u].end(), v), adj[u].end());
            adj[v].erase(remove(adj[v].begin(), adj[v].end(), u), adj[v].end());
        }
        status[u] = false;
    }
};

void searchBranch(Vertex u, vector <Vertex>& forbiddenset, ChemicalGraph& graph, ChemicalGraph& branch){
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
    ChemicalGraph graph;
    vector <int> depth;

    RootedTree(Vertex _root = 0, ChemicalGraph _graph = {}){
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

    void buildFromCoreVertex(Vertex _root, vector <Vertex> CoreSet, ChemicalGraph& _graph){
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

ChemicalGraph RemoveNullVertex(ChemicalGraph& graph){
    int EAN = graph.EffectiveAtomNum();
    vector <Vertex> efflist = graph.calcEffectiveAtom();
    ChemicalGraph new_graph(graph, efflist); 
    return new_graph;
}

RootedTree RemoveNullVertex(RootedTree& RT){
    int EAN = RT.graph.EffectiveAtomNum();
    vector <Vertex> efflist = RT.graph.calcEffectiveAtom();
    ChemicalGraph new_graph(RT.graph, efflist);
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

vector <Vertex> getBranchTree(ChemicalGraph& graph, int k){
    ChemicalGraph _graph = graph;
    _graph.CutLeaf(k);
    vector <Vertex> BranchTree = {};
    for (Vertex u = 0; u < _graph.numAtom; u++){
        if (_graph.status[u] == true){
            BranchTree.push_back(u);
        }
    }
    return BranchTree;
}        

vector <RootedTree> getFringeTrees(ChemicalGraph& graph, int k){
    ChemicalGraph _graph = graph;
    vector <Vertex> BranchTree = getBranchTree(_graph, k);
    vector <RootedTree> FringeTrees = {};
    for (auto u : BranchTree){
        RootedTree ct;
        ct.buildFromCoreVertex(u, BranchTree, graph);
        if (ct.graph.EffectiveAtomNum() > 1){
            FringeTrees.push_back(ct);
        }
    }
    return FringeTrees;
}

vector <RootedTree> getFringeTrees_H(ChemicalGraph& graph, ChemicalGraph& graph_no_H, int k){
    // ChemicalGraph _graph = graph;
    vector <Vertex> BranchTree = getBranchTree(graph_no_H, k);
    vector <RootedTree> FringeTrees = {};
    for (auto u : BranchTree){
        RootedTree ct;
        ct.buildFromCoreVertex(u, BranchTree, graph);
        // if (ct.graph.EffectiveAtomNum() > 1){
            FringeTrees.push_back(ct);
        // }
    }
    return FringeTrees;
}
