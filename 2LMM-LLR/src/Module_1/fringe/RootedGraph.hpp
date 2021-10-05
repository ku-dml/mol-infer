#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>

#include "commonused.hpp"
#include "TopologyGraph.hpp"

typedef size_t Color;

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

class TreeSignature{
public:
    vector <int> delta;
    vector <int> mu;
	vector <int> chg;

    TreeSignature(){}

    TreeSignature(const vector <int>& _delta, const vector <int>& _mu, 
		const vector <int>& _chg){
        delta = _delta;
        mu = _mu;
		chg = _chg;
    }

    bool is_equal(TreeSignature& TS){
        if (lexicographical_comparison(delta, TS.delta) == 0 && 
			lexicographical_comparison(mu, TS.mu) == 0 &&
			lexicographical_comparison(chg, TS.chg) == 0){
            return true;
        } else {
            return false;
        }
    }

    void extend(const TreeSignature& _TS){
        delta.insert(delta.end(), _TS.delta.begin(), _TS.delta.end());
        mu.insert(mu.end(), _TS.mu.begin(), _TS.mu.end());
		chg.insert(chg.end(), _TS.chg.begin(), _TS.chg.end());
    }

    void print(){
        cout << "C(K) = ((";
        for (int i = 0; i < delta.size(); i++){
            cout << delta[i] << " ";
        }
        cout << "),(";
        for (int i = 0; i < mu.size(); i++){
            cout << mu[i] << " ";
        }
		cout << "),(";
		for (int i = 0; i < chg.size(); i++) {
			cout << chg[i] << " ";
		}
        cout << "))" << endl;

    }
};

bool TS_cmp(const TreeSignature& a, const TreeSignature& b){
    return ((lexicographical_comparison(a.delta, b.delta) > 0) ||
            (lexicographical_comparison(a.delta, b.delta) == 0 &&
                lexicographical_comparison(a.mu, b.mu) > 0)  ||
			(lexicographical_comparison(a.delta, b.delta) == 0 &&
			lexicographical_comparison(a.mu, b.mu) == 0 &&
			lexicographical_comparison(a.chg, b.chg) > 0));
}

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

namespace std{
    template <>
    class hash <TreeSignature>{
    public:
        size_t operator()(const TreeSignature& aa) const {
            size_t res = 0;
            hash_combine <vector <int>>(res, aa.delta);
            hash_combine <vector <int>>(res, aa.mu);
			hash_combine <vector <int>>(res, aa.chg);
            return res;
        }
    };
}

bool operator==(const TreeSignature& a, const TreeSignature& b){
    return (a.delta == b.delta && a.mu == b.mu && a.chg == b.chg);
}
