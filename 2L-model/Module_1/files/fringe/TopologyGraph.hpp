#pragma once

#include <string>
#include <vector>
#include <algorithm>

#include "ChemicalGraph.hpp" 
#include "commonused.hpp"

class TopologyGraph{
    
public:
    string name;
    int numVertex;

    vector <vector <Vertex>> adj;
    vector <bool> status;

    TopologyGraph(){
        name = "";
        numVertex = 0;
        adj = {};
        status = {};
    }

    TopologyGraph(const TopologyGraph& Tgraph){
        name = Tgraph.name;
        numVertex = Tgraph.numVertex;
        adj = Tgraph.adj;
        status = Tgraph.status;
    }

    TopologyGraph(const ChemicalGraph& graph){
        name = graph.CID;
        numVertex = graph.numAtom;
        adj = graph.adj;
        status = graph.status;
    }

    TopologyGraph(const TopologyGraph& Tgraph, vector <Vertex>& efflist){
        numVertex = efflist.size();
        initialize(numVertex);
        name = Tgraph.name;
        for (Vertex u = 0; u < numVertex; u++){
            status[u] = true;
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
        vector <Vertex> tmp;
        for (Vertex u = 0; u < num; u++){
            adj.push_back(tmp);
            status[u] = false;
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
