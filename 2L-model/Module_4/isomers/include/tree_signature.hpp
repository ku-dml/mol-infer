#ifndef __TREE_SIGNATURE_HPP__INCLUDED
#define __TREE_SIGNATURE_HPP__INCLUDED

#include <vector>
#include <string>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>

using namespace std;

// #include "tools.hpp"

vector <string> split(const string& input, char delimiter){
    istringstream stream(input);
    string field;
    vector <string> result;
    while (getline(stream, field, delimiter)){
        result.push_back(field);
    }
    return result;
}

class TreeSignature{
public:
    vector <string> delta;  // delta is (atom, dep)-seq
    vector <int> mu;

    TreeSignature(){}

    TreeSignature(const vector <string>& _delta, const vector <int>& _mu){
        delta = _delta;
        mu = _mu;
    }

    // bool is_equal(const TreeSignature& TS){
    //     if (lexicographical_comparison(delta, TS.delta) == 0 && lexicographical_comparison(mu, TS.mu) == 0){
    //         return true;
    //     } else {
    //         return false;
    //     }
    // }

    void extend(const TreeSignature& _TS){
        delta.insert(delta.end(), _TS.delta.begin(), _TS.delta.end());
        mu.insert(mu.end(), _TS.mu.begin(), _TS.mu.end());
    }
};

void read_TS(const string& csv_filename, map <size_t, TreeSignature>& TS_map){
    TS_map.clear();

    ifstream infile;
    string line;

    try {
        infile.open(csv_filename);
    } catch (const exception& e){
        return;
    }

    if (!infile.is_open()){
        return;
    }

    if (infile.eof()){
        return;
    }

    // getline(infile, line);

    while (getline(infile, line)){
        auto lines = split(line, ',');

        if (lines.size() < 3){
            break;
        }

        TreeSignature TS;

        size_t id = stoi(lines[0]);

        TS.delta.clear();
        TS.mu.clear();

        stringstream st2;
        st2 << lines[1];

        string str_tmp;
        int tmp2;
        while (st2 >> str_tmp >> tmp2){
            TS.delta.push_back(str_tmp);
            TS.delta.push_back(to_string(tmp2));
        }

        stringstream st3;
        st3 << lines[2];

        int tmp3;
        while (st3 >> tmp3){
            TS.mu.push_back(tmp3);
        }

        TS_map.emplace(id, TS);
    }
}

#endif