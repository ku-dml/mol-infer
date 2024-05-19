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

// To prepare
// valence of chemical element
// new, 0430
unordered_map<string, int> valence_map = {
    {"C", 4},  {"O", 2}, {"N", 3}, {"S", 2},
    {"Cl", 1}, {"F", 1}, {"P", 5}, {"Si", 4},  {"H", 1}};
unordered_map<string, size_t> mass_map {
    {"H", 10},    {"He", 40},   {"Li", 70},   {"Be", 90},   {"B", 108},
    {"C", 120},   {"N", 140},   {"O", 160},   {"F", 190},   {"Ne", 200},
    {"Na", 230},  {"Mg", 240},  {"Al", 270},  {"Si", 280},  {"P", 310},
    {"S", 320},   {"Cl", 355},  {"Ar", 400},  {"K", 390},   {"Ca", 400},
    {"Sc", 450},  {"Ti", 479},  {"V", 509},   {"Cr", 520},  {"Mn", 549},
    {"Fe", 558},  {"Co", 589},  {"Ni", 587},  {"Cu", 635},  {"Zn", 654},
    {"Ga", 697},  {"Ge", 726},  {"As", 749},  {"Se", 790},  {"Br", 800},
    {"Kr", 838},  {"Rb", 855},  {"Sr", 876},  {"Y", 889},   {"Zr", 912},
    {"Nb", 929},  {"Mo", 956},  {"Tc", 989},  {"Ru", 1010}, {"Rh", 1029},
    {"Pd", 1064}, {"Ag", 1079}, {"Cd", 1124}, {"In", 1148}, {"Sn", 1187},
    {"Sb", 1218}, {"Te", 1276}, {"I", 1270},  {"Xe", 1313}, {"Cs", 1329},
    {"Ba", 1273}, {"La", 1389}, {"Ce", 1401}, {"Pr", 1409}, {"Nd", 1442},
    {"Pm", 1469}, {"Sm", 1504}, {"Eu", 1520}, {"Gd", 1573}, {"Tb", 1589},
    {"Dy", 1625}, {"Ho", 1649}, {"Er", 1673}, {"Tm", 1689}, {"Yb", 1731},
    {"Lu", 1750}, {"Hf", 1785}, {"Ta", 1809}, {"W", 1838},  {"Re", 1862},
    {"Os", 1902}, {"Ir", 1922}, {"Pt", 1951}, {"Au", 1970}, {"Hg", 2006},
    {"Tl", 2044}, {"Pb", 2072}, {"Bi", 2090}, {"Po", 2090}, {"At", 2100},
    {"Rn", 2220}, {"Fr", 2230}, {"Ra", 2260}, {"Ac", 2270}, {"Th", 2320},
    {"Pa", 2310}, {"U", 2380},  {"Np", 2370}, {"Pu", 2441}, {"Am", 2431},
    {"Cm", 2471}, {"Bk", 2471}, {"Cf", 2511}, {"Es", 2521}, {"Fm", 2571},
    {"Md", 2581}, {"No", 2591}, {"Lr", 2601}, {"Rf", 2611}, {"Db", 2621},
    {"Sg", 2631}, {"Bh", 2621}, {"Hs", 2770}, {"Mt", 2780}, {"Ds", 2810},
    {"Rg", 2840}, {"Cn", 2880}, {"Nh", 2930}, {"Fl", 2980}, {"Mc", 2990},
    {"Lv", 3020}, {"Ts", 3100}, {"Og", 3140}};
unordered_map <string, int> map_atomic_number{{"e*", 9999}, {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, 
                                              {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10},
                                              {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15},
                                              {"S", 16}, {"Cl", 17}, {"Ar", 18}, {"K", 19}, {"Ca", 20},
                                              {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25},
                                              {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30},
                                              {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35},
                                              {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40},
                                              {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44}, {"Rh", 45},
                                              {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50},
                                              {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54}, {"Cs", 55},
                                              {"Ba", 56}, {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60},
                                              {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64}, {"Tb", 65},
                                              {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70},
                                              {"Lu", 71}, {"Hf", 72}, {"Ta", 73}, {"W", 74}, {"Re", 75},
                                              {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80},
                                              {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85},
                                              {"Rn", 86}, {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90},
                                              {"Pa", 91}, {"U", 92}, {"Np", 93}, {"Pu", 94}, {"Am", 95},
                                              {"Cm", 96}, {"Bk", 97}, {"Cf", 98}, {"Es", 99}, {"Fm", 100},
                                              {"Md", 101}, {"No", 102}, {"Lr", 103}, {"Rf", 104}, {"Db", 105},
                                              {"Sg", 106}, {"Bh", 107}, {"Hs", 108}, {"Mt", 109}, {"Ds", 110},
                                              {"Rg", 111}, {"Cn", 112}, {"Nh", 113}, {"Fl", 114}, {"Mc", 115},
                                              {"Lv", 116}, {"Ts", 117}, {"Og", 118}, {"*", 999}
                                            };

vector <string> split(const string& input, char delimiter){
    istringstream stream(input);
    string field;
    vector <string> result;
    while (getline(stream, field, delimiter)){
        result.push_back(field);
    }
    return result;
}

//////////
// add for computing TS, 0609
template <class T>
int lexicographical_comparison(const vector <T>& a, const vector <T>& b){
    int min_len = min(a.size(), b.size());
    for (int i = 0; i < min_len; i++){
        if (a[i] > b[i]){
            return -1;
        }
        if (a[i] < b[i]){
            return 1;
        }
    }
    if (a.size() > b.size())
        return -1;
    if (a.size() < b.size())
        return 1;
    return 0;
    // -1 a > b
    // 0 a = b
    // 1 a < b
}

int lexicographical_comparison(vector <vector <int>> a, vector <vector <int>> b){
    int min_len = min(a.size(), b.size());
    for (int i = 0; i < min_len; i++){
        int temp = lexicographical_comparison(a[i], b[i]);
        if (temp == -1){
            return -1;
        }
        if (temp == 1){
            return 1;
        }
    }
    if (a.size() > b.size())
        return -1;
    if (a.size() < b.size())
        return 1;
    return 0;
}

bool lexicographical_equal(int n, int m, int a[], int b[]){
    if (n != m){
        return false;
    }

    for (int i = 0; i < n; i++){
        if (a[i] != b[i]){
            return false;
        }
    }
    return true;
}

template <class T>
bool operator==(const vector <T>& a, const vector <T>& b){
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); i++){
        if (a[i] != b[i]){
            return false;
        }
    }
    return true;
}


class TreeSignature{
public:
    vector <string> delta_2;  // delta is (atom, dep)-seq
    vector <int> delta;
    vector <int> mu;
    // new, 0430
    vector <int> chg; 

    TreeSignature(){}

    TreeSignature(const vector <string>& _delta_2, const vector <int>& _mu, const vector <int>& _chg){
        delta_2 = _delta_2;

        delta.clear();
        for (auto& d2 : delta_2){
            delta.push_back(map_atomic_number.at(d2));
        }

        mu = _mu;
        // new, 0430
        chg = _chg;
    }

    bool is_equal(const TreeSignature& TS){
        if (lexicographical_comparison(delta, TS.delta) == 0 && 
            lexicographical_comparison(mu, TS.mu) == 0 &&
            lexicographical_comparison(chg, TS.chg) == 0){
            return true;
        } else {
            return false;
        }
    }

    void extend(const TreeSignature& _TS){
        delta_2.insert(delta_2.end(), _TS.delta_2.begin(), _TS.delta_2.end());
        delta.insert(delta.end(), _TS.delta.begin(), _TS.delta.end());
        mu.insert(mu.end(), _TS.mu.begin(), _TS.mu.end());
        // new, 0430
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

void get_root_H(const TreeSignature& TS, int &num_H) {
	for (size_t i = 0; i < TS.mu.size(); ++i) {
		if ((TS.delta_2[2 * i + 2] == "H1") && 
			(TS.delta_2[2 * i + 3] == "1")) {// check if the current atom is H
			++num_H;
		}
	}
}

void read_TS(const string& csv_filename, map <size_t, TreeSignature>& TS_map){
    TS_map.clear();
    // fc_map.clear();

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

        if (lines.size() < 2){
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
            TS.delta_2.push_back(str_tmp);
            TS.delta_2.push_back(to_string(tmp2));

            if (valence_map.find(str_tmp) == valence_map.end()){
              stringstream st_tmp;
              st_tmp << str_tmp[str_tmp.size() - 1];
              int val_tmp;
              st_tmp >> val_tmp;
              string atom_symbol = "";
              size_t j = 0;
              while ((str_tmp[j] < '0') || (str_tmp[j] > '9')){
                atom_symbol += str_tmp[j];
                ++j;
              }
              valence_map.emplace(str_tmp, val_tmp);
              mass_map.emplace(str_tmp, mass_map.at(atom_symbol));
              map_atomic_number.emplace(str_tmp, val_tmp * 200 + map_atomic_number.at(atom_symbol));
            }

            TS.delta.push_back(map_atomic_number.at(str_tmp));
            TS.delta.push_back(tmp2);
        }

        stringstream st3;
        //st3 << lines[2];

        if (lines.size() >= 3){
            st3 << lines[2];

            int tmp3;
            while (st3 >> tmp3){
                TS.mu.push_back(tmp3);
            }
        }

        // new, 0430
        stringstream st4;

        if (lines.size() >= 4){
            st4 << lines[3];

            int tmp4;
            while (st4 >> tmp4){
                TS.chg.push_back(tmp4);
            }
        }

        TS_map.emplace(id, TS);

        // fc_map.emplace(id. TS_map.size());
    }
}

#endif
