// This file mainly contains the definition and functions for the data structure used in the code.

#ifndef __DATA_STRUCTURE_HPP__INCLUDED
#define __DATA_STRUCTURE_HPP__INCLUDED

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_map>

using namespace std;

const vector <vector <size_t>> Dcset = {{1, 2, 1}, {1, 2, 2}, {1, 2, 3}, {1, 3, 1}, {1, 3, 2}, {1, 4, 1}, {2, 2, 1}, {2, 2, 2},
                                        {2, 2, 3}, {2, 3, 1}, {2, 3, 2}, {2, 4, 1}, {3, 3, 1}, {3, 3, 2}, {3, 4, 1}, {4, 4, 1}};

const vector <vector <unsigned short>> _prt = {{}, {}, {}, {0, 0, 1, 1, 0, 4, 4, 0, 7, 7}, {0, 0, 1, 1, 1, 0, 5, 5, 5, 0, 9, 9, 9}};

template <class T>
using vector_2D = vector <vector <T>>;

template <class T>
using vector_3D = vector <vector_2D <T>>;

template <class T>
using vector_4D = vector <vector_3D <T>>;

template <class T>
inline void hash_combine(std::size_t& seed, const T& v){
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2); 
}

struct Atom{
    string NAME; // the chemical name of an atom
    size_t LABEL;   // the index of an atom  
    size_t VALENCE; // valence           
    size_t MASS;
};

bool operator<(const Atom& lhs, const Atom& rhs){
    return (lhs.MASS < rhs.MASS);
}

bool operator>(const Atom& lhs, const Atom& rhs){
    return (lhs.MASS > rhs.MASS);
}

bool operator==(const Atom& lhs, const Atom& rhs){
    return (lhs.MASS == rhs.MASS);
}

struct adj_conf{
    Atom a;
    Atom b;
    size_t k;
};

class root_status{
public:
    size_t color;
    size_t val;            // valence used for the fringe tree
    size_t deg;

    root_status(){}

    root_status(const size_t _color,
        const size_t _v,
        const size_t _d){
        color = _color;
        val = _v;
        deg = _d;
    }
};

bool operator<(const root_status& lhs, const root_status& rhs){
    if (lhs.val < rhs.val){
        return true;
    }
    if (lhs.val > rhs.val){
        return false;
    }
    if (lhs.color < rhs.color){
        return true;
    }
    if (lhs.color > rhs.color){
        return false;
    }
    return (lhs.deg < rhs.deg);
}


bool operator>(const root_status& lhs, const root_status& rhs){
    if (lhs.val > rhs.val){
        return true;
    }
    if (lhs.val < rhs.val){
        return false;
    }
    if (lhs.color > rhs.color){
        return true;
    }
    if (lhs.color < rhs.color){
        return false;
    }
    return (lhs.deg > rhs.deg);
}

bool operator==(const root_status& lhs, const root_status& rhs){
    return (lhs.color == rhs.color && lhs.val == rhs.val && lhs.deg == rhs.deg);
}

// A class intended to store status of a two-end path

class half_path_status{
public:
    root_status lrs;
    root_status rrs;

    half_path_status(){}

    half_path_status(
        const root_status& _l,
        const root_status& _r){
        lrs = _l;
        rrs = _r;
    }
};

bool operator==(const half_path_status& lhs, const half_path_status& rhs){
    return (lhs.lrs == rhs.lrs && lhs.rrs == rhs.rrs);
}

using Color = int;
using Mult = int;
using MultCol = pair <Mult, Color>;
using resource_vector_1D = vector <unsigned short>;

bool operator<(const MultCol& lhs, const MultCol& rhs){
    if (lhs.first < rhs.first){
        return true;
    } else if (lhs.first == rhs.first){
        return (lhs.second < rhs.second);
    } else {
        return false;
    }
}

bool operator>(const MultCol& lhs, const MultCol& rhs){
    if (lhs.first > rhs.first){
        return true;
    } else if (lhs.first == rhs.first){
        return (lhs.second > rhs.second);
    } else {
        return false;
    }
}


bool operator==(const MultCol& lhs, const MultCol& rhs){
    return (lhs.first == rhs.first && lhs.second == rhs.second);
}

class resource_vector{
public:
    vector_3D <unsigned short> CVE_in;
    vector_3D <unsigned short> CVE_out;
    vector <unsigned short> bc_num_in;
    vector <unsigned short> bc_num_out;   
    vector <unsigned short> deg_in;
    vector <unsigned short> deg_out;

    resource_vector(){
        bc_num_in = vector <unsigned short> (Dcset.size(), 0);
        bc_num_out = vector <unsigned short> (Dcset.size(), 0);
        deg_in = vector <unsigned short> (5, 0);
        deg_out = vector <unsigned short> (5, 0);
    }

    resource_vector(size_t M){
        bc_num_in = vector <unsigned short> (Dcset.size(), 0);
        bc_num_out = vector <unsigned short> (Dcset.size(), 0);
        deg_in = vector <unsigned short> (5, 0);
        deg_out = vector <unsigned short> (5, 0);
        CVE_in = vector_3D <unsigned short> (M + 1);
        for (size_t i = 0; i <= M; ++i){
            CVE_in[i] = vector_2D <unsigned short> (M + 1);
            for (size_t j = 0; j <= M; ++j){
                CVE_in[i][j] = vector <unsigned short> (4, 0);
            }
        }
        CVE_out = vector_3D <unsigned short> (M + 1);
        for (size_t i = 0; i <= M; ++i){
            CVE_out[i] = vector_2D <unsigned short> (M + 1);
            for (size_t j = 0; j <= M; ++j){
                CVE_out[i][j] = vector <unsigned short> (4, 0);
            }
        }
    }

    resource_vector(const vector_3D <unsigned short>& _CVE_in,
        const vector_3D <unsigned short>& _CVE_out,
        const vector <unsigned short>& _bc_num_in,
        const vector <unsigned short>& _bc_num_out,
        const vector <unsigned short>& _deg_in,
        const vector <unsigned short>& _deg_out
    ){
        CVE_in = _CVE_in;
        CVE_out = _CVE_out;
        bc_num_in = _bc_num_in;
        bc_num_out = _bc_num_out;
        deg_in = _deg_in;
        deg_out = _deg_out;
    }

    resource_vector(const size_t& M, const resource_vector_1D& rv1D){
        bc_num_in = vector <unsigned short> (Dcset.size(), 0);
        bc_num_out = vector <unsigned short> (Dcset.size(), 0);
        deg_in = vector <unsigned short> (5, 0);
        deg_out = vector <unsigned short> (5, 0);
        CVE_in = vector_3D <unsigned short> (M + 1);
        for (size_t i = 0; i <= M; ++i){
            CVE_in[i] = vector_2D <unsigned short> (M + 1);
            for (size_t j = 0; j <= M; ++j){
                CVE_in[i][j] = vector <unsigned short> (4, 0);
            }
        }
        CVE_out = vector_3D <unsigned short> (M + 1);
        for (size_t i = 0; i <= M; ++i){
            CVE_out[i] = vector_2D <unsigned short> (M + 1);
            for (size_t j = 0; j <= M; ++j){
                CVE_out[i][j] = vector <unsigned short> (4, 0);
            }
        }

        size_t ind = 0;
        for (size_t i = 1; i <= M; ++i){
            for (size_t j = 1; j <= i; ++j){
                for (size_t k = 0; k <= 3; ++k){
                    if (i == j && k != 0){
                        CVE_in[i][j][k] = rv1D[ind] * 2;
                    } else {
                        CVE_in[i][j][k] = rv1D[ind];
                        CVE_in[j][i][k] = rv1D[ind];
                    }
                    ++ind;
                }
            }
        }
        for (size_t i = 1; i <= M; ++i){
            for (size_t j = 1; j <= i; ++j){
                for (size_t k = 0; k <= 3; ++k){
                    if (i == j && k != 0){
                        CVE_out[i][j][k] = rv1D[ind] * 2;
                    } else {
                        CVE_out[i][j][k] = rv1D[ind];
                        CVE_out[j][i][k] = rv1D[ind];
                    }
                    ++ind;
                }
            }
        }
        for (size_t i = 0; i < Dcset.size(); ++i){
            bc_num_in[i] = rv1D[ind];
            ++ind;
        }
        for (size_t i = 0; i < Dcset.size(); ++i){
            bc_num_out[i] = rv1D[ind];
            ++ind;
        }
        for (size_t i = 1; i <=4; ++i){
            deg_in[i] = rv1D[ind];
            ++ind;
        }
        for (size_t i = 1; i <=4; ++i){
            deg_out[i] = rv1D[ind];
            ++ind;
        }
    }

    void initialize(size_t M){
        bc_num_in = vector <unsigned short> (Dcset.size(), 0);
        bc_num_out = vector <unsigned short> (Dcset.size(), 0);
        deg_in = vector <unsigned short> (5, 0);
        deg_out = vector <unsigned short> (5, 0);
        CVE_in = vector_3D <unsigned short> (M + 1);
        for (size_t i = 0; i <= M; ++i){
            CVE_in[i] = vector_2D <unsigned short> (M + 1);
            for (size_t j = 0; j <= M; ++j){
                CVE_in[i][j] = vector <unsigned short> (4, 0);
            }
        }
        CVE_out = vector_3D <unsigned short> (M + 1);
        for (size_t i = 0; i <= M; ++i){
            CVE_out[i] = vector_2D <unsigned short> (M + 1);
            for (size_t j = 0; j <= M; ++j){
                CVE_out[i][j] = vector <unsigned short> (4, 0);
            }
        }
    }

    void print(){
        size_t M = CVE_in.size() - 1;
        for (size_t i = 1; i <= M; ++i){
            for (size_t j = 1; j <= M; ++j){
                for (size_t k = 0; k <= 3; ++k){
                    cout << "in (i, j, k) = (" << i << ", " << j << ", " << k << "), " << CVE_in[i][j][k] << endl;
                }
            }
        }
        for (size_t i = 1; i <= M; ++i){
            for (size_t j = 1; j <= M; ++j){
                for (size_t k = 0; k <= 3; ++k){
                    cout << "out (i, j, k) = (" << i << ", " << j << ", " << k << "), " << CVE_out[i][j][k] << endl;
                }
            }
        }
        for (size_t i = 0; i < Dcset.size(); ++i){
            cout << "in (" << Dcset[i][0] << ", " << Dcset[i][1] << ", " << Dcset[i][2] << ") = "<< bc_num_in[i] << endl;
        }
        for (size_t i = 0; i < Dcset.size(); ++i){
            cout << "out (" << Dcset[i][0] << ", " << Dcset[i][1] << ", " << Dcset[i][2] << ") = "<< bc_num_out[i] << endl;
        }
        for (size_t i = 1; i <=4; ++i){
            cout << "in deg(" << i << ") = " << deg_in[i] <<endl;
        }
        for (size_t i = 1; i <=4; ++i){
            cout << "out deg(" << i << ") = " << deg_out[i] <<endl;
        }
        cout << endl;
    }
};

resource_vector_1D rv2rv1D(const resource_vector& rv){
    size_t M = rv.CVE_in.size() - 1;
    resource_vector_1D ans;
    for (size_t i = 1; i <= M; ++i){
        for (size_t j = 1; j <= i; ++j){
            for (size_t k = 0; k <= 3; ++k){
                if (i == j && k != 0){
                    ans.push_back(rv.CVE_in[i][j][k] / 2);
                } else {
                    ans.push_back(rv.CVE_in[i][j][k]);
                }
            }
        }
    }
    for (size_t i = 1; i <= M; ++i){
        for (size_t j = 1; j <= i; ++j){
            for (size_t k = 0; k <= 3; ++k){
                if (i == j && k != 0){
                    ans.push_back(rv.CVE_out[i][j][k] / 2);
                } else {
                    ans.push_back(rv.CVE_out[i][j][k]);
                }
            }
        }
    }
    for (size_t i = 0; i < Dcset.size(); ++i){
        ans.push_back(rv.bc_num_in[i]);
    }
    for (size_t i = 0; i < Dcset.size(); ++i){
        ans.push_back(rv.bc_num_out[i]);
    }
    for (size_t i = 1; i <= 4; ++i){
        ans.push_back(rv.deg_in[i]);
    }
    for (size_t i = 1; i <= 4; ++i){
        ans.push_back(rv.deg_out[i]);
    }
    return ans;
}

bool operator==(const resource_vector& lhs, const resource_vector& rhs){
    size_t M = lhs.CVE_in.size();

    for (size_t i = 0; i < M; ++i){
        for (size_t j = 0; j < M; ++j){
            for (size_t k = 0; k <= 3; ++k){
                if (lhs.CVE_in[i][j][k] != rhs.CVE_in[i][j][k]){
                    return false;
                }
            }
        }
    }
    for (size_t i = 0; i < M; ++i){
        for (size_t j = 0; j < M; ++j){
            for (size_t k = 0; k <= 3; ++k){
                if (lhs.CVE_out[i][j][k] != rhs.CVE_out[i][j][k]){
                    return false;
                }
            }
        }
    }
    
    for (size_t i = 0; i < Dcset.size(); ++i){
        if (lhs.bc_num_in[i] != rhs.bc_num_in[i]){
            return false;
        }
    }
    for (size_t i = 0; i < Dcset.size(); ++i){
        if (lhs.bc_num_out[i] != rhs.bc_num_out[i]){
            return false;
        }
    }

    for (size_t i = 1; i <= 4; ++i){
        if (lhs.deg_in[i] != rhs.deg_in[i]){
            return false;
        }
    }
    for (size_t i = 1; i <= 4; ++i){
        if (lhs.deg_out[i] != rhs.deg_out[i]){
            return false;
        }
    }

    return true;
}

class map_rv{
public:
    vector <unsigned short> rs;
    vector <size_t> w_ind;
    vector <unsigned short> mul;

    map_rv(){}

    map_rv(
        const unsigned short& _rs,
        const size_t& _w_ind
    ){
        rs = vector <unsigned short> ({_rs});
        w_ind = vector <size_t> ({_w_ind});
        mul.clear();
    }

    map_rv(const map_rv& _l, const map_rv& _r, const unsigned short& k, const bool ord){
        rs = _l.rs;
        w_ind = _l.w_ind;
        mul = _l.mul;

        mul.push_back(k);

        if (ord){
            rs.insert(rs.end(), _r.rs.begin(), _r.rs.end());
            w_ind.insert(w_ind.end(), _r.w_ind.begin(), _r.w_ind.end());
            mul.insert(mul.end(), _r.mul.begin(), _r.mul.end());
        } else {
            rs.insert(rs.end(), _r.rs.rbegin(), _r.rs.rend());
            w_ind.insert(w_ind.end(), _r.w_ind.rbegin(), _r.w_ind.rend());
            mul.insert(mul.end(), _r.mul.rbegin(), _r.mul.rend());
        }
        
    }

    map_rv(const map_rv& _m, const bool ord){
        rs.clear();
        w_ind.clear();
        mul.clear();

        if (ord){
            rs.insert(rs.end(), _m.rs.begin(), _m.rs.end());
            w_ind.insert(w_ind.end(), _m.w_ind.begin(), _m.w_ind.end());
            mul.insert(mul.end(), _m.mul.begin(), _m.mul.end());
        } else {
            rs.insert(rs.end(), _m.rs.rbegin(), _m.rs.rend());
            w_ind.insert(w_ind.end(), _m.w_ind.rbegin(), _m.w_ind.rend());
            mul.insert(mul.end(), _m.mul.rbegin(), _m.mul.rend());
        }
    }

};

template <class T>
bool operator<=(const vector <T>& lhs, const vector <T>& rhs){
    if (lhs.size() < rhs.size()){
        return true;
    } else if (lhs.size() > rhs.size()){
        return false;
    }
    for (size_t i = 0; i < lhs.size(); ++i){
        if (lhs[i] > rhs[i]){
            return false;
        }
    }
    return true;
}

bool operator<(const resource_vector& lhs, const resource_vector& rhs){
    if (lhs.CVE_in < rhs.CVE_in){
        return true;
    }
    if (lhs.CVE_in > rhs.CVE_in){
        return false;
    }
    if (lhs.CVE_out < rhs.CVE_out){
        return true;
    }
    if (lhs.CVE_out > rhs.CVE_out){
        return false;
    }
    if (lhs.bc_num_in < rhs.bc_num_in){
        return true;
    }
    if (lhs.bc_num_in > rhs.bc_num_in){
        return false;
    }
    if (lhs.bc_num_out < rhs.bc_num_out){
        return true;
    }
    if (lhs.bc_num_out > rhs.bc_num_out){
        return false;
    }
    if (lhs.deg_in < rhs.deg_in){
        return true;
    }
    if (lhs.deg_in > rhs.deg_in){
        return false;
    }
    if (lhs.deg_out < rhs.deg_out){
        return true;
    } else {
        return false;
    }
}

namespace std{
    template <>
    class hash <Atom>{
    public:
        size_t operator()(const Atom& aa) const {
            size_t res = std::hash<string>{} (aa.NAME);
            hash_combine<size_t>(res, aa.LABEL);
            hash_combine<size_t>(res, aa.VALENCE);
            hash_combine<size_t>(res, aa.MASS);
            return res;
        }
    };

    template <>
    class hash <root_status>{
    public:
        size_t operator()(const root_status& aa) const {
            size_t res = std::hash<size_t>{} (aa.color);
            hash_combine<size_t>(res, aa.val);
            hash_combine<size_t>(res, aa.deg);
            return res;
        }
    };

    template <>
    class hash <half_path_status>{
    public:
        size_t operator()(const half_path_status& aa) const {
            size_t res = std::hash<root_status>{}(aa.lrs);
            hash_combine<root_status>(res, aa.rrs);
            return res;
        }
    };

    template <class T>
    class hash <vector <T>>{
    public:
        size_t operator()(const vector <T>& aa) const {
            size_t res = 0;
            for (size_t i = 0; i < aa.size(); ++i){
                hash_combine<T>(res, aa[i]);
            }
            return res;
        }
    };

    template <>
    class hash <resource_vector>{
    public:
        size_t operator()(const resource_vector& aa) const {
            size_t res = std::hash <vector_3D <unsigned short>>{} (aa.CVE_in);
            hash_combine<vector_3D <unsigned short>> (res, aa.CVE_out);
            hash_combine<vector <unsigned short>> (res, aa.bc_num_in);
            hash_combine<vector <unsigned short>> (res, aa.bc_num_out);
            hash_combine<vector <unsigned short>> (res, aa.deg_in);
            hash_combine<vector <unsigned short>> (res, aa.deg_out);
            return res;
        }
    };

    template <>
    class hash <map_rv>{
    public:
        size_t operator()(const map_rv& aa) const {
            size_t res = std::hash <vector <unsigned short>>{} (aa.rs);
            hash_combine<vector <size_t>> (res, aa.w_ind);
            hash_combine<vector <unsigned short>> (res, aa.mul);
            return res;
        }
    }; 

    template <class T1, class T2>
    class hash <pair <T1, T2>>{
    public:
        size_t operator()(const pair <T1, T2>& aa) const {
            size_t res = std::hash <T1>{} (aa.first);
            hash_combine<T2> (res, aa.second);
            return res;
        }
    }; 
}

class input_info{ // input information
public:
    size_t num_kind_atoms;      //the number of the kinds of atoms
    vector < Atom > atoms;  //the information of atoms
    size_t num_kind_edges;    // the number of the kinds of edges, symmetric edges share one kind
    size_t diameter;   // the diameter of final trees, the length of the output paths 
    size_t dmax; //max degree in the generated graph
    size_t n;
    
    unordered_map <vector <size_t>, unsigned short> bc_ind_map;  // a map to store ddegree configuration.
    resource_vector rv;
    resource_vector_1D rv1D;
    size_t rv1D_size;

    vector <unsigned short> nb;

    size_t k1;  // delta1.  ,maximum possible valur
    size_t k2;  // delta2.  ,
    size_t k3;  // delta3.  ,fixed value

    input_info(){
        num_kind_atoms = 0;
        num_kind_edges = 0;
        diameter = 0;
        dmax = 0;
        atoms.clear();
        bc_ind_map.clear();
        n = 0;

        k1 = 0;
        k2 = 0;
    }

    int find_index(string& st){
        for (size_t i = 1; i <= num_kind_atoms; ++i){
            if (atoms[i].NAME == st){
                return i;
            }
        }
        return -1;
    }

    void init_bc_ind_map(){
        for (size_t i = 0; i < Dcset.size(); ++i){
            bc_ind_map[Dcset[i]] = i;
        }
    }

    void calc_nb(){
        nb = vector <unsigned short> (num_kind_atoms + 1, 0);

        for (size_t i = 1; i <= num_kind_atoms; ++i){
            for (size_t j = 1; j <= num_kind_atoms; ++j){
                for (size_t k = 1; k <= 3; ++k){
                    if (i == j){
                        nb[i] += rv.CVE_out[i][j][k];
                    } else {
                        nb[i] += rv.CVE_out[i][j][k];
                    }
                }
            }
        }
    }

    void print(){
        cout << "Atoms:" << endl;
        for (size_t i = 0; i <= num_kind_atoms; ++i){
            cout << atoms[i].NAME << " " << atoms[i].MASS << " " << atoms[i].VALENCE << endl;
        }

        rv.print();

        cout << "diameter:" << diameter << endl;
        cout << "k1:" << k1 << endl;
        cout << "k2:" << k2 << endl;
    }

};

void read_input(const string& infilename, input_info& IN_INFO) {
    //open the file
    ifstream infile;
    try{
        infile.open(infilename);
    } catch (const exception& e){
        cerr << "Couldn't open '" << infilename
            << "' for reading!" << endl;
        throw(e);
    }

    if (! infile.is_open()){
        cerr << "Error, infile not initialized!" << endl;
        throw(-1);
    }

    //read the file
    string line;
    stringstream st;

    //read the information of atoms
    getline(infile, line);
    st << line;
    st >> IN_INFO.num_kind_atoms;
    st.clear();

    IN_INFO.init_bc_ind_map();

    size_t M = IN_INFO.num_kind_atoms;
    IN_INFO.rv.initialize(M);

    Atom epsilon;
    epsilon.NAME = " ";
    epsilon.MASS = 0;
    IN_INFO.atoms.push_back(epsilon);

    size_t tmp_in = 0;

    unordered_map <string, unsigned short> atom_freq_temp_in;
    atom_freq_temp_in.clear();
    unordered_map <string, unsigned short> atom_freq_temp_out;
    atom_freq_temp_out.clear();

    for( size_t i = 0; i < M; ++i ){
        getline(infile, line); 
        Atom atom;
        st << line;
        st >> atom.NAME;
        st >> atom.MASS;
        st >> atom.VALENCE;
        st >> atom_freq_temp_in[atom.NAME];
        st >> atom_freq_temp_out[atom.NAME];

        IN_INFO.atoms.push_back(atom);
        IN_INFO.n += atom_freq_temp_in[atom.NAME] + atom_freq_temp_out[atom.NAME];
        tmp_in += atom_freq_temp_in[atom.NAME];
        st.clear();
    }

    sort(IN_INFO.atoms.begin(), IN_INFO.atoms.end());

    for (size_t i = 1; i <= M; ++i){
        IN_INFO.rv.CVE_in[i][i][0] = atom_freq_temp_in[IN_INFO.atoms[i].NAME];
        IN_INFO.rv.CVE_out[i][i][0] = atom_freq_temp_out[IN_INFO.atoms[i].NAME];
    }

    // read the information of adjacency configurations
    getline(infile, line);
    st << line;
    st >> IN_INFO.num_kind_edges;
    st.clear();
  
    for( size_t i = 0; i < IN_INFO.num_kind_edges; ++i ){
        adj_conf ac;
        getline(infile, line);
        string a, b;
        st << line;
        st >> a;
        st >> b;
        size_t label_a = IN_INFO.find_index(a);
        size_t label_b = IN_INFO.find_index(b);
        unsigned short k, freq_in, freq_out;
        st >> k;
        st >> freq_in;
        IN_INFO.rv.CVE_in[label_a][label_b][k] += freq_in;
        IN_INFO.rv.CVE_in[label_b][label_a][k] += freq_in;
        st >> freq_out;
        IN_INFO.rv.CVE_out[label_a][label_b][k] += freq_out;
        IN_INFO.rv.CVE_out[label_b][label_a][k] += freq_out;

        st.clear();
    }

    // read degree information
    for (size_t i = 1; i <= 4; ++i){
        unsigned short d_in, d_out;
        getline(infile, line);
        st << line;
        st >> d_in >> d_out;
        IN_INFO.rv.deg_in[i] = d_in;
        IN_INFO.rv.deg_out[i] = d_out;
        st.clear();
    }

    // read the diameter
    getline(infile, line);
    st << line;
    st >> IN_INFO.diameter;
    st.clear();
    // read dmax
    getline(infile, line);
    st << line;
    st >> IN_INFO.dmax;
    st.clear();

    // read bond configuration
    while (getline(infile, line)){
        st << line;
        unsigned short d1, d2, m, val_in, val_out;
        st >> d1 >> d2 >> m >> val_in >> val_out;
        IN_INFO.rv.bc_num_in[IN_INFO.bc_ind_map.at({d1, d2, m})] = val_in;
        IN_INFO.rv.bc_num_out[IN_INFO.bc_ind_map.at({d1, d2, m})] = val_out;
        st.clear();
    }

    IN_INFO.rv1D = rv2rv1D(IN_INFO.rv);
    IN_INFO.rv1D_size = IN_INFO.rv1D.size();

    IN_INFO.calc_nb();

    // IN_INFO.k2 = (IN_INFO.diameter - 6) / 2;
    // IN_INFO.k1 = IN_INFO.diameter - 6 - IN_INFO.k2;
    size_t k3 = tmp_in - IN_INFO.diameter + 3;

    if (k3 < 1){
        cerr << "Error, the input file is not for 3-paths!" << endl;
        throw(-1);
    }

    IN_INFO.k3 = k3 - 1;
    IN_INFO.k1 = IN_INFO.diameter - 5 - IN_INFO.k3;
    IN_INFO.k2 = IN_INFO.k3;

    // cout << "|Lambda_in| = " << tmp_in << " delta3 = " << IN_INFO.k3 << endl;

    infile.close();
}

size_t rv1D_ind_CVE_in(size_t a1, size_t a2, size_t k){
    if (a1 >= a2){
        return (2 * a1 * (a1 - 1) + 4 * (a2 - 1) + k);
    } else {
        return (2 * a2 * (a2 - 1) + 4 * (a1 - 1) + k);
    }
}

size_t rv1D_ind_bc_in(const input_info& IN_INFO, 
    size_t d1, size_t d2, size_t k){

    size_t M = IN_INFO.num_kind_atoms;
    size_t ans = 4 * M * (M + 1);
    if (d1 <= d2){
        ans += IN_INFO.bc_ind_map.at({d1, d2, k});
    } else {
        ans += IN_INFO.bc_ind_map.at({d2, d1, k});
    }
    return ans;
}

class FringeTree {
public:
    vector <MultCol> seq;       // Canonical sequence of multiplicity, color
    size_t num_verts;           // Number of non-eps vertices
    vector <unsigned short> degree;
    resource_vector rv;
    size_t root_val;
    size_t root_degree;

    // default constructor
    FringeTree() : seq(5), num_verts(0) {
        // pass
    }

    // constructor
    FringeTree(size_t n, size_t dmax) :
        seq(dmax + 1), num_verts(0), root_val(0) 
    {
        rv = resource_vector(n);
    }

    // copy constructor
    FringeTree(const FringeTree &rhs) :
        seq(rhs.seq),
        num_verts(rhs.num_verts),
        rv(rhs.rv),
        degree(rhs.degree),
        root_degree(rhs.root_degree),
        root_val(rhs.root_val)
    {
        // pass
    }

    void calc_degree(const size_t& dmax, const size_t& root_degree = 2){
        size_t m = seq.size();
        auto& prt = _prt[dmax];
        degree = vector <unsigned short>(m, 0);
        degree[0] = root_degree;
        for (size_t i = 1; i < m; ++i){
            if (seq[i].first != 0){
                ++degree[i];
                ++degree[prt[i]];
            }
        }
    }

    void calc_bc(const input_info& IN_INFO,
        const size_t& dmax, const size_t& root_degree = 1, const size_t& start_index = 1
    ){    
        rv.bc_num_out = vector <unsigned short> (Dcset.size(), 0);
        size_t m = seq.size();
        auto& prt = _prt[IN_INFO.dmax];
        calc_degree(dmax, root_degree);

        for (size_t i = start_index; i < m; ++i){
            if (seq[i].first != 0){
                unsigned short d1 = degree[i];
                unsigned short d2 = degree[prt[i]];
                unsigned short mul = seq[i].first;
                if (d1 <= d2){
                    ++rv.bc_num_out[IN_INFO.bc_ind_map.at({d1, d2, mul})];
                } else {
                    ++rv.bc_num_out[IN_INFO.bc_ind_map.at({d2, d1, mul})];
                }
            }
        }
    }

    size_t calc_root_val(const size_t& dmax){
        root_val = 0;
        size_t m = seq.size();
        auto& prt = _prt[dmax];
        for (size_t i = 1; i < m; ++i){
            if (prt[i] == 0){
                root_val += seq[i].first;
            }
        }  
        return root_val;
    }

    void calc_deg(const size_t& dmax, const size_t& root_degree){
        calc_degree(dmax, root_degree);
        rv.deg_in = vector <unsigned short>(5, 0);
        rv.deg_out = vector <unsigned short>(5, 0);
        size_t m = seq.size();
        ++rv.deg_in[degree[0]];
        for (size_t i = 1; i < m; ++i){
            ++rv.deg_out[degree[i]];
        }
    }

};

// A trie structure to store the weight vector
class PathTrie{
public:
    unsigned short key;  
    int parent_ind; // parent index
    vector <int> children;

    PathTrie(){}

    PathTrie(unsigned short _key) : 
        key(_key)
    {
        parent_ind = -1;
        children.clear();
    }

    PathTrie(unsigned short _key, const size_t& _parent) : 
        key(_key)
    {   
        parent_ind = _parent;
        children.clear();
    }   
};

size_t find_trie_index(
    const resource_vector_1D& rv,
    vector <PathTrie>& all_trie_node
){
    size_t depth = 0;
    size_t id = 0;

    while (depth < rv.size()){
        bool flag = false;
        int next_id = -1;
        for (auto& _id : all_trie_node[id].children){
            if (all_trie_node[_id].key == rv[depth]){
                next_id = _id;
                flag = true;
                break;
            }
        }

        if (!flag){
            all_trie_node.emplace_back(rv[depth], id);
            next_id = all_trie_node.size() - 1;
            all_trie_node[id].children.push_back(next_id);
        }

        id = next_id;
        ++depth;
    }

    return id;
}

resource_vector_1D traverse(
    size_t leaf_id,
    const vector <PathTrie>& all_trie_node,
    size_t _size
){
    resource_vector_1D tmp(_size);
    size_t id = leaf_id;
    size_t tmp_ind = _size - 1;
    while (id != 0){
        tmp[tmp_ind] = all_trie_node[id].key;
        --tmp_ind;
        id = all_trie_node[id].parent_ind;
    }

    return tmp;
}

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

    Graph _g = g;
    cutLeaf_oneLayer(_g);
    cutLeaf_oneLayer(_g);

    for (Vertex u = 0; u < g.numAtom; ++u) {
        if (_g.status[u] == true){
            internal_set.push_back(u);
        }
    }
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

#endif
