// This file mainly contains the definition and functions for the data structure
// used in the code.

#ifndef __DATA_STRUCTURE_HPP__INCLUDED
#define __DATA_STRUCTURE_HPP__INCLUDED

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "tree_signature.hpp"

using namespace std;

enum class DAG_direction_status{
  NG,
  EXACT_possible,
  LOOSE_possible,
};

// A vector to store prt information of fringe tree
const vector<vector<unsigned short>> _prt = {
    {},
    {},
    {},
    {0, 0, 1, 1, 0, 4, 4, 0, 7, 7},
    {0, 0, 1, 1, 1, 0, 5, 5, 5, 0, 9, 9, 9}};

template <class T>
using vector_2D = vector<vector<T>>;

template <class T>
using vector_3D = vector<vector_2D<T>>;

template <class T>
using vector_4D = vector<vector_3D<T>>;

template <class T>
inline void hash_combine(std::size_t& seed, const T& v) {
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct Atom {
  string NAME;     // the chemical name of an atom
  size_t LABEL;    // the index of an atom
  int VALENCE;  // valence
  size_t MASS;
};

bool operator<(const Atom& lhs, const Atom& rhs) {
  return (lhs.MASS < rhs.MASS);
}

bool operator>(const Atom& lhs, const Atom& rhs) {
  return (lhs.MASS > rhs.MASS);
}

bool operator==(const Atom& lhs, const Atom& rhs) {
  return (lhs.MASS == rhs.MASS);
}

class atom_degree {
 public:
  size_t color;
  size_t deg;

  atom_degree() {
    color = 0;
    deg = 0;
  }

  atom_degree(size_t _color, size_t _deg) {
    color = _color;
    deg = _deg;
  }
};

bool operator<(const atom_degree& lhs, const atom_degree& rhs) {
  return (lhs.color < rhs.color ||
          (lhs.color == rhs.color && lhs.deg < rhs.deg));
}

bool operator>(const atom_degree& lhs, const atom_degree& rhs) {
  return (lhs.color > rhs.color ||
          (lhs.color == rhs.color && lhs.deg > rhs.deg));
}

bool operator==(const atom_degree& lhs, const atom_degree& rhs) {
  return (lhs.color == rhs.color && lhs.deg == rhs.deg);
}

class rv_edge {
 public:
  string type;
  atom_degree ad1;
  atom_degree ad2;
  size_t mul;

  rv_edge() {
    type = "";
    ad1 = atom_degree(0, 0);
    ad2 = atom_degree(0, 0);
    mul = 0;
  }

  rv_edge(string _type, const atom_degree& _ad1, const atom_degree& _ad2,
          const size_t& _mul) {
    type = _type;
    ad1 = _ad1;
    ad2 = _ad2;
    mul = _mul;
  }

  bool existence_check(size_t va, size_t vb) {
    return (mul <= min(va + 1 - ad1.deg, vb + 1 - ad2.deg));
  }
};

bool operator==(const rv_edge& lhs, const rv_edge& rhs) {
  return (lhs.type == rhs.type && lhs.ad1 == rhs.ad1 && lhs.ad2 == rhs.ad2 &&
          lhs.mul == rhs.mul);
}

class root_status {
 public:
  size_t color;
  int val;  // valence used for the fringe tree
  size_t deg;

  root_status() {}

  root_status(const size_t _color, const int _v, const size_t _d) {
    color = _color;
    val = _v;
    deg = _d;
  }
};

bool operator<(const root_status& lhs, const root_status& rhs) {
  if (lhs.val < rhs.val) {
    return true;
  }
  if (lhs.val > rhs.val) {
    return false;
  }
  if (lhs.color < rhs.color) {
    return true;
  }
  if (lhs.color > rhs.color) {
    return false;
  }
  return (lhs.deg < rhs.deg);
}

bool operator>(const root_status& lhs, const root_status& rhs) {
  if (lhs.val > rhs.val) {
    return true;
  }
  if (lhs.val < rhs.val) {
    return false;
  }
  if (lhs.color > rhs.color) {
    return true;
  }
  if (lhs.color < rhs.color) {
    return false;
  }
  return (lhs.deg > rhs.deg);
}

bool operator==(const root_status& lhs, const root_status& rhs) {
  return (lhs.color == rhs.color && lhs.val == rhs.val && lhs.deg == rhs.deg);
}

bool operator!=(const root_status& lhs, const root_status& rhs) {
  return (lhs.color != rhs.color || lhs.val != rhs.val || lhs.deg != rhs.deg);
}

class root_status_height {
 public:
  size_t color;
  int val;  // valence used for the fringe tree
  size_t deg;
  size_t height;

  root_status_height() {}

  root_status_height(const size_t _color, const int _v, const size_t _d,
                     const size_t _h) {
    color = _color;
    val = _v;
    deg = _d;
    height = _h;
  }

  root_status_height(const root_status& _rs, const size_t _h) {
    color = _rs.color;
    val = _rs.val;
    deg = _rs.deg;
    height = _h;
  }
};

bool operator<(const root_status_height& lhs, const root_status_height& rhs) {
  if (lhs.val < rhs.val) {
    return true;
  }
  if (lhs.val > rhs.val) {
    return false;
  }
  if (lhs.color < rhs.color) {
    return true;
  }
  if (lhs.color > rhs.color) {
    return false;
  }
  if (lhs.deg < rhs.deg) {
    return true;
  }
  if (lhs.deg > rhs.deg) {
    return false;
  }
  return (lhs.height < rhs.height);
}

bool operator>(const root_status_height& lhs, const root_status_height& rhs) {
  if (lhs.val > rhs.val) {
    return true;
  }
  if (lhs.val < rhs.val) {
    return false;
  }
  if (lhs.color > rhs.color) {
    return true;
  }
  if (lhs.color < rhs.color) {
    return false;
  }
  if (lhs.deg > rhs.deg) {
    return true;
  }
  if (lhs.deg < rhs.deg) {
    return false;
  }
  return (lhs.height > rhs.height);
}

bool operator==(const root_status_height& lhs, const root_status_height& rhs) {
  return (lhs.color == rhs.color && lhs.val == rhs.val && lhs.deg == rhs.deg &&
          lhs.height == rhs.height);
}

// A class intended to store status of a two-end path

class half_path_status {
 public:
  root_status lrs;
  root_status rrs;

  half_path_status() {}

  half_path_status(const root_status& _l, const root_status& _r) {
    lrs = _l;
    rrs = _r;
  }
};

bool operator==(const half_path_status& lhs, const half_path_status& rhs) {
  return (lhs.lrs == rhs.lrs && lhs.rrs == rhs.rrs);
}

using Color = int;
using Mult = int;
using MultCol = pair<Mult, Color>;
using resource_vector_1D = vector<unsigned short>;
using root_status_T = pair<root_status, root_status_height>;

bool operator<(const MultCol& lhs, const MultCol& rhs) {
  if (lhs.first < rhs.first) {
    return true;
  } else if (lhs.first == rhs.first) {
    return (lhs.second < rhs.second);
  } else {
    return false;
  }
}

bool operator>(const MultCol& lhs, const MultCol& rhs) {
  if (lhs.first > rhs.first) {
    return true;
  } else if (lhs.first == rhs.first) {
    return (lhs.second > rhs.second);
  } else {
    return false;
  }
}

bool operator==(const MultCol& lhs, const MultCol& rhs) {
  return (lhs.first == rhs.first && lhs.second == rhs.second);
}

// A class to store vectors in a not so compact way
// Note: this class only stores additive variable
class resource_vector {
 public:
  size_t bl;
  // vector_3D<unsigned short> ec_co;
  vector_3D<unsigned short> ec_int;
  // vector_3D<unsigned short> ec_ex;

  vector <unsigned short> fc;

  resource_vector() {
    bl = 0;
    // ec_co.clear();
    ec_int.clear();
    // ec_ex.clear();
    fc.clear();
  }

  resource_vector(size_t M, size_t num_fc) {  // Here M is supposed to the number of kinds of
                               // (atom, degree)
    bl = 0;
    // ec_co = vector_3D<unsigned short>(M);
    // for (size_t i = 0; i < M; ++i) {
    //   ec_co[i] = vector_2D<unsigned short>(M);
    //   for (size_t j = 0; j < M; ++j) {
    //     ec_co[i][j] = vector<unsigned short>(4, 0);
    //   }
    // }
    ec_int = vector_3D<unsigned short>(M);
    for (size_t i = 0; i < M; ++i) {
      ec_int[i] = vector_2D<unsigned short>(M);
      for (size_t j = 0; j < M; ++j) {
        ec_int[i][j] = vector<unsigned short>(4, 0);
      }
    }
    // ec_ex = vector_3D<unsigned short>(M);
    // for (size_t i = 0; i < M; ++i) {
    //   ec_ex[i] = vector_2D<unsigned short>(M);
    //   for (size_t j = 0; j < M; ++j) {
    //     ec_ex[i][j] = vector<unsigned short>(4, 0);
    //   }
    // }

    fc = vector <unsigned short>(num_fc, 0);
  }

  resource_vector(const size_t& _bl, 
    // const vector_3D<unsigned short>& _ec_co,
                  const vector_3D<unsigned short>& _ec_int,
                  // const vector_3D<unsigned short>& _ec_ex
                  const vector <unsigned short>& _fc
                  ) {
    bl = _bl;
    // ec_co = _ec_co;
    ec_int = _ec_int;
    // ec_ex = _ec_ex;
    fc = _fc;
  }

  void initialize_ec(size_t M) {
    bl = 0;
    // ec_co = vector_3D<unsigned short>(M);
    // for (size_t i = 0; i < M; ++i) {
    //   ec_co[i] = vector_2D<unsigned short>(M);
    //   for (size_t j = 0; j < M; ++j) {
    //     ec_co[i][j] = vector<unsigned short>(4, 0);
    //   }
    // }
    ec_int = vector_3D<unsigned short>(M);
    for (size_t i = 0; i < M; ++i) {
      ec_int[i] = vector_2D<unsigned short>(M);
      for (size_t j = 0; j < M; ++j) {
        ec_int[i][j] = vector<unsigned short>(4, 0);
      }
    }
    // ec_ex = vector_3D<unsigned short>(M);
    // for (size_t i = 0; i < M; ++i) {
    //   ec_ex[i] = vector_2D<unsigned short>(M);
    //   for (size_t j = 0; j < M; ++j) {
    //     ec_ex[i][j] = vector<unsigned short>(4, 0);
    //   }
    // }
    
  }

  void initialize_fc(size_t num_fc){
    fc = vector <unsigned short>(num_fc, 0);
  }
};

// resource_vector_1D rv2rv1D(const resource_vector& rv){
//     size_t M = rv.ec_co.size();
//     resource_vector_1D ans;
//     ans.push_back(rv.bl);
//     for (size_t i = 0; i < M; ++i){
//         for (size_t j = 0; j <= i; ++j){
//             for (size_t k = 1; k <= 3; ++k){
//                 if (i == j){
//                     ans.push_back(rv.ec_co[i][j][k] / 2);
//                 } else {
//                     ans.push_back(rv.ec_co[i][j][k]);
//                 }
//             }
//         }
//     }
//     for (size_t i = 0; i < M; ++i){
//         for (size_t j = 0; j < M; ++j){
//             for (size_t k = 1; k <= 3; ++k){
//                 ans.push_back(rv.ec_in[i][j][k]);
//             }
//         }
//     }
//     for (size_t i = 0; i < M; ++i){
//         for (size_t j = 0; j < M; ++j){
//             for (size_t k = 1; k <= 3; ++k){
//                 ans.push_back(rv.ec_ex[i][j][k]);
//             }
//         }
//     }

//     return ans;
// }

class map_rv {
 public:
  vector<root_status> rs;
  vector<size_t> w_ind;
  vector<unsigned short> mul;

  size_t height;

  map_rv() {
    rs.clear();
    w_ind.clear();
    mul.clear();
    height = 0;
  }

  map_rv(const root_status& _rs, const size_t& _w_ind) {
    rs = vector<root_status>({_rs});
    w_ind = vector<size_t>({_w_ind});
    mul.clear();
    height = 0;
  }

  map_rv(const root_status_height& _rsh, const size_t& _w_ind) {
    root_status _rs(_rsh.color, _rsh.val, _rsh.deg);
    rs = vector<root_status>({_rs});
    w_ind = vector<size_t>({_w_ind});
    mul.clear();
    height = _rsh.height;
  }

  map_rv(const root_status& _rs, const size_t& _w_ind,
         const unsigned short& _mul) {
    rs = vector<root_status>({_rs});
    w_ind = vector<size_t>({_w_ind});
    mul = vector<unsigned short>({_mul});
    height = 0;
  }

  map_rv(const map_rv& _l, const map_rv& _r, const unsigned short& k,
         const bool ord) {
    rs = _l.rs;
    w_ind = _l.w_ind;
    mul = _l.mul;
    height = 0;

    mul.push_back(k);

    if (ord) {
      rs.insert(rs.end(), _r.rs.begin(), _r.rs.end());
      w_ind.insert(w_ind.end(), _r.w_ind.begin(), _r.w_ind.end());
      mul.insert(mul.end(), _r.mul.begin(), _r.mul.end());
    } else {
      rs.insert(rs.end(), _r.rs.rbegin(), _r.rs.rend());
      w_ind.insert(w_ind.end(), _r.w_ind.rbegin(), _r.w_ind.rend());
      mul.insert(mul.end(), _r.mul.rbegin(), _r.mul.rend());
    }
  }

  map_rv(const map_rv& _l, const map_rv& _r, const bool ord) {
    rs = _l.rs;
    w_ind = _l.w_ind;
    mul = _l.mul;
    height = 0;

    if (ord) {
      rs.insert(rs.end(), _r.rs.begin(), _r.rs.end());
      w_ind.insert(w_ind.end(), _r.w_ind.begin(), _r.w_ind.end());
      mul.insert(mul.end(), _r.mul.begin(), _r.mul.end());
    } else {
      rs.insert(rs.end(), _r.rs.rbegin(), _r.rs.rend());
      w_ind.insert(w_ind.end(), _r.w_ind.rbegin(), _r.w_ind.rend());
      mul.insert(mul.end(), _r.mul.rbegin(), _r.mul.rend());
    }
  }

  map_rv(const map_rv& _m, const bool ord) {
    rs.clear();
    w_ind.clear();
    mul.clear();
    height = 0;

    if (ord) {
      rs.insert(rs.end(), _m.rs.begin(), _m.rs.end());
      w_ind.insert(w_ind.end(), _m.w_ind.begin(), _m.w_ind.end());
      mul.insert(mul.end(), _m.mul.begin(), _m.mul.end());
    } else {
      rs.insert(rs.end(), _m.rs.rbegin(), _m.rs.rend());
      w_ind.insert(w_ind.end(), _m.w_ind.rbegin(), _m.w_ind.rend());
      mul.insert(mul.end(), _m.mul.rbegin(), _m.mul.rend());
    }
  }

  void extend(const map_rv& _m, const bool ord) {
    if (ord) {
      rs.insert(rs.end(), _m.rs.begin(), _m.rs.end());
      w_ind.insert(w_ind.end(), _m.w_ind.begin(), _m.w_ind.end());
      mul.insert(mul.end(), _m.mul.begin(), _m.mul.end());
    } else {
      rs.insert(rs.end(), _m.rs.rbegin(), _m.rs.rend());
      w_ind.insert(w_ind.end(), _m.w_ind.rbegin(), _m.w_ind.rend());
      mul.insert(mul.end(), _m.mul.rbegin(), _m.mul.rend());
    }
  }

  // void print() {
  //   cout << "map_rv" << endl;
  //   cout << "rs : ";
  //   for (auto& a : rs) {
  //     cout << a << " ";
  //   }
  //   cout << endl;
  //   cout << "w_ind : ";
  //   for (auto& a : w_ind) {
  //     cout << a << " ";
  //   }
  //   cout << endl;
  //   cout << "mul : ";
  //   for (auto& a : mul) {
  //     cout << a << " ";
  //   }
  //   cout << endl;
  // }
};

template <class T>
bool operator<=(const vector<T>& lhs, const vector<T>& rhs) {
  if (lhs.size() < rhs.size()) {
    return true;
  } else if (lhs.size() > rhs.size()) {
    return false;
  }
  for (size_t i = 0; i < lhs.size(); ++i) {
    if (lhs[i] > rhs[i]) {
      return false;
    }
  }
  return true;
}

class map_tree {
 public:
  vector<map_rv> trees;
  vector<unsigned short> mul;
  size_t height;

  map_tree() {
    trees.clear();
    mul.clear();
    height = 0;
  }

  map_tree(const map_rv& _trees) {
    trees = vector<map_rv> {(_trees)};
    mul.clear();
    height = 0;
  }

  map_tree(const map_rv& _trees, const size_t& _h) {
    trees = vector<map_rv> {(_trees)};
    mul.clear();
    height = _h;
  }

  map_tree(const map_tree& _l, const map_tree& _r, const unsigned short& k,
           const size_t& _h, const bool ord) {
    trees = _l.trees;
    mul = _l.mul;
    height = _h;

    mul.push_back(k);

    if (ord) {
      trees.insert(trees.end(), _r.trees.begin(), _r.trees.end());
      mul.insert(mul.end(), _r.mul.begin(), _r.mul.end());
    } else {
      trees.insert(trees.end(), _r.trees.rbegin(), _r.trees.rend());
      mul.insert(mul.end(), _r.mul.rbegin(), _r.mul.rend());
    }
  }

  map_tree(const map_tree& _l, const map_tree& _r, const size_t& _h,
           const bool ord) {
    trees = _l.trees;
    mul = _l.mul;
    height = _h;

    if (ord) {
      trees.insert(trees.end(), _r.trees.begin(), _r.trees.end());
      mul.insert(mul.end(), _r.mul.begin(), _r.mul.end());
    } else {
      trees.insert(trees.end(), _r.trees.rbegin(), _r.trees.rend());
      mul.insert(mul.end(), _r.mul.rbegin(), _r.mul.rend());
    }
  }

  map_tree(const map_tree& _m, const bool ord) {
    trees.clear();
    mul.clear();
    height = _m.height;

    if (ord) {
      trees.insert(trees.end(), _m.trees.begin(), _m.trees.end());
      mul.insert(mul.end(), _m.mul.begin(), _m.mul.end());
    } else {
      trees.insert(trees.end(), _m.trees.rbegin(), _m.trees.rend());
      mul.insert(mul.end(), _m.mul.rbegin(), _m.mul.rend());
    }
  }

  void extend(const map_tree& _m, const bool ord) {
    if (_m.height > height) {
      height = _m.height;
    }

    if (ord) {
      trees.insert(trees.end(), _m.trees.begin(), _m.trees.end());
      mul.insert(mul.end(), _m.mul.begin(), _m.mul.end());
    } else {
      trees.insert(trees.end(), _m.trees.rbegin(), _m.trees.rend());
      mul.insert(mul.end(), _m.mul.rbegin(), _m.mul.rend());
    }
  }

  // void print() {
  //   cout << "map_tree" << endl;
  //   cout << "trees : ";
  //   for (auto& a : trees) {
  //     a.print();
  //   }
  //   cout << endl;
  //   cout << "tree mul : ";
  //   for (auto& a : mul) {
  //     cout << a << " ";
  //   }
  //   cout << endl;
  // }
};

class DAG_node {
 public:
  root_status rs;
  size_t l;
  size_t ind;

  DAG_node() {}

  DAG_node(const root_status& _rs, const size_t _l, const size_t& _ind) {
    rs = _rs;
    l = _l;
    ind = _ind;
  }
};

class DAG_node_T {
 public:
  root_status_T rsT;
  size_t l;
  size_t ind;

  DAG_node_T() {}

  DAG_node_T(const root_status_T& _rsT, const size_t _l, const size_t& _ind) {
    rsT = _rsT;
    l = _l;
    ind = _ind;
  }
};

class DAG_node_height {
 public:
  root_status_height rsh;
  size_t l;
  size_t ind;

  DAG_node_height() {}

  DAG_node_height(const root_status_height& _rsh, const size_t _l,
                  const size_t& _ind) {
    rsh = _rsh;
    l = _l;
    ind = _ind;
  }
};

bool operator==(const DAG_node& lhs, const DAG_node& rhs) {
  return (lhs.rs == rhs.rs && lhs.l == rhs.l && lhs.ind == rhs.ind);
}

bool operator==(const DAG_node_T& lhs, const DAG_node_T& rhs) {
  return (lhs.rsT == rhs.rsT && lhs.l == rhs.l && lhs.ind == rhs.ind);
}

bool operator==(const DAG_node_height& lhs, const DAG_node_height& rhs) {
  return (lhs.rsh == rhs.rsh && lhs.l == rhs.l && lhs.ind == rhs.ind);
}

namespace std {
template <>
class hash<Atom> {
 public:
  size_t operator()(const Atom& aa) const {
    size_t res = std::hash<string> {}(aa.NAME);
    hash_combine<size_t>(res, aa.LABEL);
    hash_combine<size_t>(res, aa.VALENCE);
    hash_combine<size_t>(res, aa.MASS);
    return res;
  }
};

template <>
class hash<atom_degree> {
 public:
  size_t operator()(const atom_degree& aa) const {
    size_t res = 0;
    hash_combine<size_t>(res, aa.color);
    hash_combine<size_t>(res, aa.deg);
    return res;
  }
};

template <>
class hash<root_status> {
 public:
  size_t operator()(const root_status& aa) const {
    size_t res = std::hash<size_t> {}(aa.color);
    hash_combine<int>(res, aa.val);
    hash_combine<size_t>(res, aa.deg);
    return res;
  }
};

template <>
class hash<root_status_height> {
 public:
  size_t operator()(const root_status_height& aa) const {
    size_t res = std::hash<size_t> {}(aa.color);
    hash_combine<int>(res, aa.val);
    hash_combine<size_t>(res, aa.deg);
    hash_combine<size_t>(res, aa.height);
    return res;
  }
};

template <>
class hash<half_path_status> {
 public:
  size_t operator()(const half_path_status& aa) const {
    size_t res = std::hash<root_status> {}(aa.lrs);
    hash_combine<root_status>(res, aa.rrs);
    return res;
  }
};

template <class T>
class hash<vector<T>> {
 public:
  size_t operator()(const vector<T>& aa) const {
    size_t res = 0;
    for (size_t i = 0; i < aa.size(); ++i) {
      hash_combine<T>(res, aa[i]);
    }
    return res;
  }
};

template <>
class hash<DAG_node> {
 public:
  size_t operator()(const DAG_node& aa) const {
    size_t res = 0;
    hash_combine<root_status>(res, aa.rs);
    hash_combine<size_t>(res, aa.l);
    hash_combine<size_t>(res, aa.ind);
    return res;
  }
};

template <>
class hash<DAG_node_T> {
 public:
  size_t operator()(const DAG_node_T& aa) const {
    size_t res = 0;
    hash_combine<root_status_T>(res, aa.rsT);
    hash_combine<size_t>(res, aa.l);
    hash_combine<size_t>(res, aa.ind);
    return res;
  }
};

template <>
class hash<DAG_node_height> {
 public:
  size_t operator()(const DAG_node_height& aa) const {
    size_t res = 0;
    hash_combine<root_status_height>(res, aa.rsh);
    hash_combine<size_t>(res, aa.l);
    hash_combine<size_t>(res, aa.ind);
    return res;
  }
};

template <>
class hash<map_rv> {
 public:
  size_t operator()(const map_rv& aa) const {
    size_t res = std::hash<vector<root_status>> {}(aa.rs);
    hash_combine<vector<size_t>>(res, aa.w_ind);
    hash_combine<vector<unsigned short>>(res, aa.mul);
    hash_combine<size_t>(res, aa.height);
    return res;
  }
};

template <>
class hash<map_tree> {
 public:
  size_t operator()(const map_tree& aa) const {
    size_t res = std::hash<vector<map_rv>> {}(aa.trees);
    hash_combine<vector<unsigned short>>(res, aa.mul);
    hash_combine<size_t>(res, aa.height);
    return res;
  }
};

template <class T1, class T2>
class hash<pair<T1, T2>> {
 public:
  size_t operator()(const pair<T1, T2>& aa) const {
    size_t res = std::hash<T1> {}(aa.first);
    hash_combine<T2>(res, aa.second);
    return res;
  }
};

template <>
class hash<rv_edge> {
 public:
  size_t operator()(const rv_edge& aa) const {
    size_t res = 0;
    hash_combine<string>(res, aa.type);
    hash_combine<atom_degree>(res, aa.ad1);
    hash_combine<atom_degree>(res, aa.ad2);
    hash_combine<size_t>(res, aa.mul);
    return res;
  }
};
}  // namespace std

class input_info {  // component information
 public:
  size_t num_kind_atoms;        // the number of the kinds of atoms
  size_t num_kind_atom_degree;  // the number of the kinds of (atom, deg)
  size_t num_kind_rv_edge;
  vector<Atom> atoms;  // the information of atoms
  vector<atom_degree> atom_deg;
  vector<rv_edge> rv_edges;
  size_t dmax;    // max degree in the generated graph
  size_t n;       // the number of  vertices in total in this component
  size_t height;  // the  height of this component
  size_t k_star;  // 2
  size_t chLB;
  size_t chUB;
  size_t delta_1;  //  the extra degree in the core part
  size_t delta_2;  //  the extra degree in the core part, only used when this
                   //  component is edge component

  size_t num_fc;   //  number of fringe configuration
  map <size_t, size_t> fc_map;  // map index in TS_map to the index in the rv

  unordered_map<atom_degree, size_t>
      atom_deg_map;  // a map to get the index of a (atom, deg),
  unordered_map<rv_edge, size_t> rv_edge_map;

  resource_vector rv;
  resource_vector_1D rv1D;
  resource_vector_1D rv1D_loose;
  size_t rv1D_size;

  vector<unsigned short> nb;

  root_status target_rs1;  // vertex information
  root_status target_rs2;  // vertex information, only used when this component
                           // is edge component

  size_t rs1_index;
  size_t rs2_index;

  vector<size_t> base_vertices;
  vector <size_t> fringe_tree_indices;

  // unordered_set<string> Gamma_co;
  unordered_set<string> Gamma_int;
  // unordered_set<string> Gamma_ex;

  input_info() {
    num_kind_atoms = 0;
    num_kind_atom_degree = 0;
    num_kind_rv_edge = 0;
    dmax = 3;
    atoms.clear();
    atom_deg.clear();
    atom_deg_map.clear();

    base_vertices.clear();
    fringe_tree_indices.clear();

    n = 0;
    k_star = 2;
    height = 0;
    chLB = 0;
    chUB = 0;

    rs1_index = 0;
    rs2_index = 0;

    num_fc = 0;
    fc_map.clear();

    // Gamma_co.clear();
    Gamma_int.clear();
    // Gamma_ex.clear();
  }

  void initialize_atom_deg() {
    atom_deg.clear();
    atom_deg_map.clear();
    num_kind_atom_degree = 0;

    for (size_t i = 0; i < atoms.size(); ++i) {
      for (size_t d = 1; d <= atoms[i].VALENCE + 3; ++d) {
        if (d > 4) continue;
        atom_degree tmp(i, d);
        atom_deg.push_back(tmp);
        atom_deg_map.emplace(tmp, num_kind_atom_degree);
        ++num_kind_atom_degree;
      }
    }
  }

  void initialize_rv_edge_map() {
    size_t M = num_kind_atom_degree;
    rv_edges.clear();
    rv_edge_map.clear();
    rv_edge null;
    rv_edges.push_back(null);

    num_kind_rv_edge = 1;

    // for (size_t i = 0; i < M; ++i) {
    //   for (size_t j = 0; j <= i; ++j) {
    //     for (size_t k = 1; k <= 3; ++k) {
    //       atom_degree ad1 = atom_deg[i];
    //       atom_degree ad2 = atom_deg[j];
    //       if (ad1.deg < 2 || ad2.deg < 2) continue;

    //       rv_edge tmp("co", ad1, ad2, k);
    //       if (tmp.existence_check(atoms[ad1.color].VALENCE,
    //                               atoms[ad2.color].VALENCE)) {
    //         rv_edges.push_back(tmp);
    //         rv_edge_map.emplace(tmp, num_kind_rv_edge);
    //         ++num_kind_rv_edge;
    //       }
    //     }
    //   }
    // }
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j <= i; ++j) {
        for (size_t k = 1; k <= 3; ++k) {
          atom_degree ad1 = atom_deg[i];
          atom_degree ad2 = atom_deg[j];
          if (ad1.deg < 2 || ad2.deg < 2) continue;

          rv_edge tmp("int", ad1, ad2, k);
          if (tmp.existence_check(atoms[ad1.color].VALENCE + 3,
                                  atoms[ad2.color].VALENCE + 3)) {
            rv_edges.push_back(tmp);
            rv_edge_map.emplace(tmp, num_kind_rv_edge);
            ++num_kind_rv_edge;
          }
        }
      }
    }
    // for (size_t i = 0; i < M; ++i) {
    //   for (size_t j = 0; j < M; ++j) {
    //     for (size_t k = 1; k <= 3; ++k) {
    //       atom_degree ad1 = atom_deg[i];
    //       atom_degree ad2 = atom_deg[j];
    //       if (ad1.deg == 1 && ad2.deg == 1) continue;

    //       rv_edge tmp("ex", ad1, ad2, k);
    //       if (tmp.existence_check(atoms[ad1.color].VALENCE,
    //                               atoms[ad2.color].VALENCE)) {
    //         rv_edges.push_back(tmp);
    //         rv_edge_map.emplace(tmp, num_kind_rv_edge);
    //         ++num_kind_rv_edge;
    //       }
    //     }
    //   }
    // }
    rv1D_size = num_kind_rv_edge;
  }

  int find_index(const string& st) const {
    for (size_t i = 0; i < num_kind_atoms; ++i) {
      if (atoms[i].NAME == st) {
        return i;
      }
    }
    return -1;
  }

  void calculate_Gamma() {
    size_t M = num_kind_atom_degree;
    // for (size_t i = 0; i < M; ++i) {
    //   for (size_t j = 0; j <= i; ++j) {
    //     for (size_t k = 1; k <= 3; ++k) {
    //       if (rv.ec_co[i][j][k] != 0) {
    //         string tmp = "";
    //         atom_degree ad1 = atom_deg[i];
    //         tmp += atoms[ad1.color].NAME + to_string(ad1.deg);
    //         atom_degree ad2 = atom_deg[j];
    //         tmp += atoms[ad2.color].NAME + to_string(ad2.deg);
    //         tmp += to_string(k);
    //         Gamma_co.insert(tmp);
    //       }
    //     }
    //   }
    // }
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < M; ++j) {
        for (size_t k = 1; k <= 3; ++k) {
          if (rv.ec_int[i][j][k] != 0) {
            string tmp = "";
            atom_degree ad1 = atom_deg[i];
            tmp += atoms[ad1.color].NAME + to_string(ad1.deg);
            atom_degree ad2 = atom_deg[j];
            tmp += atoms[ad2.color].NAME + to_string(ad2.deg);
            tmp += to_string(k);
            Gamma_int.insert(tmp);
          }
        }
      }
    }
    // for (size_t i = 0; i < M; ++i) {
    //   for (size_t j = 0; j < M; ++j) {
    //     for (size_t k = 1; k <= 3; ++k) {
    //       if (rv.ec_ex[i][j][k] != 0) {
    //         string tmp = "";
    //         atom_degree ad1 = atom_deg[i];
    //         tmp += atoms[ad1.color].NAME + to_string(ad1.deg);
    //         atom_degree ad2 = atom_deg[j];
    //         tmp += atoms[ad2.color].NAME + to_string(ad2.deg);
    //         tmp += to_string(k);
    //         Gamma_ex.insert(tmp);
    //       }
    //     }
    //   }
    // }
  }

  void gen_rv1D_loose(const vector <unsigned short>& fc_all){
    rv1D_loose = rv1D;
    for (size_t i = num_kind_rv_edge; i < rv1D_size; ++i){
      rv1D_loose[i] = fc_all[i - num_kind_rv_edge];
    }
  }

  // need to modify
  // void calc_nb(){
  //     nb = vector <unsigned short> (num_kind_atoms + 1, 0);

  //     for (size_t i = 1; i <= num_kind_atoms; ++i){
  //         for (size_t j = 1; j <= num_kind_atoms; ++j){
  //             for (size_t k = 1; k <= 3; ++k){
  //                 if (i == j){
  //                     nb[i] += rv.CVE_out[i][j][k];
  //                 } else {
  //                     nb[i] += rv.CVE_out[i][j][k];
  //                 }
  //             }
  //         }
  //     }
  // }
};

resource_vector_1D rv2rv1D(const resource_vector& rv,
                           const input_info& IN_INFO) {
  resource_vector_1D ans(IN_INFO.rv1D_size, 0);
  size_t M = IN_INFO.num_kind_atom_degree;
  ans[0] = rv.bl;
  // for (size_t i = 0; i < M; ++i) {
  //   for (size_t j = 0; j <= i; ++j) {
  //     for (size_t k = 1; k <= 3; ++k) {
  //       rv_edge tmp("co", IN_INFO.atom_deg[i], IN_INFO.atom_deg[j], k);
  //       if (IN_INFO.rv_edge_map.find(tmp) != IN_INFO.rv_edge_map.end()) {
  //         size_t ind = IN_INFO.rv_edge_map.at(tmp);
  //         if (i == j) {
  //           ans[ind] = rv.ec_co[i][j][k] / 2;
  //         } else {
  //           ans[ind] = rv.ec_co[i][j][k];
  //         }
  //       }
  //     }
  //   }
  // }
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      for (size_t k = 1; k <= 3; ++k) {
        rv_edge tmp("int", IN_INFO.atom_deg[i], IN_INFO.atom_deg[j], k);
        if (IN_INFO.rv_edge_map.find(tmp) != IN_INFO.rv_edge_map.end()){
          size_t ind = IN_INFO.rv_edge_map.at(tmp);
          if (i == j){
            ans[ind] = rv.ec_int[i][j][k] / 2;
          } else {
            ans[ind] = rv.ec_int[i][j][k];
          }
        }
      }
    }
  }
  // for (size_t i = 0; i < M; ++i) {
  //   for (size_t j = 0; j < M; ++j) {
  //     for (size_t k = 1; k <= 3; ++k) {
  //       rv_edge tmp("ex", IN_INFO.atom_deg[i], IN_INFO.atom_deg[j], k);
  //       if (IN_INFO.rv_edge_map.find(tmp) != IN_INFO.rv_edge_map.end()) {
  //         size_t ind = IN_INFO.rv_edge_map.at(tmp);
  //         ans[ind] = rv.ec_ex[i][j][k];
  //       }
  //     }
  //   }
  // }
  size_t ind = IN_INFO.num_kind_rv_edge;
  for (size_t i = 0; i < IN_INFO.num_fc; ++i){
    ans[ind] = rv.fc[i];
    ++ind;
  }
  return ans;
}

// need to modify
// size_t rv1D_ind_CVE_in(size_t a1, size_t a2, size_t k){
//     if (a1 >= a2){
//         return (2 * a1 * (a1 - 1) + 4 * (a2 - 1) + k);
//     } else {
//         return (2 * a2 * (a2 - 1) + 4 * (a1 - 1) + k);
//     }
// }

class FringeTree {
 public:
  vector<MultCol> seq;  // Canonical sequence of multiplicity, color
  size_t num_verts;     // Number of non-eps vertices
  vector<unsigned short> degree;
  resource_vector rv;
  int root_val;
  size_t root_degree;
  size_t height;

  vector <int> charge;

  // default constructor
  FringeTree() : seq(5), num_verts(0) {
    // pass
  }

  // constructor
  // FringeTree(size_t n, size_t dmax)
  //     : seq(dmax + 1), num_verts(0), root_val(0), height(0) {
  //   rv = resource_vector(n);
  // }

  // copy constructor
  // from a sequence of fringe tree H-unsuppressed get a
  // fringe tree H-suppressed
  FringeTree(const FringeTree& rhs)
      : seq(rhs.seq),
        num_verts(rhs.num_verts),
        rv(rhs.rv),
        degree(rhs.degree),
        root_degree(rhs.root_degree),
        root_val(rhs.root_val),
        height(rhs.height),
        charge(rhs.charge) {
    // pass
  }

  // constructor from TS
  FringeTree(
    const input_info& IN_INFO, 
    const TreeSignature& TS, 
    const size_t& dmax
  ){
    seq.clear();
    num_verts = TS.mu.size() + 1;
    root_degree = 0;
    height = 0;

    seq.emplace_back(0, IN_INFO.find_index(TS.delta_2[0]));

    // 0509, add charge info
    charge.clear();
    charge.push_back(TS.chg[0]);

    size_t dep2_con = 0;
	bool first_non_H = true;
    for (size_t i = 0; i < TS.mu.size(); ++i){
      // new, 0430
		  if (TS.delta_2[2 * i + 2] == "H1") // check if the current atom is H
			 continue;
      if (TS.delta_2[2 * i + 3] == "1"){
        ++root_degree;
        if (height < 1) height = 1;
        while (!first_non_H && dep2_con < dmax - 1){
          seq.emplace_back(0, 0);
          charge.push_back(TS.chg[i + 1]);
          ++dep2_con;
        }
		first_non_H = false;
        seq.emplace_back(TS.mu[i], IN_INFO.find_index(TS.delta_2[2 * i + 2]));
        charge.push_back(TS.chg[i + 1]);
        dep2_con = 0;
      } else {
        if (height < 2) height = 2;
        seq.emplace_back(TS.mu[i], IN_INFO.find_index(TS.delta_2[2 * i + 2]));
        charge.push_back(TS.chg[i + 1]);
        ++dep2_con;
      }
    }

    while (dep2_con < dmax - 1){
      seq.emplace_back(0, 0);
      charge.push_back(0);
      ++dep2_con;
    }

    // calc_root_val();
  }

  FringeTree(
    const unsigned short& color,
    const size_t& dmax
  ){
    seq.clear();
    num_verts = 1;
    root_degree = 0;
    height = 0;

    seq.emplace_back(0, color);

    charge.clear();

    while (seq.size() < dmax + 1){
      seq.emplace_back(0, 0);
      charge.push_back(0);
    }
  }

  void calc_degree(const size_t& dmax, const size_t& root_degree = 2) {
    size_t m = seq.size();
    auto& prt = _prt[dmax];
    degree = vector<unsigned short>(m, 0);
    degree[0] = root_degree;
    for (size_t i = 1; i < m; ++i) {
      if (seq[i].first != 0) {
        ++degree[i];
        ++degree[prt[i]];
      }
    }
  }

  void calc_rv(const input_info& IN_INFO, const size_t& root_degree, const size_t& id) {
    rv.initialize_ec(IN_INFO.num_kind_atom_degree);
    rv.initialize_fc(IN_INFO.num_fc);
    ++rv.fc[IN_INFO.fc_map.at(id)];

    // size_t m = seq.size();
    // auto& prt = _prt[IN_INFO.dmax];
    // calc_degree(IN_INFO.dmax, root_degree);

    // for (size_t i = 1; i < m; ++i) {
    //   if (seq[i].first != 0) {
    //     size_t d1 = degree[i];
    //     size_t d2 = degree[prt[i]];  // then the direction should be d2->d1
    //     size_t mul = seq[i].first;
    //     atom_degree ad1(seq[i].second, d1);
    //     atom_degree ad2(seq[prt[i]].second, d2);
    //     ++rv.ec_ex[IN_INFO.atom_deg_map.at(ad2)][IN_INFO.atom_deg_map.at(ad1)]
    //               [mul];
    //   }
    // }
  }

  int calc_root_val(const size_t& dmax) {
    root_val = 0;
    size_t m = seq.size();
    auto& prt = _prt[dmax];
    for (size_t i = 1; i < m; ++i) {
      if (prt[i] == 0) {
        root_val += seq[i].first;
      }
    }
    // 0505, add charge info
    root_val -= charge[0];
    return root_val;
  }

  // void build_from_TS(
  //   const input_info& IN_INFO, 
  //   const TreeSignature& TS, 
  //   const size_t& dmax
  // ){
  //   seq.clear();
  //   num_verts = TS.mu.size() + 1;
  //   root_degree = 0;
  //   height = 0;

  //   seq.emplace_back(0, IN_INFO.find_index(TS.delta[0]));

  //   size_t dep2_con = 0;
  //   for (size_t i = 0; i < TS.mu.size(); ++i){
  //     if (TS.delta[2 * i + 3] == "1"){
  //       ++root_degree;
  //       if (height < 1) height = 1;
  //       while (i != 0 && dep2_con <= dmax - 1){
  //         seq.emplace_back(0, 0);
  //         ++dep2_con;
  //       }
  //       seq.emplace_back(TS.mu[i], IN_INFO.find_index(TS.delta[2 * i + 2]));
  //       dep2_con = 0;
  //     } else {
  //       if (height < 2) height = 2;
  //       seq.emplace_back(TS.mu[i], IN_INFO.find_index(TS.delta[2 * i + 2]));
  //       ++dep2_con;
  //     }
  //   }

  //   while (dep2_con < dmax - 1){
  //     seq.emplace_back(0, 0);
  //     ++dep2_con;
  //   }

  //   // calc_root_val();
  // }
};

class DAG_edge {
 public:
  DAG_node node;

  root_status rs;
  size_t w_ind;
  unsigned short k;

  DAG_edge() {}

  DAG_edge(
      const root_status& _rs1, const size_t& _l,
      const size_t& _ind,  // first three parameters show the node in the DAG
      const root_status& _rs, const size_t& _w_ind,
      const unsigned short& _mul  // last three parameters show the information
                                  // of the edge (fringe tree)
  ) {
    node = DAG_node(_rs1, _l, _ind);
    rs = _rs;
    w_ind = _w_ind;
    k = _mul;
  }
};

class DAG_edge_T {
 public:
  DAG_node_T node;

  root_status_height rsh;
  size_t T2_ind;
  unsigned short k;

  DAG_edge_T() {}

  DAG_edge_T(const root_status_T& _rsT, const size_t& _l, const size_t& _ind,
             const root_status_height& _rsh, const size_t& _T2_ind,
             const unsigned short& _mul) {
    node = DAG_node_T(_rsT, _l, _ind);
    rsh = _rsh;
    T2_ind = _T2_ind;
    k = _mul;
  }
};

class feasible_pair_W {
 public:
  DAG_node node1;
  DAG_node node2;

  unsigned short k;

  DAG_direction_status dir;

  feasible_pair_W() {}

  feasible_pair_W(const root_status& _rs1, const root_status& _rs2,
                  const size_t& _l1, const size_t& _ind1, const size_t& _l2,
                  const size_t& _ind2, const unsigned short& _k, const DAG_direction_status& _dir) {
    node1 = DAG_node(_rs1, _l1, _ind1);
    node2 = DAG_node(_rs2, _l2, _ind2);
    k = _k;
    dir = _dir;
  }
};

class feasible_pair_T {
 public:
  DAG_node_T node1;
  DAG_node_T node2;

  unsigned short k;

  DAG_direction_status dir;

  feasible_pair_T() {}

  feasible_pair_T(const root_status_T& _rs1, const root_status_T& _rs2,
                  const size_t& _l1, const size_t& _ind1, const size_t& _l2,
                  const size_t& _ind2, const unsigned short& _k, const DAG_direction_status& _dir) {
    node1 = DAG_node_T(_rs1, _l1, _ind1);
    node2 = DAG_node_T(_rs2, _l2, _ind2);
    k = _k;
    dir = _dir;
  }
};

// A data structure to store the graph generated by one component
class component_result {
 public:
  vector<map_tree> all_seq_h;
  vector<map_tree> all_seq;

  vector <map_tree> all_seq_h_loose;
  vector <map_tree> all_seq_loose;

  vector_2D<FringeTree> All_FT_W1;
  vector_2D<FringeTree> All_FT_W2;
  vector_2D<FringeTree> All_FT_W3;
  vector_2D<FringeTree> All_FT_W2_core;

  unordered_map<pair<root_status, size_t>, vector<size_t>> map_FT_W1;
  unordered_map<pair<root_status, size_t>, vector<size_t>> map_FT_W2;
  unordered_map<pair<root_status, size_t>, vector<size_t>> map_FT_W3;
  unordered_map<pair<root_status_height, size_t>, vector<size_t>> map_FT_W2_core;

  // unordered_map <root_status_height, map <size_t, pair <vector
  // <DAG_node_height>, size_t>>> map_T2;

  // unordered_map <DAG_node, vector <DAG_edge>> map_DAG_W1;
  // unordered_map <DAG_node, vector <DAG_edge>> map_DAG_W3;
  // unordered_map <DAG_node_T, vector <DAG_edge_T>> map_DAG_T;

  // vector <feasible_pair_W> all_fp_W;
  // vector <feasible_pair_T> all_fp_T;

  size_t num_h;
  size_t num;

  // size_t fp_id;
  // vector_2D <size_t> fp_W_path1;
  // vector_2D <size_t> fp_W_path2;

  // vector_3D <size_t> fp_T_path1;
  // vector_3D <size_t> fp_T_path2;

  bool tree_component;
  bool simple_tree;

  // root_status_height rsh_simple_tree;
  // size_t ind_simple_tree;

  component_result() {
    all_seq_h.clear();
    all_seq.clear();
    all_seq_h_loose.clear();
    all_seq_loose.clear();
    All_FT_W1.clear();
    All_FT_W2.clear();
    All_FT_W3.clear();
    All_FT_W2_core.clear();
    map_FT_W1.clear();
    map_FT_W2.clear();
    map_FT_W3.clear();
    map_FT_W2_core.clear();
    num_h = 0;
    num = 0;

    // map_T2.clear();

    // map_DAG_W.clear();
    // map_DAG_T.clear();
    // all_fp_W.clear();
    // fp_W_path1.clear();
    // fp_W_path2.clear();
    // fp_T_path1.clear();
    // fp_T_path2.clear();

    // fp_id = 0;
  }

  component_result(size_t M) {
    all_seq_h.clear();
    all_seq.clear();
    all_seq_h_loose.clear();
    all_seq_loose.clear();
    All_FT_W1.clear();
    All_FT_W2.clear();
    All_FT_W3.clear();
    All_FT_W2_core.clear();
    map_FT_W1.clear();
    map_FT_W2.clear();
    map_FT_W3.clear();
    map_FT_W2_core.clear();

    All_FT_W1.resize(M);
    All_FT_W2.resize(M);
    All_FT_W3.resize(M);
    All_FT_W2_core.resize(M);

    num_h = 0;
    num = 0;

    // map_T2.clear();

    // map_DAG_W.clear();
    // map_DAG_T.clear();
    // all_fp_W.clear();
    // fp_W_path1.clear();
    // fp_W_path2.clear();
    // fp_T_path1.clear();
    // fp_T_path2.clear();

    // fp_id = 0;
  }

  void initialize(size_t M) {
    all_seq_h.clear();
    all_seq.clear();
    all_seq_h_loose.clear();
    all_seq_loose.clear();
    All_FT_W1.clear();
    All_FT_W2.clear();
    All_FT_W3.clear();
    All_FT_W2_core.clear();
    map_FT_W1.clear();
    map_FT_W2.clear();
    map_FT_W3.clear();
    map_FT_W2_core.clear();

    All_FT_W1.resize(M);
    All_FT_W2.resize(M);
    All_FT_W3.resize(M);
    All_FT_W2_core.resize(M);

    num_h = 0;
    num = 0;

    // map_T2.clear();

    // map_DAG_W.clear();
    // map_DAG_T.clear();
    // all_fp_W.clear();
    // fp_W_path1.clear();
    // fp_W_path2.clear();
    // fp_T_path1.clear();
    // fp_T_path2.clear();

    // fp_id = 0;
  }

  // void initialize_W_path(
  //     const DAG_node& _node,
  //     vector_2D <size_t>& fp_path,
  //     unordered_map <DAG_node, vector <DAG_edge>>& map_DAG
  // ){  // W3
  //     auto p = _node;
  //     while (map_DAG.at(p).size() > 0){
  //         vector <size_t> tmp = {0, 0};
  //         fp_path.push_back(tmp);
  //         p = map_DAG[p][0].node;
  //     }

  //     fp_path.push_back(vector <size_t> ({0, 0}));
  // }

  // void initialize_W(){
  //     fp_W_path1.clear();
  //     fp_W_path2.clear();

  //     initialize_W_path(all_fp_W[fp_id].node1, fp_W_path1, map_DAG_W3);
  //     initialize_W_path(all_fp_W[fp_id].node2, fp_W_path2, map_DAG_W1);
  // }

  // bool get_next_graph_W_path1(map_rv& seq, const DAG_node& _node, const bool
  // change){   // W3
  //     map_rv tmp;
  //     auto p = _node;
  //     int i = 0;
  //     vector <DAG_node> p_vector;
  //     while (map_DAG_W3.at(p).size() > 0){
  //         p_vector.push_back(p);
  //         auto& _tmp = map_DAG_W3[p][fp_W_path1[i][0]];

  //         size_t w_ind = _tmp.w_ind;
  //         FT_ind = map_FT_W2[make_pair(_tmp.rs, w_ind)][fp_W_path1[i][1]];

  //         tmp.rs.push_back(_tmp.rs);
  //         tmp.w_ind.push_back(FT_ind);
  //         tmp.mul.push_back(_tmp.k);

  //         p = map_DAG_W3[p][fp_W_path1[i][0]].node;
  //         ++i;
  //     }
  //     p_vector.push_back(p);
  //     tmp.rs.push_back(p.rs.color);
  //     tmp.w_ind.public(map_FT_W3[make_pair(p.rs.color,
  //     p.ind)][fp_W_path1[i][1]]);

  //     seq.extend(tmp, false);

  //     if (change){
  //         while (i >= 0){
  //             ++fp_W_path1[i][1];
  //             auto& p = p_vector[i];
  //             if (i == p_vector.size() - 1){
  //                 if (fp_W_path1[i][1] >= map_FT_W3[make_pair(p.rs.color,
  //                 p.ind)].size()){
  //                     --i;
  //                     fp_W_path1.pop_back();
  //                 } else {
  //                     break;
  //                 }
  //             } else {
  //                 if (fp_W_path1[i][1] >= map_FT_W2[make_pair(p.rs.color,
  //                 p.ind)].size()){
  //                     ++fp_W_path1[i][0];
  //                     fp_W_path1[i][1] = 0;
  //                     if (fp_W_path1[i][0] >= map_DAG_W3[p].size()){
  //                         --i;
  //                         fp_W_path1.pop_back();
  //                     } else {
  //                         initialize_W_path(map_DAG_W3[p][fp_W_path1[i][0]].node,
  //                         fp_W_path1, map_DAG_W3); break;
  //                     }
  //                 } else {
  //                     initialize_W_path(map_DAG_W3[p][fp_W_path1[i][0]].node,
  //                     fp_W_path1, map_DAG_W3); break;
  //                 }
  //             }
  //         }
  //     }

  //     if (i < 0){
  //         return false;
  //     } else {
  //         return true;
  //     }

  // }

  // bool get_next_graph_W_path2(map_rv& seq, const DAG_node& _node, const bool
  // change){   // W1
  //     map_rv tmp;
  //     auto p = _node;
  //     int i = 0;
  //     vector <DAG_node> p_vector;
  //     while (map_DAG_W1.at(p).size() > 0){
  //         p_vector.push_back(p);
  //         auto& _tmp = map_DAG_W1[p][fp_W_path2[i][0]];

  //         size_t w_ind = _tmp.w_ind;
  //         FT_ind = map_FT_W2[make_pair(_tmp.rs, w_ind)][fp_W_path2[i][1]];

  //         tmp.rs.push_back(_tmp.rs);
  //         tmp.w_ind.push_back(FT_ind);
  //         tmp.mul.push_back(_tmp.k);

  //         p = map_DAG_W1[p][fp_W_path2[i][0]].node;
  //         ++i;
  //     }
  //     p_vector.push_back(p);
  //     tmp.rs.push_back(p.rs.color);
  //     tmp.w_ind.public(map_FT_W1[make_pair(p.rs.color,
  //     p.ind)][fp_W_path2[i][1]]);

  //     seq.extend(tmp, true);

  //     if (change){
  //         while (i >= 0){
  //             ++fp_W_path2[i][1];
  //             auto& p = p_vector[i];
  //             if (i == p_vector.size() - 1){
  //                 if (fp_W_path2[i][1] >= map_FT_W1[make_pair(p.rs.color,
  //                 p.ind)].size()){
  //                     --i;
  //                     fp_W_path2.pop_back();
  //                 } else {
  //                     break;
  //                 }
  //             } else {
  //                 if (fp_W_path2[i][1] >= map_FT_W2[make_pair(p.rs.color,
  //                 p.ind)].size()){
  //                     ++fp_W_path2[i][0];
  //                     fp_W_path2[i][1] = 0;
  //                     if (fp_W_path2[i][0] >= map_DAG_W1[p].size()){
  //                         --i;
  //                         fp_W_path2.pop_back();
  //                     } else {
  //                         initialize_W_path(map_DAG_W1[p][fp_W_path2[i][0]].node,
  //                         fp_W_path2, map_DAG_W1); break;
  //                     }
  //                 } else {
  //                     initialize_W_path(map_DAG_W1[p][fp_W_path2[i][0]].node,
  //                     fp_W_path2, map_DAG_W1); break;
  //                 }
  //             }
  //         }
  //     }

  //     if (i < 0){
  //         return false;
  //     } else {
  //         return true;
  //     }

  // }

  // bool get_next_graph(vector <map_tree>& seq, const bool change, const size_t
  // rho = 2){
  //     if (tree_component && simple_tree){
  //         auto& indices =
  //         map_FT_W2_core.at(rsh_simple_tree).at(ind_simple_tree);

  //         map_rv tmp(rsh_simple_tree.color, indices[fp_id]);
  //         seq = map_tree(tmp, rsh_simple_tree.height);

  //         if (change){
  //             ++fp_id;
  //             if (fp_id == indices.size()){
  //                 fp_id = 0;
  //                 return false;
  //             } else {
  //                 return true;
  //             }
  //         } else {
  //             return true;
  //         }

  //     } else if (tree_component){
  //         map_rv tmp2;
  //         tmp2.rs.clear();
  //         tmp2.w_ind.clear();
  //         tmp2.mul.clear();

  //         unsigned short k = all_fp_W[fp_id].k;
  //         size_t h = all_fp_W[fp_id].node1.l + all_fp_W[fp_id].node2.l + 1 +
  //         rho;

  //         bool flag = false;
  //         flag = flag || get_next_graph_W_path2(tmp2, all_fp_W[fp_id].node2,
  //         change); if (!flag){
  //             fp_W_path2.clear();
  //             initialize_W_path(all_fp_W[fp_id].node2, fp_W_path2,
  //             map_DAG_W1);
  //         }

  //         map_rv tmp1;
  //         tmp1.rs.clear();
  //         tmp1.w_ind.clear();
  //         tmp1.mul.clear();

  //         flag = flag || get_next_graph_W_path1(tmp1, all_fp_W[fp_id].node1,
  //         !flag); if (!flag){
  //             ++fp_id;
  //             if (fp_id == all_fp_W.size()){
  //                 fp_id = 0;
  //             } else {
  //                 flag = true;
  //             }
  //             initialize_W();
  //         }

  //         map_rv tmp(tmp1, tmp2, k, false);

  //         seq = map_tree(tmp, h);

  //         return flag;
  //     } else { // edge component

  //     }
  // }
};

// A trie structure to store the weight vector
class PathTrie {
 public:
  unsigned short key;
  int parent_ind;  // parent index
  vector<int> children;

  PathTrie() {}

  PathTrie(unsigned short _key) : key(_key) {
    parent_ind = -1;
    children.clear();
  }

  PathTrie(unsigned short _key, const size_t& _parent) : key(_key) {
    parent_ind = _parent;
    children.clear();
  }
};

size_t find_trie_index(const resource_vector_1D& rv,
                       vector<PathTrie>& all_trie_node) {
  size_t depth = 0;
  size_t id = 0;

  while (depth < rv.size()) {
    bool flag = false;
    int next_id = -1;
    for (auto& _id : all_trie_node[id].children) {
      if (all_trie_node[_id].key == rv[depth]) {
        next_id = _id;
        flag = true;
        break;
      }
    }

    if (!flag) {
      all_trie_node.emplace_back(rv[depth], id);
      next_id = all_trie_node.size() - 1;
      all_trie_node[id].children.push_back(next_id);
    }

    id = next_id;
    ++depth;
  }

  return id;
}

resource_vector_1D traverse(size_t leaf_id,
                            const vector<PathTrie>& all_trie_node,
                            size_t _size) {
  resource_vector_1D tmp(_size);
  size_t id = leaf_id;
  size_t tmp_ind = _size - 1;
  while (id != 0) {
    tmp[tmp_ind] = all_trie_node[id].key;
    --tmp_ind;
    id = all_trie_node[id].parent_ind;
  }

  return tmp;
}

class component {
 public:
  input_info IN_INFO;
  component_result comp_graph;
  bool no_feasible_pair;
  bool no_feasible_pair_loose;
  bool with_core_height;
};

bool cmp_v_component(const component& a, const component& b) {
  if (a.IN_INFO.n != b.IN_INFO.n) {
    return (a.IN_INFO.n < b.IN_INFO.n);
  } else if (a.IN_INFO.height != b.IN_INFO.height) {
    return (a.IN_INFO.height < b.IN_INFO.height);
  } else if (a.IN_INFO.Gamma_int.size() != b.IN_INFO.Gamma_int.size()) {
    return (a.IN_INFO.Gamma_int.size() < b.IN_INFO.Gamma_int.size());
  } else {
    return false;
  }
}

bool cmp_e_component(const component& a, const component& b) {
  if (a.IN_INFO.n != b.IN_INFO.n) {
    return (a.IN_INFO.n < b.IN_INFO.n);
  } else if (a.IN_INFO.height != b.IN_INFO.height) {
    return (a.IN_INFO.height < b.IN_INFO.height);
  } else if (a.IN_INFO.base_vertices.size() != b.IN_INFO.base_vertices.size()) {
    return (a.IN_INFO.base_vertices.size() < b.IN_INFO.base_vertices.size());
  } else if (a.IN_INFO.Gamma_int.size() != b.IN_INFO.Gamma_int.size()) {
    return (a.IN_INFO.Gamma_int.size() < b.IN_INFO.Gamma_int.size());
  } else {
    return false;
  }
}

#endif
