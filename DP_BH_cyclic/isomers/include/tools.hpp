// This file contains some minor functions used for the code.

#ifndef INCLUDES_TOOLS_HPP_
#define INCLUDES_TOOLS_HPP_

#include <cstdlib>
#include <iostream>
#include <vector>

#include "data_structures.hpp"

using namespace std;

void print(const FringeTree& FT, const input_info& IN_INFO) {
  int count = 0;
  for (size_t i = 0; i < FT.seq.size(); ++i) {
    if (count != 0 && FT.seq[i].first == 0) {
      cout << "(0, 0), ";
    } else {
      cout << "(" << FT.seq[i].first << ", "
           << IN_INFO.atoms[FT.seq[i].second].NAME << "), ";
    }
    count += 1;
  }
  cout << endl;
}

string get_fringe_tree_string(const FringeTree& FT) {
  string ans = "";
  for (size_t i = 0; i < FT.seq.size(); ++i) {
    ans += to_string(FT.seq[i].first) + to_string(FT.seq[i].second);
  }
  return ans;
}

void printseq(const input_info& IN_INFO, vector<MultCol>& seq) {
  for (size_t i = 0; i < seq.size(); ++i) {
    cout << seq[i].first << IN_INFO.atoms[seq[i].second].NAME;
  }
  cout << endl;
}

// bool dc_compare_out(const input_info& IN_INFO,
//     const vector <unsigned short>& _bc_num){
//     // for (auto& dc_t : Dcset){
//     //     if (_dc_map.at(dc_t) > IN_INFO.dc_map.at(dc_t))
//     //         return false;
//     // }
//     for (size_t i = 0; i < Dcset.size(); ++i){
//     	if (_bc_num[i] > IN_INFO.rv.bc_num_out[i]){
//     		return false;
//     	}
//     }
//     return true;
// }

// new functions//
// computing sum or difference of two 3D vectors without adding entries
template <class T>
vector_3D<T> SUMinus3D(
    const vector_3D<T>& A, const vector_3D<T>& B,
    const bool& decision)  // decision true = compute sum; else difference
{
  size_t Lambdasize = A.size();
  vector_3D<T> Sum3Dvectors(Lambdasize);  // output
  for (size_t i = 0; i < Lambdasize; ++i) {
    Sum3Dvectors[i] = vector_2D<T>(Lambdasize);
    for (size_t j = 0; j < Lambdasize; ++j) {
      Sum3Dvectors[i][j] = vector<T>(4, 0);
    }
  }
  if (decision) {  // compute sum
    for (size_t i = 0; i < Lambdasize; ++i) {
      for (size_t j = 0; j < Lambdasize; ++j) {
        for (size_t k = 1; k <= 3; ++k) {
          Sum3Dvectors[i][j][k] = A[i][j][k] + B[i][j][k];
        }  // for k
      }    // for j
    }      // for i
  }        // if
  else {   // compute difference
    for (size_t i = 0; i < Lambdasize; ++i) {
      for (size_t j = 0; j < Lambdasize; ++j) {
        for (size_t k = 1; k <= 3; ++k) {
          Sum3Dvectors[i][j][k] = A[i][j][k] - B[i][j][k];
        }  // for k
      }    // for j
    }      // for i
  }        // else
  return Sum3Dvectors;
}

resource_vector_1D SUMinus3D(const resource_vector_1D& A,
                             const resource_vector_1D& B, const bool decision) {
  size_t Lambdasize = A.size();
  resource_vector_1D ans(Lambdasize);

  if (decision) {
    for (size_t i = 0; i < Lambdasize; ++i) {
      ans[i] = A[i] + B[i];
    }
  } else {
    for (size_t i = 0; i < Lambdasize; ++i) {
      ans[i] = A[i] - B[i];
    }
  }

  return ans;
}

resource_vector SUMinus3D(const resource_vector& A, const resource_vector& B,
                          const bool decision) {
  size_t Lambdasize = A.ec_co.size();
  resource_vector ans(Lambdasize);

  if (decision) {
    for (size_t i = 0; i < Lambdasize; ++i) {
      for (size_t j = 0; j < Lambdasize; ++j) {
        for (size_t k = 1; k <= 3; ++k) {
          ans.ec_co[i][j][k] = A.ec_co[i][j][k] + B.ec_co[i][j][k];
          ans.ec_in[i][j][k] = A.ec_in[i][j][k] + B.ec_in[i][j][k];
          ans.ec_ex[i][j][k] = A.ec_ex[i][j][k] + B.ec_ex[i][j][k];
        }  // for k
      }    // for j
    }      // for i
    ans.bl = A.bl + B.bl;
  } else {
    for (size_t i = 0; i < Lambdasize; ++i) {
      for (size_t j = 0; j < Lambdasize; ++j) {
        for (size_t k = 1; k <= 3; ++k) {
          ans.ec_co[i][j][k] = A.ec_co[i][j][k] - B.ec_co[i][j][k];
          ans.ec_in[i][j][k] = A.ec_in[i][j][k] - B.ec_in[i][j][k];
          ans.ec_ex[i][j][k] = A.ec_ex[i][j][k] - B.ec_ex[i][j][k];
        }  // for k
      }    // for j
    }      // for i
    ans.bl = A.bl - B.bl;
  }
  return ans;
}

/*
 * Test if a 3D vector A <= B
 * true if A <= B
 */
template <class T>
bool TestLessEqAB(const vector_3D<T>& A, const vector_3D<T>& B) {
  size_t Lambdasize = A.size();
  for (size_t i = 1; i < Lambdasize; ++i) {
    for (size_t j = 1; j < Lambdasize; ++j) {
      for (size_t k = 1; k <= 3; ++k) {
        if (A[i][j][k] > B[i][j][k]) return false;
      }
    }
  }
  return true;
}

bool TestLessEqAB(const resource_vector_1D& A, const resource_vector_1D& B) {
  size_t Lambdasize = A.size();
  for (size_t i = 1; i < Lambdasize; ++i) {
    if (A[i] > B[i]) return false;
  }
  return true;
}

bool TestEqAB(const resource_vector_1D& A, const resource_vector_1D& B) {
  size_t Lambdasize = A.size();
  for (size_t i = 1; i < Lambdasize; ++i) {
    if (A[i] != B[i]) return false;
  }
  return true;
}

bool TestLessEqAB(const resource_vector& A, const resource_vector& B) {
  size_t Lambdasize = A.ec_co.size();
  for (size_t i = 0; i < Lambdasize; ++i) {
    for (size_t j = 0; j < Lambdasize; ++j) {
      for (size_t k = 1; k <= 3; ++k) {
        if (A.ec_co[i][j][k] > B.ec_co[i][j][k]) return false;
        if (A.ec_in[i][j][k] > B.ec_in[i][j][k]) return false;
        if (A.ec_ex[i][j][k] > B.ec_ex[i][j][k]) return false;
      }
    }
  }
  if (A.bl > B.bl) return false;
  return true;
}

void print_rv(const input_info& IN_INFO) {
  size_t M = IN_INFO.num_kind_atom_degree;
  cout << "ec_co : " << endl;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      for (size_t k = 1; k <= 3; ++k) {
        auto& a = IN_INFO.atom_deg[i];
        auto& b = IN_INFO.atom_deg[j];
        cout << "#(" << IN_INFO.atoms[a.color].NAME << a.deg << ", "
             << IN_INFO.atoms[b.color].NAME << b.deg;
        cout << ", " << k << ") = " << IN_INFO.rv.ec_co[i][j][k] << endl;
      }
    }
  }
  cout << "ec_in : " << endl;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < M; ++j) {
      for (size_t k = 1; k <= 3; ++k) {
        auto& a = IN_INFO.atom_deg[i];
        auto& b = IN_INFO.atom_deg[j];
        cout << "#(" << IN_INFO.atoms[a.color].NAME << a.deg << ", "
             << IN_INFO.atoms[b.color].NAME << b.deg;
        cout << ", " << k << ") = " << IN_INFO.rv.ec_in[i][j][k] << endl;
      }
    }
  }
  cout << "ec_ex : " << endl;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < M; ++j) {
      for (size_t k = 1; k <= 3; ++k) {
        auto& a = IN_INFO.atom_deg[i];
        auto& b = IN_INFO.atom_deg[j];
        cout << "#(" << IN_INFO.atoms[a.color].NAME << a.deg << ", "
             << IN_INFO.atoms[b.color].NAME << b.deg;
        cout << ", " << k << ") = " << IN_INFO.rv.ec_ex[i][j][k] << endl;
      }
    }
  }
  cout << "bl = " << IN_INFO.rv.bl << endl;
}

void print_rv1D(const input_info& IN_INFO, const resource_vector_1D& rv1D) {
  cout << "bl = " << rv1D[0] << endl;
  for (size_t i = 1; i < IN_INFO.rv_edges.size(); ++i) {
    if (rv1D[i] == 0) continue;
    rv_edge tmp = IN_INFO.rv_edges[i];
    cout << "#(" << tmp.type << " " << IN_INFO.atoms[tmp.ad1.color].NAME
         << tmp.ad1.deg << " ";
    cout << IN_INFO.atoms[tmp.ad2.color].NAME << tmp.ad2.deg << " " << tmp.mul
         << ") = " << rv1D[i] << endl;
  }
}

void print_rv1D_2(const input_info& IN_INFO, const resource_vector_1D& rv1D) {
  // cout << "bl = " << rv1D[0] << endl;
  for (size_t i = 1; i < IN_INFO.rv_edges.size(); ++i) {
    if (rv1D[i] == 0) continue;
    rv_edge tmp = IN_INFO.rv_edges[i];
    // cout << "#(" << tmp.type << " " << IN_INFO.atoms[tmp.ad1.color].NAME
    //      << tmp.ad1.deg << " ";
    // cout << IN_INFO.atoms[tmp.ad2.color].NAME << tmp.ad2.deg << " " <<
    // tmp.mul
    //      << ") = " << rv1D[i] << endl;
    string at1 = IN_INFO.atoms[tmp.ad1.color].NAME;
    string at2 = IN_INFO.atoms[tmp.ad2.color].NAME;
    at1[0] = tolower(at1[0]);
    at2[0] = tolower(at2[0]);
    string type;
    if (tmp.type == "in") {
      type = "inn";
    } else {
      type = tmp.type;
    }

    cout << "\\" << rv1D[i] << "_{(\\at" << at1 << tmp.ad1.deg << ", \\at"
         << at2 << tmp.ad2.deg << ", " << tmp.mul << ")}^{\\" << type << "} + ";
  }

  cout << "$,\\\\" << endl;
}

void sort_child(vector<int>& e_ind, const vector<DAG_edge>& child,
                const unordered_map<DAG_node, DAG_count>& DAG_count) {
  vector<pair<size_t, int>> tmp;
  int ind = 0;
  for (auto& e : child) {
    auto& node = e.node;
    auto& path_num = DAG_count.at(node).path;
    tmp.push_back({path_num, ind});
    ind++;
  }
  sort(tmp.begin(), tmp.end());

  for (auto& e : tmp) {
    auto& ind = e.second;
    e_ind.push_back(ind);
  }
  return;
}

void sort_child_T(vector<int>& e_ind, const vector<DAG_edge_T>& child,
                  const unordered_map<DAG_node_T, DAG_count_T>& DAG_count) {
  vector<pair<size_t, int>> tmp;
  int ind = 0;
  for (auto& e : child) {
    auto& node = e.node;
    auto& path_num = DAG_count.at(node).path;
    tmp.push_back({path_num, ind});
    ind++;
  }
  sort(tmp.begin(), tmp.end());

  for (auto& e : tmp) {
    auto& ind = e.second;
    e_ind.push_back(ind);
  }
  return;
}

// A function to free memory.
template <class T>
void destroy(T& a) {
  T().swap(a);
}

#endif /* INCLUDES_TOOLS_HPP_ */
