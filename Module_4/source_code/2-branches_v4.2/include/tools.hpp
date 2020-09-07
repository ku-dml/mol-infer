// This file contains some minor functions used for the code.

#ifndef INCLUDES_TOOLS_HPP_
#define INCLUDES_TOOLS_HPP_

#include <iostream>
#include <cstdlib>
#include <vector>

#include "data_structures.hpp"

using namespace std;

void print(const FringeTree& FT) {
	for (size_t i = 0; i < FT.seq.size(); ++i) {
		cout << "(" << FT.seq[i].first << ", " << FT.seq[i].second << "), ";
	}
	cout << endl;
}

string get_fringe_tree_string(const FringeTree& FT){
	string ans = "";
	for (size_t i = 0; i < FT.seq.size(); ++i) {
		ans += to_string(FT.seq[i].first) + to_string(FT.seq[i].second);
	}
	return ans;
}

vector <MultCol> reverse_seq(vector <MultCol>& seq){
	vector <MultCol> P(seq.size());
	for (size_t i = 0; i < seq.size(); ++i){
		if (i == 0){
			P[i] = make_pair(0, seq[seq.size() - 1].second);
		} else {
			P[i] = make_pair(seq[seq.size() - i].first, seq[seq.size() - 1 - i].second);
		}
	}
	return P;
}

void printseq(const input_info& IN_INFO, vector <MultCol>& seq){
    for (size_t i = 0; i < seq.size(); ++i){
        cout << seq[i].first << IN_INFO.atoms[seq[i].second].NAME;
    }
    cout << endl;
}

bool dc_compare_out(const input_info& IN_INFO,
    const vector <unsigned short>& _bc_num){
    // for (auto& dc_t : Dcset){
    //     if (_dc_map.at(dc_t) > IN_INFO.dc_map.at(dc_t))
    //         return false;
    // }
    for (size_t i = 0; i < Dcset.size(); ++i){
    	if (_bc_num[i] > IN_INFO.rv.bc_num_out[i]){
    		return false;
    	}
    }
    return true;
}

//new functions//
//computing sum or difference of two 3D vectors without adding entries
template <class T>
vector_3D <T> SUMinus3D(const vector_3D<T>& A,
	const vector_3D<T>& B,
	const bool &decision)// decision true = compute sum; else difference
{
	size_t Lambdasize = A.size();
	vector_3D <T> Sum3Dvectors(Lambdasize); // output
	for (size_t i = 0; i < Lambdasize; ++i) {
		Sum3Dvectors[i] = vector_2D<T>(Lambdasize);
		for (size_t j = 0; j < Lambdasize; ++j) {
			Sum3Dvectors[i][j] = vector<T>(4, 0);
		}
	}
	if (decision) { // compute sum
		for (size_t i = 0; i < Lambdasize; ++i) {
			for (size_t j = 0; j < Lambdasize; ++j) {
				for (size_t k = 0; k <= 3; ++k) {
					Sum3Dvectors[i][j][k] = A[i][j][k] + B[i][j][k];
				} // for k
			} // for j
		} // for i
	}//if 
	else {// compute difference
		for (size_t i = 0; i < Lambdasize; ++i) {
			for (size_t j = 0; j < Lambdasize; ++j) {
				for (size_t k = 0; k <= 3; ++k) {
					Sum3Dvectors[i][j][k] = A[i][j][k] - B[i][j][k];
				} // for k
			} // for j
		} // for i
	}// else
	return Sum3Dvectors;
}

resource_vector_1D SUMinus3D(const resource_vector_1D& A, const resource_vector_1D& B, const bool decision){
	size_t Lambdasize = A.size();
	resource_vector_1D ans(Lambdasize);

	if (decision){
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

resource_vector SUMinus3D(const resource_vector& A, const resource_vector& B, const bool decision){
	size_t Lambdasize = A.CVE_in.size();
	resource_vector ans(Lambdasize - 1);

	if (decision){
		for (size_t i = 0; i < Lambdasize; ++i) {
			for (size_t j = 0; j < Lambdasize; ++j) {
				for (size_t k = 0; k <= 3; ++k) {
					ans.CVE_in[i][j][k] = A.CVE_in[i][j][k] + B.CVE_in[i][j][k];
					ans.CVE_out[i][j][k] = A.CVE_out[i][j][k] + B.CVE_out[i][j][k];
				} // for k
			} // for j
		} // for i

		for (size_t i = 0; i < Dcset.size(); ++i){
			ans.bc_num_in[i] = A.bc_num_in[i] + B.bc_num_in[i];
			ans.bc_num_out[i] = A.bc_num_out[i] + B.bc_num_out[i];
		}

		for (size_t i = 1; i <= 4; ++i){
			ans.deg_in[i] = A.deg_in[i] + B.deg_in[i];
			ans.deg_out[i] = A.deg_out[i] + B.deg_out[i];
		}
	} else {
		for (size_t i = 0; i < Lambdasize; ++i) {
			for (size_t j = 0; j < Lambdasize; ++j) {
				for (size_t k = 0; k <= 3; ++k) {
					ans.CVE_in[i][j][k] = A.CVE_in[i][j][k] - B.CVE_in[i][j][k];
					ans.CVE_out[i][j][k] = A.CVE_out[i][j][k] - B.CVE_out[i][j][k];
				} // for k
			} // for j
		} // for i

		for (size_t i = 0; i < Dcset.size(); ++i){
			ans.bc_num_in[i] = A.bc_num_in[i] - B.bc_num_in[i];
			ans.bc_num_out[i] = A.bc_num_out[i] - B.bc_num_out[i];
		}

		for (size_t i = 1; i <= 4; ++i){
			ans.deg_in[i] = A.deg_in[i] - B.deg_in[i];
			ans.deg_out[i] = A.deg_out[i] - B.deg_out[i];
		}
	}

	return ans;

}

/*
* Test if a 3D vector A <= B
* true if A <= B
*/
template <class T>
bool TestLessEqAB(const vector_3D<T>& A,
	const vector_3D<T>& B)
{ //
	size_t Lambdasize = A.size();
	for (size_t i = 1; i < Lambdasize; ++i) {
		for (size_t j = 1; j < Lambdasize; ++j) {
			for (size_t k = 0; k <= 3; ++k) {
				if (A[i][j][k] > B[i][j][k]) return false;
			}
		}
	}
	return true;
}

bool TestLessEqAB(const resource_vector_1D& A, const resource_vector_1D& B){
	size_t Lambdasize = A.size();
	for (size_t i = 1; i < Lambdasize; ++i){
		if (A[i] > B[i]) return false;
	}
	return true;
}

bool TestLessEqAB(const resource_vector& A, const resource_vector& B){
	size_t Lambdasize = A.CVE_in.size();
	for (size_t i = 1; i < Lambdasize; ++i) {
		for (size_t j = 1; j < Lambdasize; ++j) {
			for (size_t k = 0; k <= 3; ++k) {
				// cout << "M = "  << Lambdasize << "  i = " << i << " j = " << j << " k = " << k << endl;
				if (A.CVE_in[i][j][k] > B.CVE_in[i][j][k]) return false;
			}
		}
	}
	for (size_t i = 1; i < Lambdasize; ++i) {
		for (size_t j = 1; j < Lambdasize; ++j) {
			for (size_t k = 0; k <= 3; ++k) {
				// cout << "M = "  << Lambdasize << "  i = " << i << " j = " << j << " k = " << k << endl;
				if (A.CVE_out[i][j][k] > B.CVE_out[i][j][k]) return false;
			}
		}
	}
	for (size_t i = 0; i < Dcset.size(); ++i){
		if (A.bc_num_in[i] > B.bc_num_in[i]) return false;
	}
	for (size_t i = 0; i < Dcset.size(); ++i){
		if (A.bc_num_out[i] > B.bc_num_out[i]) return false;
	}
	for (size_t i = 1; i <= 4; ++i){
		if (A.deg_in[i] > B.deg_in[i]) return false;
	}
	for (size_t i = 1; i <= 4; ++i){
		if (A.deg_out[i] > B.deg_out[i]) return false;
	}
	return true;
}

void print(const resource_vector_1D& rv1D, size_t M){

	size_t ind = 0;
    for (size_t i = 1; i <= M; ++i){
        for (size_t j = 1; j <= i; ++j){
            for (size_t k = 0; k <= 3; ++k){
            	cout << rv1D[ind] << " ";
            	++ind;
            }
            cout << "     ";
        }
        cout << endl;
    }
    for (size_t i = 1; i <= M; ++i){
        for (size_t j = 1; j <= i; ++j){
            for (size_t k = 0; k <= 3; ++k){
            	cout << rv1D[ind] << " ";
            	++ind;
            }
            cout << "     ";
        }
        cout << endl;
    }
    for (size_t i = 0; i < Dcset.size(); ++i){
        cout << rv1D[ind] << " ";
        ++ind;
    }
    cout << endl;
    for (size_t i = 0; i < Dcset.size(); ++i){
        cout << rv1D[ind] << " ";
        ++ind;
    }
    cout << endl;
    for (size_t i = 1; i <= 4; ++i){
        cout << rv1D[ind] << " ";
        ++ind;
    }
    cout << endl;
    for (size_t i = 1; i <= 4; ++i){
        cout << rv1D[ind] << " ";
        ++ind;
    }
    cout << endl << endl;;
}

// A function used to free memory.
template <class T>
void destroy(T& a){
	T().swap(a);
}

#endif /* INCLUDES_TOOLS_HPP_ */
