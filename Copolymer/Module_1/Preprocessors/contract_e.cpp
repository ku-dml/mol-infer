/************************************************************
  A code to contract the e* and output link edge at the last.

 ************************************************************/

/********** include header files **********/
#include "fv_def.h"
#include "fv_common.hpp"
#include "compute_fc.hpp"

#include "./fringe/ChemicalGraph.hpp"

#define N_INDEX 0



/********** read arguments **********/
void read_arguments(int argc, char *argv[], string& input_file, string& output_file){
  if(argc!=3){
    fprintf(stderr, "usage: %s (input.sdf)(output.sdf)\n\n", argv[0]);
    fprintf(stderr, "Caution: in this model, input.sdf MUST BE completely H-defined and with one e* atom.\n\n");
    fprintf(stderr, "The program outputs:\n");
    fprintf(stderr, "\toutput.sdf: sdf file after processing\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  input_file = argv[1];
  output_file = argv[2];
    
}


/********** main **********/
int main(int argc, char *argv[]){
  
  string input_file, output_file;

  read_arguments(argc, argv, input_file, output_file);

  ifstream infile;

  try {
    infile.open(input_file);
  } catch (const exception& e){
    cerr << "Cannot open '" << input_file << "' for reading!" << endl;
    throw(e);
  }
  if (!infile.is_open()){
    cerr << "Error, infile not initialized!" << endl;
    throw(-1);
  }

  ofstream output_tmp(output_file, ios::out);
  output_tmp.close();

  // read chemical compounds one by one
  while (! infile.eof()){
    ChemicalGraph _g; // must with H
    int cond = _g.readFromFile(infile, false, true);
    if (cond == 0){
      break;
    }
    if (cond == -2){
      continue;
    }

    vector <Vertex> e_id_vector;
    e_id_vector.clear();

    Vertex e_id;

    for (Vertex u = 0; u < _g.numAtom; ++u){
      if (_g.alpha[u] == "e*"){
        e_id = u;
        e_id_vector.push_back(u);
      }
    }

    if (e_id_vector.size() == 1){
      
      if (_g.adj[e_id].size() != 2){
        cerr << "Error, the degree of e* is not 2 !!! CID : " << _g.CID << endl;
        throw(-2);
      }

      /////////////////////////////////////////////////////////////////////
      // ignore the graph that will create multiple edge
      // size_t e_id_1 = _g.adj[e_id][0];
      // size_t e_id_2 = _g.adj[e_id][1];
      // if (find(_g.adj[e_id_1].begin(), _g.adj[e_id_1].end(), e_id_2) != _g.adj[e_id_1].end()) continue;
      // ignore the graph that will create self-loop
      size_t e_id_1 = _g.adj[e_id][0];
      size_t e_id_2 = _g.adj[e_id][1];
      if (e_id_1 == e_id_2) {
        cout << "CID : "<< _g.CID << " will create self-loop after contraction!! ignored." << endl;
        continue;
      }
      /////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////
      // output SDF file
      ofstream output(output_file, ios::app);

      size_t real_n = _g.numAtom - 1;
      size_t real_m = 0;

      for (Vertex u = 0; u < _g.numAtom; ++u){
        for (Vertex v = u + 1; v < _g.numAtom; ++v){
          if (_g.beta[u][v] > 0){
            ++real_m;
          }
        }
      }
      --real_m; // erase the edge incident to e*

      output << _g.CID << "\n";
      output << ""
             << "\n";
      output << ""
             << "\n";
      output << std::setw(3) << real_n << std::setw(3) << real_m
             << "  0  0  0  0  0  0  0  0999 V2000 "
             << "\n";

      for (size_t i = 0; i < _g.numAtom; ++i) {
        string atom_symbol = _g.alpha[i];
        if (atom_symbol == "e*") continue; // ignore e*
        int chg_ele = 0;
        if (_g.charge[i] != 0) chg_ele = 4 - _g.charge[i];
        output << "    0.0000    0.0000    0.0000" << std::setw(3) << atom_symbol
           << "  0" << std::setw(3) << chg_ele << "  0  0  0  0  0  0  0  0  0  0"
           << "\n";
      }

      for (size_t i = 0; i < _g.numAtom; ++i) {
        for (size_t j = i + 1; j < _g.numAtom; ++j) {
          if (_g.beta[i][j] > 0
            && i != e_id  // ignore e*
            && j != e_id  // ignore e*
          ) {
            size_t i_real = i;
            size_t j_real = j;

            if (i > e_id) --i_real;
            if (j > e_id) --j_real;

            output << std::setw(3) << i_real + 1 << std::setw(3) << j_real + 1 << std::setw(3)
                   << _g.beta[i][j] << "  0  0  0  0"
                   << "\n";
          }
        }
      }

      size_t i_real = _g.adj[e_id][0];
      size_t j_real = _g.adj[e_id][1];
      if (i_real > e_id) --i_real;
      if (j_real > e_id) --j_real;

      output << std::setw(3) << i_real + 1 << std::setw(3) << j_real + 1 << std::setw(3)
                   << 1 << "  0  0  0  0"
                   << "\n";

      size_t n_chg = 0;
      vector <int> vector_chg;
      vector_chg.clear();

      for (size_t i = 0; i < _g.numAtom; ++i){

          if (i != e_id && _g.charge[i] != 0){
            ++n_chg;
            if (i < e_id){
              vector_chg.push_back(i + 1);
              vector_chg.push_back(_g.charge[i]);
            } else {
              vector_chg.push_back(i);
              vector_chg.push_back(_g.charge[i]);
            }
            
          }
        
      }

      if (n_chg == 0){
        output << "M  END"
               << "\n";
      } else {
        output << "M  CHG" << std::setw(3) << n_chg;
        for (auto& tmp : vector_chg){
          output << std::setw(4) << tmp;
        }
        output << "\n";
        output << "M  END"
               << "\n";
      }
      
      output << "$$$$"
             << "\n";
      output.close();
      /////////////////////////////////////////////////////////////////////
    } 
    // else if (e_id_vector.size() == 2){

    //   int v_id_1 = _g.adj[e_id_vector[0]][0];
    //   int v_id_2 = _g.adj[e_id_vector[1]][0];
    //   int e_id_1 = e_id_vector[0];
    //   int e_id_2 = e_id_vector[1];

    //   /////////////////////////////////////////////////////////////////////
    //   // output SDF file
    //   ofstream output(output_file, ios::app);

    //   size_t real_n = _g.numAtom - 2;
    //   size_t real_m = 0;

    //   for (Vertex u = 0; u < _g.numAtom; ++u){
    //     for (Vertex v = u + 1; v < _g.numAtom; ++v){
    //       if (_g.beta[u][v] > 0){
    //         ++real_m;
    //       }
    //     }
    //   }
    //   real_m -= 2; // erase the edge incident to e*

    //   output << _g.CID << "\n";
    //   output << ""
    //          << "\n";
    //   output << ""
    //          << "\n";
    //   output << std::setw(3) << real_n << std::setw(3) << real_m
    //          << "  0  0  0  0  0  0  0  0999 V2000 "
    //          << "\n";

    //   map <Vertex, int> v_map;
    //   v_map.clear();

    //   for (size_t i = 0; i < _g.numAtom; ++i) {
    //     string atom_symbol = _g.alpha[i];
    //     if (atom_symbol == "e*") continue; // ignore e*
    //     output << "    0.0000    0.0000    0.0000" << std::setw(3) << atom_symbol
    //            << "  0  0  0  0  0  0  0  0  0  0  0  0"
    //            << "\n";

    //     int v_map_size = v_map.size() + 1;
    //     v_map.emplace(i, v_map_size);
    //   }

    //   for (size_t i = 0; i < _g.numAtom; ++i) {
    //     for (size_t j = i + 1; j < _g.numAtom; ++j) {
    //       if (_g.beta[i][j] > 0 
    //         && find(e_id_vector.begin(), e_id_vector.end(), i) == e_id_vector.end()
    //         && find(e_id_vector.begin(), e_id_vector.end(), j) == e_id_vector.end()
    //       ) {
    //         size_t i_real = v_map.at(i);
    //         size_t j_real = v_map.at(j);

    //         output << std::setw(3) << i_real << std::setw(3) << j_real << std::setw(3)
    //                << _g.beta[i][j] << "  0  0  0  0"
    //                << "\n";
    //       }
    //     }
    //   }

    //   size_t i_real = v_map.at(v_id_1);
    //   size_t j_real = v_map.at(v_id_2);


    //   output << std::setw(3) << i_real << std::setw(3) << j_real << std::setw(3)
    //                << 1 << "  0  0  0  0"
    //                << "\n";

    //   size_t n_chg = 0;
    //   vector <int> vector_chg;
    //   vector_chg.clear();

    //   for (size_t i = 0; i < _g.numAtom; ++i){

    //       if (find(e_id_vector.begin(), e_id_vector.end(), i) == e_id_vector.end()
    //            && _g.charge[i] != 0){
    //         ++n_chg;
    //         vector_chg.push_back(v_map.at(i));
    //         vector_chg.push_back(_g.charge[i]);  
    //       }
        
    //   }

    //   if (n_chg == 0){
    //     output << "M  END"
    //            << "\n";
    //   } else {
    //     output << "M  CHG" << std::setw(3) << n_chg;
    //     for (auto& tmp : vector_chg){
    //       output << std::setw(4) << tmp;
    //     }
    //     output << "\n";
    //     output << "M  END"
    //            << "\n";
    //   }
      
    //   output << "$$$$"
    //          << "\n";
    //   output.close();
    //   /////////////////////////////////////////////////////////////////////

    // } 
    else {
      cout << "CID : " << _g.CID << " has " << e_id_vector.size() << " e*!!!! ignored." << endl;
    }

    

  }

  return 0;
}
