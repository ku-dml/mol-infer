# -*- coding: utf-8 -*-

"""
The program predicts the target values of compounds in argv[3],
using the ANN specified by argv[1] and argv[2]. 

Author:    Discrete Mathematics Lab,
           Department of Applied Mathematics and Physics,
           Graduate School of Informatics,
           Kyoto University
"""

import pandas as pd
import numpy as np
import sys

class ANN:
    """
    An encapsulating class to conveniently store
    information about a trained neural network
    """

    def __init__(self, weight_tensor=list(), bias_matrix=list()):
        # prepare a list to store all the nodes
        self.nodes = list()
        # prepare a list to store sets of nodes by layer
        self.layers = list()
        # a dictionary to bind pairs of vertices to weights
        self.weights = dict()
        # a dictionary to bind each vertex to a bias
        self.biases = dict()
        # a dictionary to store a list of the in neighbors
        # for each vertex not in the input layer
        self.in_neighbors = dict()

        for l, bias_layer in enumerate(bias_matrix):
            layer_nodes = list()
            for k, bias in enumerate(bias_layer):
                node = (l+1, k+1)
                self.nodes.append(node)
                self.biases[node] = bias
                layer_nodes.append(node)
            self.layers.append(layer_nodes)

        """
        Store the nodes in lists by layer
        """
        for l, layer in enumerate(self.layers):
            try:
                for node in self.layers[l+1]:
                    self.in_neighbors[node] = layer
            except IndexError:
                # We are at the last layer that
                # does not have a succeeding layer
                pass
            except ValueError:
                # This should not happen
                print("""
                Received Value Error
                in distributing ANN nodes by layers
                """)
                sys.exit()

        """
        Store the weights as a dictionary
        of ordered node pairs
        """
        for l, weight_matrix in enumerate(weight_tensor):
            for i, weight_row in enumerate(weight_matrix):
                for j, weight in enumerate(weight_row):
                    # source node in current layer
                    u = (l+1, i+1) # mind the +1 offset in naming
                    # target node in next layer
                    v = (l+2, j+1)
                    self.weights[(u, v)] = weight

    def activation(self, x, func=None):
        """
        Define an activation function for the ANN
        If no particular function is supplied, then
        ReLU is used by default
        """
        if func:
            # This means that some activation function was supplied
            return func(x)
        else:
            # Otherwise, we use ReLU
            return max(0, x)

    @property
    def input_layer(self):
        """
        Simply return the first layer
        """
        return self.layers[0]

    @property
    def output_layer(self):
        """
        Return the last layer
        """
        return self.layers[-1]

    @property
    def hidden_layers(self):
        """
        Return an iterator of the internal layers,
        i.e. Layers other than index 0 and -1
        """
        for l in self.layers[1:-1]:
            yield l

    @property
    def non_input_nodes(self):
        """
        Return a generator of all the nodes in all
        layers except the first one
        """
        for layer in self.layers[1:]:
            for node in layer:
                yield node

    def propagate(self, fv):
        """
        Predict using trained model

        Arguments:
            fv {array_like} -- feature vector (input layer)
        Returns:
            pred_val {scalar object} -- predicted value (output layer)

        NOTE
        input layer:
            y[0] = fv

        hidden layer and output layer:
            x[i+1] = dot(w[i+1], y[i]) + b[i+1]
            y[i+1] = activation(x[i+1])
                i = 1, 2, ..., num_layer - 1
        """

        # list of nodes in current layer
        cur_layer = self.input_layer

        # vector in current layer {numpy.ndarray}
        vec_y = np.array(fv)

        for next_layer in self.layers[1:]:
            # prepare weight matrix {numpy.ndarray}
            weight_matrix = np.array([
                                [self.weights[(u, v)] for u in cur_layer]
                                for v in next_layer
                            ])

            # prepare bias vector {numpy.ndarray}
            bias_vector = np.array([
                            self.biases[v]
                            for v in next_layer
                        ])
            
            vec_x = np.dot(weight_matrix, vec_y) + bias_vector
            vec_y = np.array([
                        self.activation(x)
                        for x in vec_x
                    ])
            # update current layer
            cur_layer = next_layer

        return list(vec_y)


def line2num(line, numtype=float, sep=" "):
    line = line.strip()
    return [numtype(num) for num in line.split(sep)]


def read_trained_ANN(weight_filename, bias_filename):
    """
    Given filenames to files that
    contain the values for the weights and biases
    of a trained ANN, read the values form the files,
    contruct an ANN and return it
    """
    try:
        with \
        open(weight_filename, "r") as wf, \
        open(bias_filename, "r") as bf:
            layer_sizes = line2num(wf.readline(), int)
            weight_tensor = list()
            for ell in layer_sizes[:-1]:
                weight_matrix = list()
                for _ in range(ell):
                    weight_line = wf.readline()
                    weight_matrix.append(line2num(weight_line))
                weight_tensor.append(weight_matrix)
            bias_matrix = list()
            """
            The input layer has no bias, and
            we define it to be 0
            """
            bias_matrix.append([0]*layer_sizes[0])
            for ell in layer_sizes[1:]:
                bias_row = list()
                for _ in range(ell):
                    bias_line = bf.readline().strip()
                    bias_row.append(float(bias_line))
                bias_matrix.append(bias_row)
    except:
        print(
        """
        An error occured when trying to read ANN data from files
        {}
        {}
        Now will terminate
        """.format(weight_filename, bias_filename))
        sys.exit()

    return weight_tensor, bias_matrix


def main(argv):
    """
    assumptions:
        - argv[0] = this script's name
        - argv[1] = filename for weights of NN
        - argv[2] = filename for biases of NN
        - argv[3] = filename for fv format file (csv)
        - argv[4] = filename for writing the result
    
    Make a prediction for each fv in the fv file
    """
    try:
        w, b = read_trained_ANN(argv[1], argv[2])
        ann = ANN(w, b)
    except:
        sys.stderr.write('error: failed to read ANN data.\n')
        exit(1)

    try:
        df = pd.read_csv(argv[3], sep = ",")
        CIDs = df.values[:, :1]
        fvs = df.values[:, 1:]
        ans = list()
    except:
        sys.stderr.write('error: failed to read csv file.\n')
        exit(1)

    try:
        for fv in fvs:
            ans.append(ann.propagate(fv))
    except:
        sys.stderr.write('error: failed to make prediction.\n')
        exit(1)
        
    try:
        f = open(argv[4],'w')
        f.write('CID,a\n')
        for (cid,a) in zip(CIDs,ans):
            f.write('{},{}\n'.format(int(cid[0]), a[0]))
        f.close()
    except:
        sys.stderr.write('error: failed to write the prediction result.\n')
        exit(1)

    N = len(CIDs)
    print('The prediction results for {} compounds in {} are written to {}.\n'.format(N,argv[3],argv[4]))
    

if len(sys.argv) != 5:
    sys.stderr.write('usage {}: (ANN_weight)(ANN_bias)(INPUT.csv)(OUTPUT.csv)\n'.format(sys.argv[0]))
    exit(1)
    
main(sys.argv)
