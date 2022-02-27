# -*- coding: utf-8 -*-

"""
This file implements
 - a class to store the architecture of a multilayer perceptron (ANN)
 - functions to initialize vatiables and build a pulp MILP model for
   inverting an ANN of a given architecture

Author:    Discrete Mathematics Lab,
           Department of Applied Mathematics and Physics,
           Graduate School of Informatics,
           Kyoto University
"""

# standard imports
import sys
import pandas as pd
import numpy as np

# import pulp for MILP modeling
import pulp


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
            
            if next_layer == self.output_layer:
                vec_y = np.array([x for x in vec_x])
            else:
                vec_y = np.array([
                            self.activation(x)
                            for x in vec_x
                        ])
                
            # update current layer
            cur_layer = next_layer

        return list(vec_y)


def initialize_constants(ann, training_data):
    """
    Given a trained ANN and a matrix of training data
    where each row is a single sample,
    calculate the ranges of values that can be passed
    to each node in the ANN
    """

    # cast the training data as a numpy array
    td = np.array(training_data)
    # Check if the training_data's feature vectors have exactly the same
    # number of elements as nodes in the input layer
    if td.shape[1] != len(ann.input_layer):
        print("###############")
        print("""
        Vector length mismatch between training data
        and ANN input layer!
        """)
        sys.exit()
    # prepare dictionaries to store values
    a_low = dict()
    a_high = dict()
    a = dict()
    b_low = dict()
    b_high = dict()
    b_hat = dict()
    b = dict()
    c = dict()
    z13 = dict()
    # first, calculate b_low, b_high values for the
    # nodes of the input layer
    for node, column in zip(ann.input_layer, td.T):
        b_low[node] = min(column)
        b_high[node] = 1000
        b[node] = (b_low[node], 0, b_high[node])

    # calculate the ranges for the remaining layers
    for v in ann.non_input_nodes:
        a_low[v], a_high[v] = 0, 0
        for u in ann.in_neighbors[v]:
            w_uv = ann.weights[(u, v)]
            if w_uv > 0:
                a_low[v] += w_uv*b_low[u]
                a_high[v] += w_uv*b_high[u]
            else:
                a_high[v] += w_uv*b_low[u]
                a_low[v] += w_uv*b_high[u]
        # if ann.biases[v] >= 0:
        #     a_high[v] += ann.biases[v]
        # else:
        #     a_low[v] += ann.biases[v]
        a_high[v] += ann.biases[v]
        a_low[v] += ann.biases[v]

        # zip the a_low, a_high variables into a single triplet
        a[v] = (a_low[v], 0, a_high[v])

        b_low[v] = 0

        ##################
        if v in ann.output_layer:
            b_low[v] = -1000
        ##################

        if a_high[v] > 0:
            b_high[v] = a_high[v]
        else:
            b_high[v] = 0
        # zip the b_low, b_high variables into a single triplet
        b[v] = (b_low[v], 0, b_high[v])
        b_hat[v] = 2*a_high[v] - a_low[v]
        c[v] = (0, 1)
        z13[v] = (1, 0)

        # print(v, " a_high[v] = ", a_high[v], " a_low[v] = ", a_low[v])

    return a, b, b_hat, c, z13


def initialize_lp_variables(ann, a, b, forbidden_node, property_name="def"):
    """
    Given a trained ANN,
    initialize variables for each node.
    Dictionaries a and b give the ranges for the variables
    """
    low = 0
    high = 2

    """
    A dictionary to store x_v variable for each node v not in the input layer.
    x_v is real-valued and a[v][low] <= x_v <= a[v][high]
    """
    x = {node: pulp.LpVariable("ann_x{}_{}".format(node, property_name),
                               a[node][low],
                               a[node][high])
         for node in ann.non_input_nodes}

    """
    A dictionary to store y_v variable for each node v in ann.
    y_v is real-valued and 0 <= y_v <= b[v][high]
    """
    y = {node: pulp.LpVariable("ann_y{}_{}".format(node, property_name),
                               b[node][low],
                               b[node][high])
         for node in ann.nodes if node not in forbidden_node}
    """
    A dictionary to store binary z_v variable
    for each node v not in the input layer.
    """
    z = {node: pulp.LpVariable("ann_z{}_{}".format(node, property_name), 0, 1, pulp.LpBinary)
         for node in ann.non_input_nodes}

    return x, y, z


def build_MILP_ReLU(
    model    : pulp.LpProblem,
    ann      : ANN, 
    variables: tuple, 
    constants: tuple, 
    target   : tuple,
    eps      : float,
    forbidden_node,
    property_name = "def" 
    )-> pulp.LpProblem:
    """
    Given matrices of variables and constants,
    construct an MILP model in PuLP according to the
    MILP formulation for ReLU activation functions,
    Akutsu and Nagamochi, 2019
    """

    if type(target) == int or type(target) == float:
        target = [target, ]
    """
    First, check if the last (output) layer of ann
    has the same size as the target data
    """
    if len(ann.output_layer) != len(target):
        print("""
        Error: The size of the output layer and the
        target data do not match!
        The program will now exit.
        """)
        sys.exit()

    x, y, z2 = variables  # unzip the variables
    a, b, b_hat, c, z13 = constants  # unzip the constants
    """
    For convenience, zip the constants from z13 and variables from z2
    into a single structure
    """
    z = dict()
    for node in ann.non_input_nodes:
        z[node] = (z13[node][0], z2[node], z13[node][1])

    # Constraint on the activation of ann
    for v in ann.non_input_nodes:
        in_u = ann.in_neighbors[v]
        w_v = [ann.weights[(u, v)] for u in in_u if u not in forbidden_node]
        in_y = [y[u] for u in in_u if u not in forbidden_node]

        model += \
            x[v] == pulp.lpDot(in_y, w_v) + ann.biases[v], \
            "Output_node_{}_{}".format(v, property_name)

        if v in ann.output_layer:
            model += y[v] == x[v], "ReLU_y{}_{}".format(v, property_name)
        else:
        ###################
            model += \
                x[v] - a[v][1] <= (a[v][2] - a[v][0]) * z[v][1], \
                "ReLU_x{}_ub_{}".format(v, property_name)

            if a[v][2] > 0:
                model += \
                    x[v] - a[v][1] >= (a[v][0] - a[v][2]) * (1 - z[v][1]), \
                    "ReLU_x{}_lb_{}".format(v, property_name)

                for i in [0, 1]:
                    model += \
                        y[v] <= c[v][i] * (x[v] - a[v][i]) + b[v][i] \
                        + b_hat[v] * (1 + z[v][i + 1] - z[v][i]), \
                        "ReLU_y{}_{}_ub_{}".format(v, i, property_name)

                    model += \
                        y[v] >= c[v][i] * (x[v] - a[v][i]) + b[v][i] \
                        - b_hat[v] * (1 + z[v][i + 1] - z[v][i]), \
                        "ReLU_y{}_{}_lb_{}".format(v, i, property_name)
            else:
                model += z[v][1] == 0, "ReLU_x{}_lb_{}".format(v, property_name)
                model += y[v] == 0, "ReLU_y{}_{}".format(v, property_name)

        #######################


    for out_node, tv in zip(ann.output_layer, target):
        #######################
        if tv >= 0:
            model += \
                y[out_node] >= tv * (1 - eps), \
                "lower_bound_target_{}_{}".format(out_node, property_name)

            model += \
                y[out_node] <= tv * (1 + eps), \
                "upper_bound_target_{}_{}".format(out_node, property_name)
        else:
            model += \
                y[out_node] >= tv * (1 + eps), \
                "lower_bound_target_{}_{}".format(out_node, property_name)

            model += \
                y[out_node] <= tv * (1 - eps), \
                "upper_bound_target_{}_{}".format(out_node, property_name)
        #######################

    # Finally, return the built model
    return model


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


def read_fv_descriptors(training_data_filename):
    """
    Read the textual names for the descriptors
    in the feature vector, as the first line from
    the supplied training data file.
    Return a list of these names as strings
    """
    try:
        data_frame = pd.read_csv(training_data_filename, sep=",")
    except:
        print("""
        Error reading the file {}
        with pandas.
        """.format(training_data_filename))
        sys.exit()
    return list(data_frame.columns)


def read_training_data(training_data_filename):
    """
    Given a set of feature vectors as training data,
    return a matrix that contains one feature vector
    per row
    """
    try:
        data_frame = pd.read_csv(training_data_filename, sep=",")
    except:
        print("""
        Error reading the file {}
        with pandas.
        """.format(training_data_filename))
        sys.exit()
    # print(data_frame.values) # testing

    try:
        table = data_frame.values  # numpy array
        table = table[:, 1:]
    except:
        print("""
        Exception in converting the dataframe
        to a numpy array, file
        {}
        """.format(training_data_filename))
        sys.exit()
    # Success
    return table        
 
def get_input_layer_variables(ann, variables, descriptors, forbidden_node):
    """
    Given a list of descriptor names and a tuple of
    variables of the MILP_Relu model,
    create and return a dictionary that for each descriptor
    gives the corresponding variable
    """
    # unpack the variables
    x, y, z = variables
    # Initialize an empty dictionary
    sol = dict()
    # Iterate over the input layer and the list of descriptors
    # at the same time
    for v, name in zip(ann.input_layer, descriptors):
        if v not in forbidden_node:
            sol[name] = y[v]
    return sol

# Testing
def test():
    w, b = read_trained_ANN(
        "../../test_files/fv4_Kow_weight_5.txt", 
        "../../test_files/fv4_Kow_bias_5.txt" )
    
    ann = ANN(w, b)
    #     for key, val in ann.in_neighbors.items():
    #         print("{} : {}".format(key, val))
    
    training_data = read_training_data("../../test_files/fv4_Kow.txt")
    descriptors = read_fv_descriptors("../../test_files/fv4_Kow.txt")
    # print(descriptors)
    
    consts = initialize_constants(ann, training_data)
    a, b, b_hat, c, z  = consts
    #     print("a:\n", a)
    #     print("b:\n", b)
    #     print("b_hat:\n", b_hat)
    #     print("c:\n", c)
    #     print("z:\n", z)
    #     znp = np.array(z)
    #     print(znp)

    LpVars = initialize_lp_variables(ann, a, b)
    #     print(LpVars)
    model = pulp.LpProblem("Test_A", pulp.LpMinimize)
    LpModel = build_MILP_ReLU(model, ann, LpVars, consts, 5.0, 0.02)
    # print(LpModel)
    # LpModel.writeLP("test_LP_file")

    print("Start solving")
    LpModel.solve()
    print("Solved")
    x, y, z = LpVars
    sol_vars = get_input_layer_variables(ann, LpVars, descriptors)
        
    sol_y = [var.value() for var in sol_vars.values()] 
    print(sol_y)
    res = ann.propagate(sol_y)
    print(res)
    
# test()
