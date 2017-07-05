# -*- coding: utf-8 -*-
"""
Created on Sun Jul  2 19:55:19 2017

@author: Turtleflying
"""
import numpy as np


def Measurement(mu):
    if 3*np.sqrt(np.pi)/2 > mu > np.sqrt(np.pi)/2:
        mu = mu - np.sqrt(np.pi)
        error = -1
    elif -3*np.sqrt(np.pi)/2 < mu < -np.sqrt(np.pi)/2:
        mu = mu + np.sqrt(np.pi)
        error = -1
    else:
        mu = mu
        error = +1
    return mu, error


import networkx as nx
import time
from scipy.stats import norm
import matplotlib.pyplot as plt
for iter1 in range(5):  # This loop is over the dimension of the code

    rate_list = []      # a list used to store the error rate
    sigma_list = []     # a list used to stor the standard deviation corresponding to error rate
    for iter2 in range(
            5):  # This loop is over the standard deviation sigma of the shift errors
        start_time = time.time()
        logical_error1 = 0    # This variable is used to count the number of wrong error correction
        logical_error2 = 0
        sum = 100
        dim = 8 + iter1*4           # The dimension of the code
        sigma = 0.55 + iter2*0.005  # The standard deviation of the shift errors
        #sigma = 0.4
        print('dimension', dim)
        print('sigma', sigma)
        print('sum =', sum)
        print('-------------')

        for iter in range(sum):

            #---------------------------------------------------------------------------#
            #---------------------------------------------------------------------------#
            #---------------------------------------------------------------------------#

         # Initialize the code with no weight
            G1 = nx.grid_2d_graph(dim, dim, True)
            G2 = nx.grid_2d_graph(dim, dim, True)

            for (u, v) in G1.edges():
                shift_error = np.random.normal(0, sigma)
                measure, error = Measurement(shift_error)

                rate = norm.pdf(measure,
                                0,
                                sigma) / (norm.pdf(measure,
                                                   0,
                                                   sigma) + norm.pdf(measure + np.sqrt(np.pi),
                                                                     0,
                                                                     sigma) + norm.pdf(measure - np.sqrt(np.pi),
                                                                                       0,
                                                                                       sigma))

                if rate < 0.5:
                    G1.edge[u][v]['weight'] = 0
                else:
                    G1.edge[u][v]['weight'] = np.log2(rate) - np.log2(1 - rate)

                G2.edge[u][v]['weight'] = error
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#

# If the product of four neighbouring edges' weights  is -1,
# the vertice will be regarded as a defect

            syndrome_list = []
            for u in G1.nodes():
                neighbors = nx.all_neighbors(G1, u)

                syndrome = 1
                for v in neighbors:
                    syndrome *= G2.edge[u][v]['weight']
                if syndrome < 0:
                    syndrome_list.append(u)
 #---------------------------------------------------------------------------#
 #---------------------------------------------------------------------------#
 #---------------------------------------------------------------------------#
 #---------------------------------------------------------------------------#

#  Initialize another graph to match all those defects with minimum weight
            G3 = nx.Graph()
            G3.add_nodes_from(syndrome_list)
            for u in syndrome_list:
                for v in syndrome_list:
                    if u == v:
                        pass
                    else:
                        w = nx.algorithms.shortest_paths.weighted.dijkstra_path_length(
                            G1, u, v, weight='weight')

                        # This '50-w' will transform the max-weight matching
                        # into minimum weight matching
                        G3.add_edge(u, v, weight=50 - w)

            matching_list = nx.max_weight_matching(G3)

#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#

# After minimum weight matching, we correct each qubit on the path which
# connects each pair of defect

            for u in syndrome_list:

                start = u
                end = matching_list[start]
                syndrome_list.remove(end)
                path = nx.dijkstra_path(G1, start, end, weight='weight')

                length = len(path)
                for i in range(length-1):
                    u = path[i]
                    v = path[i+1]
                    G2.edge[u][v]['weight'] *= (-1)

#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#

# To check whether we do the right correction, which means no logical error
            check = 1
            for i in range(dim):
                for j in range(dim):
                    i2 = (i + 1) % dim
                    check *= G2[(i, j)][(i2, j)]['weight']
                if check < 0:
                    logical_error1 += 1
                    break
                check = 1

            check = 1
            for j in range(dim):
                for i in range(dim):
                    j2 = (j + 1) % dim
                    check *= G2[(i, j)][(i, j2)]['weight']
                if check < 0:
                    logical_error2 += 1
                    break
                check = 1

        error_rate = (logical_error1 + logical_error2) / sum
        rate_list.append(error_rate)
        sigma_list.append(sigma)
        print('error rate',  logical_error1 / sum)
        print('error rate2', logical_error2 / sum)
        end_time = time.time()
        print(end_time - start_time)
        print('-----------------------------------')
    x = ['g--', 'b--', 'r--', 'y--', 'w--', 'k--']
    plt.plot(sigma_list, rate_list, x[iter1], label='original')

plt.show()
