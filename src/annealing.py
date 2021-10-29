'''
Annealing class for size-bases constrained simulated annealing of community detection

@author Daniel Kaiser. Created on 2021-03-16 17:16 EST.
'''
# =========================== SCRIPT INFO =============================
__author__ = "Daniel Kaiser"
__credits__ = ["Daniel Kaiser", "Sangpil Youm"]
__maintainer__ = "Daniel Kaiser"
__email__ = "kaiserd@iu.edu"
__version__ = "0.2.3"


# ========================= IMPORTS & GLOBAL ==========================
# ------- System Imports ---------
import sys
import os

# ------- General Imports ---------
import numpy as np
import pandas as pd
import random
from time import perf_counter

# ------- Network Imports ----------
import igraph as ig

# ------- Plotting Imports ---------
import matplotlib.pyplot as plt
import pandas as pd


# -------- Module Imports ---------


# ~~~~~ Debug ~~~~~~~


# ====================== FUNCTIONS & CLASSES ==========================
class Annealing:
    '''
    Class to do size-based constrained simulated annealing for network community detection.\n

    Only valid on single-layer simple undirected networks. \n

    Currently valid for .gml, .edgelist, XXX files.
    '''

    def __init__(self, graph, prob_ms, m_size):
        '''
        Class constructor.

        Inputs :
            graph : str or igraph GraphBase object
                File location or igraph GraphBase object to be annealed.

        Returns :
            None
        '''

        if type(graph) == str:
            self.graph = ig.Graph.Read(graph) # Graph object
        else:
            self.graph = graph # Graph object

        self.n = self.graph.vcount() # Number of nodes

        # ! DEBUG
        self.mod_list = [-1, ]
        self.mod_diff = []

        self.prob_ms = prob_ms
        self.current = 0
        self.comm_formation = None
        self.svg_group = [None for i in range(m_size)]

    def helper_init_partition(self):
        '''
        Documentation here, I guess
        '''
        # Declare initial partition
        init_num_groups = np.random.choice(range(1,self.max_size))
        #change here
        #init_num_groups = int(np.ceil(self.n / self.max_size)) # To ensure within max size
        perm = np.random.permutation(self.n)
        init_partition = [perm[x] % init_num_groups for x in range(self.n)] # initially a bunch of max size
        init_groups = [[node for node, idx in enumerate(init_partition) if idx == k] for k in range(init_num_groups)]
        new_groups = []
        for c, group in enumerate(init_groups):
            #change here
            #num_move = np.random.randint(0,len(group)-self.min_size)
            num_move = np.random.randint(0,len(group)-self.min_size+1)
            new_groups.append(group[:num_move])
            init_groups[c] = group[num_move:]

        groups = init_groups + new_groups
        self.comm = {}
        for group_idx, group in enumerate(groups):
            for node in group:
                self.comm[node] = group_idx

        self.comm = {node: self.comm[node] for node in range(self.n)}
        self.modularity = ig.VertexClustering(self.graph, self.comm.values()).q # For movement calculations
        self.mod_list.append(self.modularity)

    def helper_global_merge(self, temp, selected_comm_idx, selected_comm):
        comm = self.comm.copy()

        # Select destination community
        new_selected_comm_idx = np.random.randint(low=0, high=max(self.comm.values())+1)
        new_selected_comm = [i for i, j in self.comm.items() if j == new_selected_comm_idx]

        # Check the merge wouldn't be too large
        if (len(new_selected_comm) + len(selected_comm) <= self.max_size) \
            & (selected_comm_idx != new_selected_comm_idx):

            # Try moving nodes
            for node in selected_comm:
                comm[node] = new_selected_comm_idx

            # Check modularity
            new_modularity = ig.VertexClustering(self.graph, comm.values()).q

            # If the move is better, assign move
            if new_modularity > self.modularity:
                self.modularity = new_modularity

            # If not better, check temperature and perhaps assign move anyways
            elif np.random.rand() < np.exp(-abs(self.modularity - new_modularity)/temp):
                self.modularity = new_modularity

            # If move not accepted, undo it
            else:
                comm = self.comm

        self.comm = comm

    def helper_global_split(self, temp, selected_comm_idx, selected_comm):
        comm = self.comm.copy()

        # Find cut and check it is nontrivial
        subgraph = ig.VertexClustering(graph=self.graph, membership=comm).subgraph(selected_comm_idx)
        if len(subgraph.es) > 0:
            mc = subgraph.mincut()
            if len(mc.partition) == 1:
                return

            # Do partition adjustment in split
            new_comm = mc.partition[1]
            new_comm_idx = max(comm.values()) + 1
            for node in new_comm:
                comm[node] = new_comm_idx

            # Retrieve modularity change and check state
            new_modularity = ig.VertexClustering(self.graph, comm.values()).q

            # If the move is better, keep it
            if new_modularity > self.modularity:
                self.modularity = new_modularity

            # If not better, check temperature and perhaps assign move anyways
            elif np.random.rand() < np.exp(-abs(self.modularity - new_modularity)/temp):
                self.modularity = new_modularity

            # If move not accepted, undo it
            else:
                comm = self.comm

            # Track (new) communities
            self.comm = comm
            return
        else:
            return

    def local_move(self, temp):
        comm = self.comm.copy()

        ######  TODO LOOK FOR QUICKER WAYS TO DO THIS ######
        selected_comm_idx = np.random.randint(low=0, high=max(self.comm.values())+1)
        selected_comm = [node for node, node_comm_idx in self.comm.items()
                            if node_comm_idx == selected_comm_idx]

        # If there is room in constraint to do the local move
        if self.min_size < len(selected_comm) <= self.max_size:
            # Select destination community
            new_selected_comm_idx = np.random.randint(low=0, high=max(self.comm.values())+1)
            new_selected_comm = [i for i, j in self.comm.items() if j == new_selected_comm_idx]

            # Check destination has room
            if len(new_selected_comm) < self.max_size:
                # Get random node index to move
                node = np.random.choice(selected_comm)

                # Check move modularity
                comm[node] = new_selected_comm_idx
                new_modularity = ig.VertexClustering(self.graph, comm.values()).q
                thresh = np.exp(-abs(self.modularity - new_modularity)/temp)

                # ! Debug
                self.mod_diff.append(new_modularity - self.modularity)

                # If the move is better, assign move
                if new_modularity >= self.modularity:
                    self.modularity = new_modularity

                # If not better, check temperature and perhaps assign move anyways

                elif np.random.rand() < thresh:
                    self.modularity = new_modularity

                # If move not accepted, undo it
                else:
                    comm = self.comm

        # Note comm is only different if a move was done
        self.comm = comm
        return

    def global_move(self, temp, p_m = 0.495, p_s = 0.495):
        p = np.random.rand()
        selected_comm_idx = np.random.randint(low=0, high=max(self.comm.values())+1)
        selected_comm = [i for i, j in self.comm.items() if j == selected_comm_idx]

        # Do merge move
        if p < p_m:
            self.helper_global_merge(temp=temp, selected_comm_idx=selected_comm_idx, selected_comm=selected_comm)

        # Do split move
        elif p_m <= p < p_m + p_s:
            self.helper_global_split(temp=temp, selected_comm_idx=selected_comm_idx, selected_comm=selected_comm)

        # Do local move
        elif p_s <= p:
            self.local_move(temp)

    def anneal(self, min_size, max_size,
        temp = 100.0, cooling_factor = 0.995, stop_epsilon = 1e-8, stop_condition=5,
        p_m = 0.40, p_s = 0.60):
        # Book-keeping
        self.min_size = int(min_size)
        self.max_size = int(max_size)
        max_temp = temp

        # Get random initial partition
        self.helper_init_partition()

        # Annealing process
        # len(set(self.mod_list[-stop_condition:])) > 1
        valid = True

        while temp > stop_epsilon and valid == True  :

            prob = np.random.rand()
            # Do n global moves
            if prob > self.prob_ms :
                for _ in range(self.n):
                    self.global_move(temp=temp, p_m=p_m, p_s=p_s)

            # Do n^2 local moves
            else :
                for _ in range(self.n*2):
                    self.local_move(temp=temp)

            self.mod_list.append(self.modularity)

            # get new mapping group
            set_group = self.comm_config(self.comm)
            #chceck uniqueness
            valid = self.unique(set_group)

            if valid :
                self.add_comm()


            # Cool system
            temp *= cooling_factor

            # ! Debug
            pct_cooled = 100*(max_temp - temp)/(max_temp-stop_epsilon)
            if np.random.rand() <= 0.01:
                print("System temperature {:.2f} percent cooled".format(pct_cooled))
                print("System has {:.4f} modularity".format(self.modularity))

        return self.comm, self.modularity

    def add_comm(self):

        # max number of hash

        new_nodes = self.comm_config(self.comm)
        self.svg_group[self.current] = new_nodes
        #print(self.svg_group)

        self.comm_formation = None
        self.current +=1


    def comm_config(self, comm):

        comm_list = list(comm.values())
        comm_arr = np.array(comm_list)
        len_comm = len(comm_list)
        mapping = -1 * np.ones(len_comm)

        idx = 0
        for i in range(len_comm) :
            if mapping[comm_arr[i]] <0 :
                mapping[comm_arr[i]] = idx
                idx +=1

        comm_set = mapping[comm_arr].tolist()
        self.comm_formation = str(comm_set)

        return self.comm_formation

    def unique(self, new_comm):

        for groups in self.svg_group[-1000:] :

            if groups == None :
                break

            if groups == new_comm :
                return False

        return True



if __name__ == '__main__':
    # Karate
    '''
    G = ig.Graph.Read('/Users/macbookpro/Project/size-csa/data/real_networks/karate.gml', format='gml')
    ann = Annealing(G,0.10)
    p, q = ann.anneal(min_size=1, max_size=18, stop_condition=5, stop_epsilon=1e-4, cooling_factor=0.90)
    print(p.values())
    print(max(p.values())+1, q)
    '''

    # Dolphins

    '''
    G = ig.Graph.Read('/Users/macbookpro/Project/size-csa/data/real_networks/dolphins.gml', format='gml')
    ann = Annealing(G,0.90,10000)
    p, q = ann.anneal(min_size=10, max_size=20, stop_condition=5, stop_epsilon=1e-4, cooling_factor=0.9)
    print(max(p.values())+1, q)
    '''

    # Les Mis
    '''
    G = ig.Graph.Read('/Users/macbookpro/Project/size-csa/data/real_networks/lesmis.gml', format='gml')
    ann = Annealing(G,0.10, 10000)
    p, q = ann.anneal(min_size=1, max_size=20, stop_condition=5, stop_epsilon=1e-4, cooling_factor=0.995)
    print(max(p.values())+1, q)
    '''
    num_group= []
    modul = []
    result_set = {}
    size_l = []
    mu_l = []

    for size in [500] :
        for mu in [0.1]:
            for i in range(5) :
                print("try {}...".format(i))
                G = ig.Graph.Read('/Users/macbookpro/Project/size-csa/data/benchmarks/ntwk'+str(size)+'_'+str(mu)+'.gml', format='gml')
                ann = Annealing(G, 0.10,1000000)
                time1 = perf_counter()
                p, q = ann.anneal(min_size=1, max_size=500, stop_condition=10, stop_epsilon=1e-4, cooling_factor=0.99)
                time2 = perf_counter()
                print(max(p.values()) + 1, q)

                final_set = {}
                f_nodes = list(p.values())
                for node, group in enumerate(f_nodes) :
                    if group not in final_set.keys() :
                        final_set[group] = [node]
                    else :
                        final_set[group].append(node)

                print(final_set)
                print("Running time {}".format(time2-time1))

                num_group.append(max(p.values()) + 1)
                modul.append(q)
                size_l.append(size)
                mu_l.append(mu)

    result_set['num_group'] = num_group
    result_set['modularity'] = modul
    result_set['size'] = size_l
    result_set['mu'] = mu_l

    df = pd.DataFrame(data = result_set)
    df.to_csv('/Users/macbookpro/Project/size-csa/data/benchmarks/benchmark_500_0_1_result.csv',index = False)