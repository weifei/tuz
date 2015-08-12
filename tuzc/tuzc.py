import numpy as np
import networkx as nx
import random
import align
import pandas as pd
from gurobipy import *
from util import *
from collections import OrderedDict


class MUGraph(nx.MultiDiGraph):

    def __init__(self, data=None, **attr):
        super(MUGraph, self).__init__(data, **attr)
        self.dsts_evolution = None
        self.coding_matrix = None
        self.directed_simple_G = None
        self.sources = None
        self.destinations = None
        self.alignment_nodes = {}
        self.ordered_edges = None
        if data is not None:
            self.set_indices()

    def set_sources(self, srcs):
        self.sources = srcs

    def set_destinations(self, dsts):
        self.destinations = dsts

    def set_indices(self):
        if not nx.is_directed_acyclic_graph(self):
            raise ValueError('The graph is not DAG')
        if not nx.is_connected(self.to_undirected()):
            raise ValueError('The graph is not connected')

        self.directed_simple_G = nx.DiGraph(self)
        self.ordered_nodes = nx.topological_sort(self)
        for idx, node in enumerate(self.ordered_nodes):
            self.node[node]['index'] = idx

        self.ordered_edges = OrderedDict({})

        index = 0
        for tail in self.ordered_nodes:
            for head in sorted(self[tail]):
                self.ordered_edges[(tail, head)] = index
                self.directed_simple_G[tail][head]['capacity'] = len(self[tail][head])
                for idx in sorted(self[tail][head]):
                    self[tail][head][idx]['index'] = index
                    index = index + 1

        # reset data structures
        self.coding_matrix = None
        self.dsts_evolution = None
        self.alignment_nodes = []


    def get_edge_index(self, edge):
        tail = edge[0]
        head = edge[1]

        if len(edge) == 3:
            idx = edge[2]
            if idx >= self.number_of_edges(tail, head):
                raise ValueError('Edge ' + str(edge) + 'does not exist!')
            return self.ordered_edges[(tail, head)] + idx

        return sorted([self[tail][head][ix]['index'] for ix in self[tail][head].keys()])

    def get_edge_from_index(self, idx):
        for edge, start in reversed(self.ordered_edges.items()):
            if start <= idx:
                return (edge[0], edge[1], idx-start)

    def get_edges_indices(self, edges):
        index = []

        for edge in edges:
            tail = edge[0]
            head = edge[1]
            if len(edge) == 3 and edge[3] < self.number_of_edges(tail, head):
                idx = [edge[2]]
            else:
                idx = self[tail][head].keys()
            index.extend(sorted([self[tail][head][ix]['index'] for ix in idx]))

        return index

    def get_vertex_index(self, vertex):
        return self.node[vertex]['index']

    def get_vertex_from_index(self, idx):
        return self.ordered_nodes[idx]

    def set_coding_field(self, F):
        if F == 'Real':
            zeroElement=0.0
            identityElement=1.0
            add=operator.__add__
            sub=operator.__sub__
            mul=operator.__mul__
            div = operator.__div__
            eq = operator.__eq__
            str=lambda x:`x`,
            equalsZero = None
        else:
            F = F

    def get_all_coding_matrices(self):
        """
            this is to be called after getting the coding solution,

            it returns all the coding matrices for each step and their rank
        """
        class tran_matrix():
            def __init__(self):
                self.S = None
                self.T = None
                self.H1 = None
                self.H2 = None
                self.G1 = None
                self.G2 = None
                self.H = None
                self.H2G2 = None

            def __init__(self, S, T, H1, H2, G1, G2, H, H2G2):
                self.S = S
                self.T = T
                self.H1 = H1
                self.H2 = H2
                self.G1 = G1
                self.G2 = G2
                self.H = H
                self.H2G2 = H2G2

        s = self.sources
        s_edges = [self.out_edges(x) for x in s]
        S = [list(set(self.get_edges_indices(src))) for src in s_edges]

        matrices = []
        for dsts in self.dsts_evolution[::-1]:
            T = [sorted(x) for x in dsts]
            H2_raw = self.coding_matrix[T[0]].ix[S[1]].values
            G2_raw = self.coding_matrix[T[1]].ix[S[1]].values
            H1 = numpy_to_generic(self.coding_matrix[T[0]].ix[S[0]].values)
            H2 = numpy_to_generic(H2_raw)
            G1 = numpy_to_generic(self.coding_matrix[T[1]].ix[S[0]].values)
            G2 = numpy_to_generic(G2_raw)
            H = numpy_to_generic(self.coding_matrix[T[0]].values)
            H2G2 = numpy_to_generic(np.hstack((H2_raw, G2_raw)))
            temp = tran_matrix(S, T, H1, H2, G1, G2, H, H2G2)
            matrices.append(temp)

        return matrices

    def get_coding_solution(self, F='Real'):
        """
            F: the underlying coding field
        """
        self.set_coding_field(F)
        self.run_destination_reduction()
        self.run_recursive_coding()
        return self.get_solution_summary()

    def run_destination_reduction(self):
        """
            instance method fro destination reduction algorithm

            self: The network itself, including specifications of
                - 2 Sources, obtained from G.graph['srcs']
                - 2 Destiantions, obtained from G.graph['dsts']

            F: the underlying field for operation

            Return: a structure that specifies the reduction procedure

        """
        s = self.sources
        t = self.destinations

        # get the union set of source edges
        s_edges = [self.out_edges(x) for x in s]
        s_edges_idx = [set(self.get_edges_indices(src)) for src in s_edges]
        S = set.union(*s_edges_idx)

        # now starting from the sources
        T_curr = [set(self.in_edges(x)) for x in t]
        T_curr_idx = [set(self.get_edges_indices(x)) for x in T_curr]
        T_union = set.union(*[set(x) for x in T_curr])
        T_union_idx = set.union(*T_curr_idx)
        self.dsts_evolution = [T_curr_idx]

        while not T_union_idx.issubset(S):
            # get the edges out of the max-ordered node
            tail_nodes = set([x[0] for x in T_union])
            tail_nodes = tail_nodes - set(s)
            v = max(tail_nodes, key=lambda p: self.node[p]['index'])

            E = set([edge for edge in T_union if edge[0] == v])
            E_curr = [x.intersection(E) for x in T_curr]

            T_next = []
            Uv = self.in_edges(v)
            for (Ex, Tx) in zip(E_curr, T_curr):
                if Ex:
                    T_next.append((Tx - Ex).union(Uv))
                else:
                    T_next.append(Tx)

            T_curr = T_next
            T_curr_idx = [set(self.get_edges_indices(x)) for x in T_curr]
            self.dsts_evolution.append(T_curr_idx)

            T_union = set.union(*[set(x) for x in T_curr])
            T_union_idx = set.union(*T_curr_idx)


    def run_recursive_coding(self):
        """
            self: The network itself, including specifications of
                - 2 Sources, obtained from self.graph['srcs']
                - 2 Destiantions, obtained from self.graph['dsts']

            dsts_iters: Destations over iterations, the list of destinations
                created over the iterations of the destination reduction algorithm

            Return: A linear coding matrix from sources to destinations specified in G

        """
        # get the source edges (for now we assume that they equal to the min-cuts)
        s = self.sources
        t = self.destinations
        self.alignment_nodes = {}
        dsts_iter = self.dsts_evolution

        # get the union set of source edges
        s_edges = [self.out_edges(x) for x in s]
        s_edges_idx = [list(set(self.get_edges_indices(src))) for src in s_edges]

        S = []
        for x in s_edges_idx:
            S = S + list(x)

        edge_num = len(self.edges())

        # create the empty dataframe that representation the ENTIRE solution
        # zero_array = np.zeros((len(S), edge_num))
        # code_matrix = pd.DataFrame(zero_array, index=list(S), columns=range(edge_num))
        code_matrix = pd.DataFrame(index=list(S), columns=range(edge_num), dtype='int')
        for col in code_matrix.columns:
            code_matrix[col] = [int(0)] * len(code_matrix.index)

        # initialize identity transfer matrix for the source edges
        for edge in S:
            code_matrix[edge][edge] = F.unity

        dsts_evo = dsts_iter[::-1]
        T_curr = dsts_evo[0]

        for dst in dsts_evo[1:]:
            # print code_matrix
            # print dst
            # load the next destination sets
            T_next = dst

            # get groups of edges for easy indexing
            O1 = T_next[0] - T_curr[0]
            O2 = T_next[1] - T_curr[1]
            I1 = list(T_curr[0] - T_next[0])
            I2 = list(T_curr[1] - T_next[1])
            I = list(set(I1).union(set(I2)))
            O = list(O1)
            B = list(O2 - O1)
            U1 = list(T_next[0].intersection(T_curr[0]))
            U2 = list(T_next[1].intersection(T_curr[1]))

            if len(I) == 0:
                T_curr = T_next
                continue

            # get incoming edge matrix
            # convert it to genericmatrix
            #print "input edges: ", I
            in_matrix = code_matrix[I].values

            # Step 1: code for B
            if B:
                # print "Random Coding for Edges for T2 only edges", B
                code_matrix[B] = generate_rc_matrix(in_matrix, len(B))

            if not O:
                T_curr = T_next
                continue

            # Step 2
            # code for O
            # edge by edge coding for alignment

            while O:
                # pick a random edge from O
                # curr_o = random.choice(O)
                if align.check_alignment_criterion(code_matrix, U1, U2, I1,
                                                   I2, B, s_edges_idx):
                    curr_o = max(O)
                    # print "alignment on : ", curr_o
                    temp = self.get_edge_from_index(I[0])
                    if temp[1] in self.alignment_nodes:
                        self.alignment_nodes[temp[1]].append(curr_o)
                    else:
                        self.alignment_nodes[temp[1]] = [curr_o]
                    # code_matrix[O] = align.get_alignment_solution(
                    # code_matrix, s_edges_idx, U1, I1, len(O))
                    code_matrix[curr_o] = align.get_alignment_solution(
                        code_matrix, s_edges_idx, U1, I1, 1)
                else:
                    curr_o = min(O)
                    # print "random coding edges: ", curr_o
                    code_matrix[curr_o] = generate_rc_matrix(in_matrix, 1)
                    # add the coded edge into U1 and/or U2

                if curr_o in O1:
                    U1.append(curr_o)
                    U1.sort()
                if curr_o in O2:
                    U2.append(curr_o)
                    U2.sort()
                # remove the edge from O
                O = [x for x in O if x != curr_o]
#
            T_curr = T_next

        self.coding_matrix = code_matrix

    def get_cut_sets(self):
        c11 = self.get_min_cut(self.sources[0], self.destinations[0])
        c22 = self.get_min_cut(self.sources[1], self.destinations[1])
        c21 = self.get_min_cut(self.sources[1], self.destinations[0])
        gns = self.get_GNS_cut_reduced()
        cuts = [c11, c22, c21]
        return (gns, cuts)

    def get_min_cut(self, s, t):
        # get min cut between two nodes
        return nx.minimum_cut(self.directed_simple_G, s, t)

    def get_final_grank(self):
        return self.get_grank(self.dsts_evolution[0])

    def get_grank(self, T):
        s = self.sources

        # get the union set of source edges
        s_edges = [self.out_edges(x) for x in s]
        s_edges_idx = [list(set(self.get_edges_indices(src))) for src in s_edges]

        T = [list(t) for t in T]
        return grank(self.coding_matrix, s_edges_idx, T)

    def get_solution_summary(self):
        s = self.sources

        # get the union set of source edges
        s_edges = [self.out_edges(x) for x in s]
        s_edges = [list(set(self.get_edges_indices(src))) for src in s_edges]

        t = [list(x) for x in self.dsts_evolution[0]]
        H2_raw = self.coding_matrix[t[0]].ix[s_edges[1]].values
        G2_raw = self.coding_matrix[t[1]].ix[s_edges[1]].values

        H1 = numpy_to_generic(self.coding_matrix[t[0]].ix[s_edges[0]].values)
        H2 = numpy_to_generic(H2_raw)
        G1 = numpy_to_generic(self.coding_matrix[t[1]].ix[s_edges[0]].values)
        G2 = numpy_to_generic(G2_raw)
        r_H1 = H1.Rank()
        r_H2 = H2.Rank()
        r_G1 = G1.Rank()
        r_G2 = G2.Rank()

        rks = [r_H1, r_H2, r_G1, r_G2]

        mats = [H1, H2, G1, G2]

        H = numpy_to_generic(self.coding_matrix[t[0]].values)
        H2G2 = numpy_to_generic(np.hstack((H2_raw, G2_raw)))
        g1 = H.Rank()
        g2 = H2G2.Rank()
        g3 = r_H2
        final_grank = g1 + g2 - g3

        sum_rate_bound = min(final_grank, r_H1 + r_G2)

        gmats = [H, H2G2, H2]
        grks = [g1, g2, g3]

        summary = "Neutralization happened at nodes: " + \
                  str(self.alignment_nodes) + '\n'
        summary = summary + "Final Grank = " + str(final_grank) + '\n'
        summary = summary + "Sum Rate bound = " + str(sum_rate_bound) + '\n'
        summary = summary + "Rank H1= " + str(r_H1) + '\n'
        summary = summary + "Rank G2= " + str(r_G2) + '\n'
        summary = summary + "Rank H2= " + str(r_H2)

        class soln_summary():
            def __init__(self):
                self.grank = None
                self.sum_rate_bound = None
                self.gmats = None
                self.grks = None
                self.mats = None
                self.rks = None
                self.summary = None

            def __init__(self, grank, sum_rate, gmats, grks, mats, rks, summary):
                self.grank = grank
                self.sum_rate_bound = sum_rate
                self.gmats = gmats
                self.grks = grks
                self.mats = mats
                self.rks = rks
                self.summary = summary

        summ_struct = soln_summary(final_grank, sum_rate_bound, gmats, grks, mats, rks, summary)

        return summ_struct

    def set_random_gn_graph(self, num_nodes, num_edges=None, degree_s=[None, None]):
        """
            Set the graph to a "layered" random dag
        """
        # remove the current graph first
        self.clear()

        if num_edges is None:
            num_edges = int((num_nodes ** 2 / 4))
        #
        # first create a GN graph

        #G = nx.gn_graph(num_nodes)
        G = nx.gn_graph(num_nodes)
        H = nx.DiGraph()
        for u, v in G.edges():
            H.add_edge(v, u, weight=1)
        G = H
        for u, v in G.edges():
            G[u][v]['weight'] = 1

        nodes = nx.topological_sort(G)
        num_edges = num_edges - G.number_of_edges()
        for i in range(num_edges):
            u_idx = random.choice(range(len(nodes)-1))
            u = nodes[u_idx]
            v = random.choice(nodes[u_idx+1:])
            if (u,v) in G.edges():
                G[u][v]['weight'] += 1
            else:
                G.add_edge(u, v, weight=1)

        self.set_random_session(G, degree_s)


    def set_random_session(self, G, degree_s):
        """
            Get a base graph and pick a random 2-unicast session on the graph
            choose degrees properly if needed
        """
        sorted_nodes = nx.topological_sort(G)
        num_nodes = G.number_of_nodes()

        # create sources and destinations of each of the sections
        # name the nodes to be the last 4 numbers
        srcs = [num_nodes, num_nodes + 1]
        dsts = [num_nodes + 2, num_nodes + 3]

        end_idx = int(0.3 * len(sorted_nodes))
        end_idx = max(end_idx, 2)
        for i in range(2):
            s = srcs[i]
            t = dsts[i]
            reachables = []
            iter_num = 0

            while len(reachables) == 0:
                iter_num += 1
                if iter_num > 100:
                    end_idx = end_idx * 2

                # pick an entry point from the first 30%
                entry_point = random.choice(sorted_nodes[:end_idx])
                # print "Source ", i
                # print "candidates: ", sorted_nodes[:end_idx]
                # print "entry point: ", entry_point
                # print "all nodes: ", G.nodes()

                # pick a random point from the reachables
                reachables = nx.shortest_path(G, entry_point)
                del reachables[entry_point]
                #print "reachables: ", reachables
                reachables = reachables.keys()

            exit_point = random.choice(reachables)
            #print "exit_point: ", exit_point

            if degree_s[i]:
                G.add_edge(s, entry_point, weight=degree_s[i])
                G.add_edge(exit_point, t, weight=degree_s[i])
            else:
                # figure out the out_degree of entry point
                out_degree = np.sum(G[u][v]['weight'] for u,v in G.out_edges(entry_point))
                G.add_edge(s, entry_point, weight=out_degree)

                # figure out the int_degree of exit point
                in_degree = np.sum(G[u][v]['weight'] for u,v in G.in_edges(exit_point))
                G.add_edge(exit_point, t,  weight=in_degree)

        edges = G.edges()
        for u, v in edges:
            par_num = int(G[u][v]['weight'])
            for i in range(par_num):
                self.add_edge(u, v)

        # set indices etc
        self.set_sources(srcs)
        self.set_destinations(dsts)
        self.set_indices()
        #print "number of nodes: " + str(self.number_of_nodes())
        #print "number of edges: " + str(self.number_of_edges())


    def set_random_dag(self, num_nodes, num_edges=None, degree_s=[None, None]):
        """
            Make G a random DAG of n nodes
            if twounicast is True, set the sources and destination such that we have a (1,N)
            configuration.

            Do I want to fix min-cut out degree at the sources?

        """
        # first create a DAG with num_nodes number of nodes by randomly pack
        # adjacency matrix which is lower triangular
        # this may not be the best way
        # how to generate MultiDigraph?
        if num_edges is None:
            num_edges = int((num_nodes ** 2 / 4))

        while True:
            adj_matrix = np.zeros((num_nodes, num_nodes))

            positions = []
            for i in range(num_nodes):
                for j in range(i+1, num_nodes):
                    positions.append((i,j))

            # create random edges
            rand_idx = [random.randint(0, len(positions)-1) for x in range(num_edges)]
            for i in rand_idx:
                a, b = positions[i]
                adj_matrix[a][b] += 1

            # Done with adjacency matrix
            G = nx.from_numpy_matrix(adj_matrix, create_using=nx.DiGraph())

            if nx.is_connected(nx.Graph(G)):
                break

        self.set_random_session(G, degree_s)

    def get_GNS_cut(self):
        """
            Test Version, calculating GNS cut for two-unicast-z

            Using Gurobi integer optimization

        """
        # we build the optimization around the casted digraph instead of multidigraph
        # for simplicity
        G = self.directed_simple_G
        s_1 = self.sources[0]
        s_2 = self.sources[1]
        t_1 = self.destinations[0]
        t_2 = self.destinations[1]
        edges = G.edges()
        nodes = G.nodes()

        try:

            # Great an gurobi instance of the optimization model
            m = Model("GNS")
            m.setParam('OutputFlag', False)

            x_v = {}
            # vertex variables  for s_1, t_1 cut
            for v in nodes:
                x_v[v] = m.addVar(vtype=GRB.BINARY)

            x_e = {}
            # edge variables for s_1, t_1 cut
            for (u,v) in edges:
                x_e[u,v] = m.addVar(vtype=GRB.BINARY)

            y_v = {}
            # vertex variables  for s_2, t_2 cut
            for v in nodes:
                y_v[v] = m.addVar(vtype=GRB.BINARY)

            y_e = {}
            # edge variables for s_2, t_2 cut
            for (u,v) in edges:
                y_e[u,v] = m.addVar(vtype=GRB.BINARY)

            z_v = {}
            # vertex variables  for s_2, t_1 cut
            for v in nodes:
                z_v[v] = m.addVar(vtype=GRB.BINARY)

            z_e = {}
            # edge variables for s_2, t_1 cut
            for (u,v) in edges:
                z_e[u,v] = m.addVar(vtype=GRB.BINARY)

            e = {}
            # GNS indicator variable
            for (u,v) in edges:
                e[u,v] = m.addVar(vtype=GRB.BINARY, obj=G[u][v]['capacity'])

            # Done with decision variable creation
            # update model
            m.update()

            # Constraints
            # 1. Constraints for s_1 - t_1 cut
            for (u,v) in edges:
                if (u,v) == (s_1, t_1):
                    m.addConstr(x_e[u,v] >= 1)
                elif u == s_1:
                    m.addConstr(x_v[v] + x_e[u,v] >= 1)
                elif v == t_1:
                    m.addConstr(-x_v[u] + x_e[u,v] >= 0)
                else:
                    m.addConstr(x_v[v] - x_v[u] + x_e[u,v] >= 0)

            # 2. Constraints for s_2 - t_2 cut
            for (u,v) in edges:
                if (u,v) == (s_2, t_2):
                    m.addConstr(y_e[u,v] >= 1)
                elif u == s_2:
                    m.addConstr(y_v[v] + y_e[u,v] >= 1)
                elif v == t_2:
                    m.addConstr(-y_v[u] + y_e[u,v] >= 0)
                else:
                    m.addConstr(y_v[v] - y_v[u] + y_e[u,v] >= 0)

            # 3. Constraints for s_2 - t_1 cut
            for (u,v) in edges:
                if (u,v) == (s_2, t_1):
                    m.addConstr(z_e[u,v] >= 1)
                elif u == s_2:
                    m.addConstr(z_v[v] + z_e[u,v] >= 1)
                elif v == t_1:
                    m.addConstr(-z_v[u] + z_e[u,v] >= 0)
                else:
                    m.addConstr(z_v[v] - z_v[u] + z_e[u,v] >= 0)

            # 4. Constraints for e[u,v] >= max(x_e[u,v], y_e[u,v], z_e[u,v])
            for (u,v) in edges:
                m.addConstr(e[u,v] >= x_e[u,v])
                m.addConstr(e[u,v] >= y_e[u,v])
                m.addConstr(e[u,v] >= z_e[u,v])

            m.optimize()

            if m.status == GRB.status.OPTIMAL:
                print "Min GNS cut value = " + str(m.objVal)
                print "GNS cut edges:"

                for u,v in edges:
                    if e[u,v].x != 0:
                        print (u,v)
                print "s1-t1 cut edges in GNS:"
                for u,v in edges:
                    if x_e[u,v].x != 0:
                        print (u,v)

                print "s2-t2 cut edges in GNS:"
                for u,v in edges:
                    if y_e[u,v].x != 0:
                        print (u,v)

                print "s2-t1 cut edges in GNS:"
                for u,v in edges:
                    if z_e[u,v].x != 0:
                        print (u,v)
            else:
                # something went wrong...err...
                print "Something was wrong"

        except GurobiError:
            print ('Error report from Gurobi')

    def get_GNS_cut_reduced(self):
        """
            Reduced version of the previous one, testing
            calculating GNS cut for two-unicast-z

            Using Gurobi integer optimization

            returns: min GNS cut value, min GNS cut set edges (u,v,capacity)

        """
        # we build the optimization around the casted digraph instead of multidigraph
        # for simplicity
        G = self.directed_simple_G
        s_1 = self.sources[0]
        s_2 = self.sources[1]
        t_1 = self.destinations[0]
        t_2 = self.destinations[1]
        edges = G.edges()
        nodes = G.nodes()

        try:

            # Great an gurobi instance of the optimization model
            m = Model("GNS")
            m.setParam('OutputFlag', False)

            x_v = {}
            # vertex variables  for s_1, t_1 cut
            for v in nodes:
                x_v[v] = m.addVar(vtype=GRB.BINARY)

            y_v = {}
            # vertex variables  for s_2, t_2 cut
            for v in nodes:
                y_v[v] = m.addVar(vtype=GRB.BINARY)

            z_v = {}
            # vertex variables  for s_2, t_1 cut
            for v in nodes:
                z_v[v] = m.addVar(vtype=GRB.BINARY)

            e = {}
            # GNS indicator variable
            for (u,v) in edges:
                e[u,v] = m.addVar(vtype=GRB.BINARY, obj=G[u][v]['capacity'])

            # Done with decision variable creation
            # update model
            m.update()

            # Constraints
            # 1. Constraints for s_1 - t_1 cut
            for (u,v) in edges:
                if (u,v) == (s_1, t_1):
                    m.addConstr(e[u,v] >= 1)
                elif u == s_1:
                    m.addConstr(x_v[v] + e[u,v] >= 1)
                elif v == t_1:
                    m.addConstr(-x_v[u] + e[u,v] >= 0)
                else:
                    m.addConstr(x_v[v] - x_v[u] + e[u,v] >= 0)

                if (u,v) == (s_2, t_2):
                    m.addConstr(e[u,v] >= 1)
                elif u == s_2:
                    m.addConstr(y_v[v] + e[u,v] >= 1)
                elif v == t_2:
                    m.addConstr(-y_v[u] + e[u,v] >= 0)
                else:
                    m.addConstr(y_v[v] - y_v[u] + e[u,v] >= 0)

                if (u,v) == (s_2, t_1):
                    m.addConstr(e[u,v] >= 1)
                elif u == s_2:
                    m.addConstr(z_v[v] + e[u,v] >= 1)
                elif v == t_1:
                    m.addConstr(-z_v[u] + e[u,v] >= 0)
                else:
                    m.addConstr(z_v[v] - z_v[u] + e[u,v] >= 0)

            m.optimize()

            if m.status == GRB.status.OPTIMAL:
                #print "Min GNS cut value = " + str(m.objVal)
                #print "GNS cut edges:"
                cut_set_edges = []
                for u,v in edges:
                    if e[u,v].x != 0:
                        #print (u,v), str(G[u][v]['capacity'])
                        cut_set_edges.append((u,v, G[u][v]['capacity']))
                return (m.objVal, cut_set_edges)
            else:
                # something went wrong...err...
                print "Something was wrong"
                return None, None

        except GurobiError:
            print ('Error report from Gurobi')

    def get_routing_solution(self):
        """
            Compute the routing capacity of the underlying graph.

            Here we just use the underlying DiGraph, instead of the MultiDiGraph

            returns: optimal sum rate, rate r1, rate r2

        """
        G = self.directed_simple_G
        s1 = self.sources[0]
        s2 = self.sources[1]
        t1 = self.destinations[0]
        t2 = self.destinations[1]

        try:
            m = Model('routing')
            m.setParam('OutputFlag', False)

            # variables,
            # We have one variable per edge per session
            # e is the dict of dict for the variables
            e = {}
            r = {}
            for i in [1,2]:
                e[i] = {}
                r[i] = m.addVar()
                for u,v in G.edges():
                    e[i][u,v] = m.addVar(lb=0)

            m.update()

            obj = quicksum(r.values())
            m.setObjective(obj, GRB.MAXIMIZE)

            # constraints
            # 1. conservations of flow at all intermediate nodes
            # 2. capacity constraints for each edge

            for u,v in G.edges():
                m.addConstr(e[1][u,v] + e[2][u,v] <= G[u][v]['capacity'])

            m.addConstr(quicksum(e[1][u,v] for u,v in G.out_edges(s1)) == r[1])
            m.addConstr(quicksum(e[2][u,v] for u,v in G.out_edges(s2)) == r[2])
            m.addConstr(quicksum(e[1][u,v] for u,v in G.out_edges(s2)) == 0)
            m.addConstr(quicksum(e[2][u,v] for u,v in G.out_edges(s1)) == 0)
            m.addConstr(quicksum(e[1][u,v] for u,v in G.in_edges(t1)) == r[1])
            m.addConstr(quicksum(e[2][u,v] for u,v in G.in_edges(t2)) == r[2])

            for n in G.nodes():
                if n not in [s1, s2, t1, t2]:
                    for i in [1, 2]:
                        inflow = quicksum(e[i][u,v] for u,v in G.in_edges(n))
                        outflow = quicksum(e[i][u,v] for u,v in G.out_edges(n))
                        m.addConstr(inflow == outflow)

            m.optimize()

            if m.status == GRB.status.OPTIMAL:
                for u, v in G.edges():
                    G[u][v]['Routing'] = {}
                    G[u][v]['Routing'][1] = e[1][u,v].x
                    G[u][v]['Routing'][2] = e[2][u,v].x
                return (m.objVal, r[1].x, r[2].x)
            else:
                # something went wrong...err...
                print "Something was wrong, no optimal solution obtained"
                return None, None, None

        except GurobiError:
            Print ('Error Report from Gurobi')
