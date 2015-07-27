#!/usr/bin/env python

import os
import numpy as np
from networkx.readwrite import json_graph
import json
from tuzcc import tuzcc

def load_graph(case_file):
    G = tuzc.MUGraph()
    with open(case_file, 'r') as f:
        gstring = json.load(f)
    Gs = json_graph.node_link_graph(gstring)
    node_num = Gs.number_of_nodes()
    s1 = node_num - 4
    s2 = node_num - 3
    t1 = node_num - 2
    t2 = node_num - 1
    edges = Gs.edges()
    for u, v in edges:
        par_num = int(Gs[u][v]['capacity'])
        for i in range(par_num):
            G.add_edge(u, v)

    G.set_sources([s1, s2])
    G.set_destinations([t1, t2])
    G.set_indices()
    return G

def run_all(G):
    coding_soln = G.get_coding_solution()
    gns, cuts = G.get_cut_sets()
    routing_soln = G.get_routing_solution()
    return coding_soln, routing_soln, gns, cuts

def print_results(G):
    c, r, g, d = run_all(G)
    print g
    print d
    print c.summary
    print r

def verify_results(dir):
    files = os.listdir(dir)
    files = [f for f in files if f.endswith('.json')]
    files = [os.path.join(os.path.dirname(dir+'/'), f) for f in files]
    for f in files:
        print f
        G = load_graph(f)
        while True:
            print_results(G)
            c = raw_input('Continue?')
            if c == 'y':
                break

def random_graph_coding_test(test_name, node_num, edge_num, case_num=1000, degrees=[None, None]):
    case_idx = {}
    for i in range(6):
        case_idx[i] = []

    for i in range(case_num):
        print "test number: ", i
        G = tuzc.MUGraph()
        G.set_random_gn_graph(node_num, num_edges=edge_num, degree_s=degrees)
        coding_soln, routing_soln, gns, cuts = run_all(G)

        gns_sum_rate_bound = np.min([gns[0], cuts[0] + cuts[1]])
        r_H1 = coding_soln.rks[0]
        r_G2 = coding_soln.rks[3]
        c11 = cuts[0]
        c22 = cuts[1]
        rr_1 = routing_soln[1]
        rr_2 = routing_soln[2]

        # check conditions
        # Compare with cuts
        # 0: sum rate < GNS (when GNS is tight)
        # 1: r_H1 < c11
        # 2: r_G2 < c22
        # Compare with routing
        # 3: sum rate < routing
        # 4: r_H1 < routing rate for 1
        # 5: r_G2 < routing rate for 2
        flags = [False] * 6
        if gns_sum_rate_bound == gns[0] and coding_soln.sum_rate_bound < gns[0]:
            flags[0] = True
            case_idx[0].append(i)
        if r_H1 < c11:
            flags[1] = True
            case_idx[1].append(i)
        if r_G2 < c22:
            flags[2] = True
            case_idx[2].append(i)
        if coding_soln.sum_rate_bound < routing_soln[0]:
            flags[3] = True
            case_idx[3].append(i)
        if r_H1 < rr_1:
            flags[4] = True
            case_idx[4].append(i)
        if r_G2 < rr_2:
            flags[5] = True
            case_idx[5].append(i)

        if any(flags):
            # write the graph into a json file
            json_str = json_graph.node_link_data(G.directed_simple_G)
            filename = './' + test_name + '/' + "case_" + str(i) + ".json"
            resultfilename  = './' + test_name + '/' + "case_" + str(i) + ".txt"

            with open(resultfilename, 'w') as f:
                f.write("++++++++++++ Instance " + str(i) + " ++++++++++++++\n")
                f.write("GNS Cut: " + str(gns[0]) + "\n")
                f.write("Src-Dst Cuts: " + str(cuts) + "\n")
                f.write("Routing Solution: " + str(routing_soln) + "\n")
                f.write(coding_soln.summary)

            with open(filename, 'w') as f:
                json.dump(json_str, f)
                #f.write(json_str)
                #f.write(flags)

    reportFile = './' + test_name + '/' + 'test_summary'
    scenarios = {}
    scenarios[0] = "sum rate < GNS (when GNS is tight) \n"
    scenarios[1] = "rank of H1 < c11 \n"
    scenarios[2] = "rank of G2 < c22 \n"
    scenarios[3] = "sum rate < routing sum rate \n"
    scenarios[4] = "rank of H1 < r11 for routing \n"
    scenarios[5] = "rank of G2 < r22 for routing \n"

    with open(reportFile, 'w') as f:
        for j in range(len(flags)):
            f.write(scenarios[j])
            f.write(str(case_idx[j]))
            f.write('\n')

if __name__ == "__main__":
    # prompt for arguments

    prompt = '> '
    print "Test name: "
    test_name = raw_input(prompt)

    print "number of nodes: "
    node_num = int(raw_input(prompt))

    print "number of edges (default is (number of node)^2 / 20): "
    edge_num = raw_input(prompt)
    if not edge_num:
        edge_num = int(node_num ** 2)
    else:
        edge_num = int(edge_num)

    print "source degree constraints type: "
    print "0 (default): no constraints."
    print "1: s1-t1 cut is 1"
    print "2: s2-t2 cut is 1"
    cut_type = int(raw_input(prompt))
    if cut_type == 1:
        degrees = [1, None]
    elif cut_type == 2:
        degrees = [None, 1]
    else:
        degrees = [None, None]

    print "number of random cases to be tested: "
    case_num = int(raw_input(prompt))
    test_dir = './' + test_name
    if not os.path.exists(test_dir):
        os.makedirs(test_dir)

    # if len(sys.argv) < 3:
    #     print "Usage: test test_name number_of_nodes degree_constraints number_of_case_runs number_of_edges"
    #     print "Default number of edge: (number of node)^2 / 20"
    #     print "Default number of cases: 1000"
    #     print "Default degree constraint: None"
    #     sys.exit(1)
    #
    #
    # test_name = sys.argv[1]
    # test_dir = './' + test_name
    # if not os.path.exists(test_dir):
    #     os.makedirs(test_dir)
    # node_num = int(sys.argv[2])
    # if len(sys.argv) >= 4:
    #     edge_num = int(sys.argv[4])
    # else:
    #     edge_num = int(node_num ** 2)
    #
    # if len(sys.argv) >= 5:
    #     case_num = int(sys.argv[4])
    # else:
    #     case_num = 100
    #
    # if len(sys.argv) >= 6:
    #     pass
    # else:
    #     degrees = [None, None]
    #
    random_graph_coding_test(test_name, node_num, edge_num, case_num, degrees)
