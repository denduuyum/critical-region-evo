"""
CA-evo.py:
A copy of CA_ES_DCNP_FINAL_idle.py.
Colletive Attacker based Evolutionary algorithm
for Distance Based Critical node problem.
A version of CA_ES_DCNP_FINAL after parameter tuning.
"""
import networkit as nk
from sys import stdin, stdout
import sys
import matplotlib.pyplot as plt
import numpy as np
import random as rd
import math
import copy
import UnionFind as UF
import graph_tools
import argparse
import os
from concurrent.futures import ThreadPoolExecutor
import threading
import time

VERSION = 'v1.0'
OUTPUT_DIR = "./output/DCNP"
time_exceeded = False

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    # print(*args, file=sys.stdout, **kwargs)


def H_function(G, nodes, uf):
    """
    @G  A graph from which the nodes in @nodes are removed
    @nodes  Deleted nodes
    @uf  Union find structure made out of @G
    
    Sorts the nodes in @nodes according to the following formula after restoring that node.
      H(u) = SUM_{component \in G} |component| * (|component| - 1) / 2
    It essentially sums quadratic forms of the component sizes in the graph.
    Union-Find is used for a fast calculation of the sums.
    """
    score = {}

    H = 0
    for sz in uf.sizes():
        H += sz * (sz - 1) / 2
    
    for u in nodes:
        t_test = [u]
        total = uf.size(u)
        sub = uf.size(u) * (uf.size(u) - 1) / 2
        for v in G.iterNeighbors(u):
            x = uf.find(v)
            if uf.same(u, x) == False and x not in t_test:
                t_test.append(x)
                total += uf.size(x)
                sub += uf.size(x) * (uf.size(x) - 1) / 2

        score[u] = H + total * (total - 1) / 2 - sub
    return score

def reinsertion(G, SG, n, cut_nodes, deleted_nodes, budget = None):
    """
    @G  Unmodified original graph
    @SG A subgraph of @G from which @cut_nodes are deleted
    @n Number of the nodes of @G (E.g. Maximum ID of the nodes)
    @cut_nodes contains the deleted nodes of @SG
    @deleted_nodes to which the final deleted nodes are to be appended
    @budget the budget of the problem
    
    When the size of @cut_nodes are more than the @budget, E = @budget - length(@cut_node) nodes need
    to be added back (restored). It reinserts the nodes 
    """

    nn = SG.numberOfNodes() + len(cut_nodes)
    uf_bk = UF.UnionFind(n)
    for e in SG.iterEdges():
        uf_bk.union(e[0], e[1])
            
    ssize = H_function(G, cut_nodes, uf_bk)
    # Get permutation of sorted cut nodes by their neighboring set sizes
    p = [i for i in range(len(cut_nodes))]
    p = sorted(p, key = lambda i: ssize[cut_nodes[i]])

    left = len(cut_nodes) - budget
    actual_reinserted = min(2*left, len(cut_nodes))
    # print('Left: ', left)
    reinserted = []
    for x in p[:actual_reinserted]:
        u = cut_nodes[x]
        reinserted.append(u)
        SG.restoreNode(u)
        for v in G.iterNeighbors(u):
            if SG.hasNode(v):
                SG.addEdge(u, v)
    bc = nk.centrality.Betweenness(SG)
    bc.run()
    bc_data = bc.scores()
    reinserted = sorted(reinserted, key = lambda x : bc_data[x], reverse = True)
    nodes = reinserted[:actual_reinserted - left]
    for u in nodes:
        SG.removeNode(u)
    for x in p[actual_reinserted:]:
        deleted_nodes.append(cut_nodes[x])
    for u in nodes:
        deleted_nodes.append(u)


def CA_extraction(G, b, c2, dist, logn_pow):
    """
    Formerly known find_topK
    Returns sqrt(b) nodes with the highest d-BC
    """
    bc = nk.centrality.EstimateBetweenness(G, logn_pow, _K = dist)
    bc.run()
    bc_data = bc.scores()

    topK = max(1, int(round(c2 * math.sqrt(b), 0)))
    rnodes = []
    p = [u for u in G.iterNodes()]
    p = sorted(p, key = lambda u : bc_data[u], reverse = True)
    for i in range(0, topK):
        if i >= len(p):
            break
        rnodes.append(p[i])

    return rnodes

def elite_gene_extraction(G, n, nlog, budget, dist, T, c2, CAs = None):
    BG = nk.Graph(G)
    s_clock = time.perf_counter()
    deleted_nodes = []
    bc_iter = 0
    cz_alpha = 0
    
    while len(deleted_nodes) < 2 * budget:
        if time.perf_counter() - s_clock > T:
            time_exceeded = True
            break
        eprint()
        eprint('Iteration: ', bc_iter + 1, ' -------------------------------------------------------')
        bc_iter += 1
        eprint('Number of deleted nodes:', len(deleted_nodes))

        logn_pow = 1
        t = math.log2(n - len(deleted_nodes))
        for i in range(nlog):
            logn_pow *= t
        logn_pow = int(logn_pow)


        
        rnodes = CA_extraction(G, budget, c2, dist, logn_pow)
        CA = []
        for v in rnodes:
            if G.hasNode(v):
                G.removeNode(v)
                CA.append(v)

        deleted_nodes += CA
        if CAs != None:
            CAs.append(CA)

    return deleted_nodes

def initialization(G, n, nlog, npopulation, CL, budget, dist, T, c2):
    CAs = []
    elite_gene_extraction(nk.Graph(G), n, 2, budget, dist, T, c2, CAs)
    P = []
    for i in range(npopulation):
        P.append(copy.deepcopy(rd.choices(CAs, k = CL)))

    return P, CAs

def recombination(P, CAs, CL, BG):
    offspring = []
    n_preserved = CL
    for i in range(len(P)):
        selected = rd.choices(P, k = 2)
        genes1 = copy.deepcopy(rd.choices(selected[0], k = math.floor(n_preserved / 2)))
        genes2 = copy.deepcopy(rd.choices(selected[1], k = math.floor(n_preserved / 2)))
        offspring.append(genes1 + genes2)
    return offspring

def repair(BG, n, P, CL, budget, dist):
    G = nk.Graph(BG)
    deleted_nodes = set()
    for CA in P:
        graph_tools.remove_nodes(G, CA)
        for u in CA:
            deleted_nodes.add(u)
    deleted_nodes = list(deleted_nodes)
    #eprint('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx: ', G.numberOfNodes())

    # print('# del nodes: ', len(deleted_nodes), end = ' ')
    if len(deleted_nodes) > budget:
        final_deleted_nodes = []
        reinsertion(nk.Graph(BG), G, n, deleted_nodes, final_deleted_nodes, budget = budget)
    else:
        final_deleted_nodes = deleted_nodes


    for CA in P:
        for x in list(CA):
            if x not in final_deleted_nodes:
                CA.remove(x)

    if len(final_deleted_nodes) < budget:
        left = budget - len(final_deleted_nodes)
        bc = nk.centrality.Betweenness(G)
        bc.run()
        bc_data = bc.scores()
        nodes = [u for u in G.iterNodes()]
        nodes = sorted(nodes, key = lambda x: bc_data[x], reverse = True)
        CA = rd.sample(nodes[:2*left], k = left)
        P.append(CA)
        final_deleted_nodes += CA

def objective_function(G, D, dist):
    for u in D:
        if G.hasNode(u):
            G.removeNode(u)

    D = 0
    for u in G.iterNodes():
        bfs = nk.distance.BFS(G, u, storePaths=False, storeNodesSortedByDistance=True, _K = dist)
        bfs.run()
        nodes = bfs.getNodesSortedByDistance()
        D += len(nodes) - 1
    
    return D / 2
    
def EA(G, n, tmax, npopulation, budget, dist, c1, c2, idle_max):
    # IDLE_ITER = max(100, 10 * budget)
    # IDLE_ITER = 2000
    s_clock = time.perf_counter()
    CL = max(int(c1 * math.sqrt(budget)*math.log(n,10)), 5)
    # print('dna_size: ', dna_size)
    # print('gene_size: ', int(round(math.sqrt(budget), 0)))

    P, CAs = initialization(G, n, 2, npopulation, CL, budget, dist, tmax, c2)
    P_del_nodes = []
    P_obj = []
    
    i_cnt = 0
    best = (n + 1) * (n + 1) * (n + 1)
    best_del_nodes = None
    best_time = 0

    for person in P:
        repair(nk.Graph(G), n, person, CL, budget, dist)
        del_nodes = set()
        for x in person:
            for y in x:
                del_nodes.add(y)

        P_del_nodes.append(del_nodes)
        obj = objective_function(nk.Graph(G), del_nodes, dist)
        P_obj.append(obj)
        if obj < best:
            best = obj
            best_del_nodes = list(del_nodes)
            best_time = time.perf_counter() - s_clock
            #print('ooooooooooooooooooooooooooooooooooooooooooooooooooo Best:', time.perf_counter() - s_clock, i_cnt, len(set(del_nodes)), best)

    i_idle_cnt = 0
    psz = [i for i in range(len(P))]
    psz = sorted(psz, key = lambda x : P_obj[x], reverse = True)
        
    while time.perf_counter() - s_clock < tmax and i_idle_cnt < idle_max:
        i_cnt += 1
        i_idle_cnt += 1
        # eprint('########################### iteration ' + str(i_cnt) + ' #################################################################')
        S = recombination(P, CAs, CL, nk.Graph(G))
        for offspring in S:
            repair(nk.Graph(G), n, offspring, CL, budget, dist)
            del_nodes = set()
            for x in offspring:
                for y in x:
                    del_nodes.add(y)

            obj = objective_function(nk.Graph(G), del_nodes, dist)
            if obj < best:
                i_idle_cnt = 0
                best = obj
                best_del_nodes = list(del_nodes)
                best_time = time.perf_counter() - s_clock
                #print('ooooooooooooooooooooooooooooooooooooooooooooooooooo Best:', best_time, i_cnt, len(set(del_nodes)), best)

            for i in range(len(P)):
                if P_obj[psz[i]] > obj:
                    if i > 0:
                        P_obj[psz[i - 1]] = P_obj[psz[i]]
                        P[psz[i - 1]] = P[psz[i]]
                    
                    P_obj[psz[i]] = obj
                    P[psz[i]] = offspring
                else:
                    break

    return best, best_del_nodes, best_time, time.perf_counter() - s_clock,

def arg_setup():
    parser = argparse.ArgumentParser(description='Betweenness separator.')
    parser.add_argument('-T', '--tmax', dest='tmax', type = int, 
                    help='The run time of the algorithm', default = 3600)
    parser.add_argument('-N', '--npopulation', dest='npopulation', type = int, 
                    help='The run time of the algorithm', default = 3)
    parser.add_argument('-n', '--nruns', dest='runs', type = int, 
                    help='Number of different runs', default = 1)
    parser.add_argument('-b', '--budget', dest='budget', type = int, required = True,
                    help='Number of nodes to delete')
    parser.add_argument('-D', '--dist', dest='dist', type = int, 
                    help='Distance parameter', default = 3)

    parser.add_argument('-c1', '--coef', dest='c1', type = float, 
                    help='Coefficient before chromosome length (CL) (dna size)', default  = 2)

    parser.add_argument('-i', '--idle_max', dest='idle_max', type = int, 
                    help='number of idles size', default=100)

    parser.add_argument('-c2', '--gscoef', dest='c2', type = float, 
                    help='Coefficient before gene size', default = 1)


    
    parser.add_argument('file', help = 'The file that contains graph information')
    return parser

def main(): 
    # Let's run in 4 threads
    nk.setNumberOfThreads(1)
    
    parser = arg_setup()
    args = parser.parse_args()

    with open(args.file, 'r') as fin:
        G, n = graph_tools.read_graph_mapped(fin)
        
    graph_tools.clean_graph(G)

    BG = nk.Graph(G)


    fname = 'CA_ES_DCNP_FINAL_dna_size'
    fname += str(args.tmax) + '_' + str(args.npopulation) + '_' + str(int(5 * math.sqrt(n))) + '_' + str(args.budget) + '_' + str(args.dist) + '_'
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fname = OUTPUT_DIR + '/' + fname + args.file.split('/')[-1]  + '.out'
    wrst = best = None
    best_del_nodes = None
    best_time = None
    best_total_time = None
    avg = 0
    avg_time = 0
    avg_total_time = 0
    
    with open(fname, 'w') as f:
        f.write('# obj value, best_time, total_time: deleted nodes\n')
        for i in range(args.runs):
            c_best, c_best_del_nodes, c_best_time, c_total_time = EA(G, n, args.tmax, args.npopulation, args.budget, args.dist, args.c1, args.c2, args.idle_max)

            if time_exceeded:
                f.write("Time exceeded: Halts now\n")
                break
            f.write(str(c_best) + ', ' + str(c_best_time) + ', ' + str(c_total_time) + ': ' + str(c_best_del_nodes) + '\n')
            f.flush()
            if best == None:
                best, best_del_nodes, best_time, best_total_time = c_best, c_best_del_nodes, c_best_time, c_total_time
                wrst = c_best
            elif best > c_best:
                best, best_del_nodes, best_time, best_total_time = c_best, c_best_del_nodes, c_best_time, c_total_time
            elif wrst < c_best:
                wrst = c_best
                
            avg += c_best
            avg_time += c_best_time
            avg_total_time += c_total_time
            G = nk.Graph(BG)

        if not time_exceeded:
            avg /= args.runs
            avg_time /= args.runs
            avg_total_time /= args.runs
            f.write(str(best) + ' ' + str(avg) + ' ' + str(wrst))
            f.write(str(round(best_time, 2)) + ' ' + str(round(best_total_time, 2)))
            f.write(str(round(avg_time, 2)) + ' ' + str(round(avg_total_time, 2)))
            f.flush()
    if not time_exceeded:        
        print('----------------------------------------------------------')
        print('CA-evo is executed with the following parameters.')
        print('file= ',args.file, ' nrun=',args.runs, 'T=', args.tmax, ' N=',args.npopulation, ' Chromosome lenght=', 
              max(int(args.c1*math.sqrt(args.budget)*math.log(n,10)), 5), 
              'gene size=', int(args.c2*round(math.sqrt(args.budget), 0)), 
              'budget=',args.budget,
              'c1=', args.c1,
              'c2=', args.c2,
              'idle_max=', args.idle_max)
        print('----------------------------------------------------------')
        print(int(best), round(avg,1), int(wrst), end = ' ')
        #print(round(best_time, 2), round(best_total_time, 2), end = ' ')
        #print(round(avg_time, 2), round(avg_total_time, 2))
        print(round(avg_total_time, 2))
        print('----------------------------------------------------------')


if __name__ == '__main__':
    sys.setrecursionlimit(int(1e9))
    main()
