"""
Modified version of CR-greedy.py with middle reinsertion and estimate dBC. 
Critical region: sqrt(n) node with highest BC score and re-insert some of them.
Greedy Algorithm: sequentially removes CAs until meeting the constraint

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
import time
import array as arr
from networkit import components 

OUTPUT_DIR = "./output/"

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

def restore_node(G, SG, u):
    SG.restoreNode(u)
    for v in G.iterNeighbors(u):
        if SG.hasNode(v):
            SG.addEdge(u, v)

def repair_operator(G, SG, n, solution, dist, budget):

    # the last CR is extracted by CR-extraction that uses approximate dBC
    CR_size = max(int(math.sqrt(budget)), 1) 
    CR = CR_extraction_app(SG, n, CR_size)
    # To cummulate CRs
    solution += CR

    num_2restore = len(solution) - budget + 1

    eprint('----------------repair operator-----------------')
    eprint('len of S', len(solution))
    eprint('num to restore: ', num_2restore)    

    # computing H-function value for nodes in solution
    uf_bk = UF.UnionFind(n)
    for e in SG.iterEdges():
        uf_bk.union(e[0], e[1])
    ssize = H_function(G, solution, uf_bk)
    # Get permutation of sorted cut nodes by their neighboring set sizes
    p = [i for i in range(len(solution))]
    p = sorted(p, key = lambda i: ssize[solution[i]])    
    repaired_solution = []
    for x in p[:num_2restore]:
        u = solution[x]
        SG.restoreNode(u)
        for v in G.iterNeighbors(u):
            if SG.hasNode(v):
                SG.addEdge(u, v)
    
    for x in p[num_2restore:]:
        u = solution[x]
        repaired_solution.append(u)

    # the final extract 
    CR = CR_extraction_est(SG, n, 1)
    # To cummulate CRs
    repaired_solution += CR


    eprint('Len of deleted nodes before final re: ', len(solution))
    eprint('Len of deleted nodes after final re: ', len(repaired_solution))

    return repaired_solution
    
def print_graphs(G,n):
    # for u in G.iterNodes():
    #     print(u,':')
    for x in range(0, n):
        print(x,':',end=' ')
        for v in G.iterNeighbors(x):
            print(v, end=' ')
        print()

def print_component_info(G):
    cc = components.ConnectedComponents(G)
    cc.run()
    compSizes = cc.getComponentSizes()
    numCC = len(compSizes)
    maxCC = max(compSizes.values())
    # eprint("#cc = %d,largest = %d"%(numCC,maxCC))
    
    p = [i for i in range(len(compSizes))]
    p = sorted(p, key = lambda i: compSizes[i], reverse = True)
    eprint('\n\n\n***Largest three components in the remaining networks with estimate bc')
    r = min(numCC, 5)
    for x in range(0,r):
        eprint(compSizes[p[x]], end =' ')
    eprint('\n\n')

def CR_extraction_est(G, n, k):
    
    num_sample = 25*(round(math.log(n,2),0)**2)

    bc = nk.centrality.EstimateBetweenness(G, num_sample, parallel=True, _K = 3)
    bc.run()
    bc_data = bc.scores()

    CR_size = k
    CR = []

    p = [u for u in G.iterNodes()]
    p = sorted(p, key = lambda u : bc_data[u], reverse = True)
    
    for i in range(0, CR_size):
        if i >= len(p):
            break
        CR.append(p[i])

    # CR = CR_extraction(G, n, budget, bc_data)
    if len(CR) == 0:
        eprint('CR-extract is failed')
        exit()


    # removal of CR from G
    graph_tools.remove_nodes(G, CR)

    return CR

def CR_extraction_app(G, n, k):
    
    bc = nk.centrality.ApproxBetweenness(G, epsilon=0.04, delta = 0.1)
    bc.run()
    bc_data = bc.scores()

    CR_size = k
    CR = []

    p = [u for u in G.iterNodes()]
    p = sorted(p, key = lambda u : bc_data[u], reverse = True)
    
    for i in range(0, CR_size):
        if i >= len(p):
            break
        CR.append(p[i])


    if len(CR) == 0:
        eprint('CR-extract is failed')
        exit()


    # removal of CR from G
    graph_tools.remove_nodes(G, CR)

    return CR

def repair_appbc(G, SG, n, solution, CR_size, budget):
    eprint('----------------repair operator-----------------')
    num_2restore = (len(solution) + CR_size) - budget
    
    eprint('len of S', len(solution))
    eprint('num2restore: ', num_2restore)


    # the last CR-extract

    # the last CR is extracted by CR-extraction that uses approximate dBC
    CR = CR_extraction_app(SG, n, CR_size)
    # To cummulate CRs
    solution += CR
    eprint('Number of nodes in solution:', len(solution))
    
    # num_2restore = min(num_2left + 1, len(solution))
    # nn = SG.numberOfNodes() + len(solution)
    uf_bk = UF.UnionFind(n)
    for e in SG.iterEdges():
        uf_bk.union(e[0], e[1])
    ssize = H_function(G, solution, uf_bk)

    # Get permutation of sorted cut nodes by their neighboring set sizes
    p = [i for i in range(len(solution))]
    p = sorted(p, key = lambda i: ssize[solution[i]])
    # eprint('reinserted: ', num_2restore)
    # print('Len of cut nodes: ', len(solution))

    solution2 = []
    reinserted = {}
    for x in p[:]:
        if x < num_2restore + 1:
            print(x, solution[x])
            u = solution[x]
            reinserted[u] = 1
            SG.restoreNode(u)
            for v in G.iterNeighbors(u):
                if SG.hasNode(v):
                    SG.addEdge(u, v)
        else:
            u = solution[x]
            solution2.append(u)

    CR = CR_extraction_est(G, n, 1)
    solution2 += CR
    
    eprint('Number of nodes in solution:', len(solution2), len(solution))
    # exit()
    return solution2

def CR_greedy(G, n, budget, dist, nsample):
    s_clock = time.perf_counter()
    #initialize variables
    BG = nk.Graph(G)
    cummulative_CR = []
    solution = []
    i_cnt = 0    
    CR_size = max(int(math.sqrt(budget)), 1)  

    # main loop for constructing a solution by sequential call of CR-extraction that uses estimate dBC
    while len(cummulative_CR) < budget - CR_size:
    
        eprint('\nIteration: ', i_cnt + 1, ' -------------------------------------------------------')
        i_cnt += 1

        s_clock_local = time.perf_counter()
        #Extraction of CR
        CR = CR_extraction_est(G, n, CR_size)
        cummulative_CR += CR
        eprint('Number of deleted nodes:', len(cummulative_CR), 'time:      ', round(time.perf_counter()-s_clock, 4))        


    # repair the constructed solution if infeasible
    if (len(cummulative_CR) < budget):
        solution_fnl = repair_operator(nk.Graph(BG), G, n, cummulative_CR, dist, budget = budget)
        eprint('Validation of solution (length of solution ==? budget) : ',len(solution_fnl),' ==? ', budget)
    elif (len(cummulative_CR) == budget):
        solution = copy.deepcopy(cummulative_CR)    
    else:
        print('Error in solution construction')
        exit()

    # print_component_info(G)

    e_clock = time.perf_counter()
    elapsed_time = e_clock - s_clock

    return solution_fnl, elapsed_time

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

def my_read_graph_mapped(fin):
    n = int(fin.readline())
    G = nk.Graph(n)
    d = {}
    id = 0
    
    while True:
        try:
            line = fin.readline()
            
        except:
            break

        line = line.split()
        if len(line) == 0:
            break

        x = int(line[0][:-1])
        id += 1

        arr = [int(y) for y in line[1:]]
        for y in arr:
            if (x < y):
                G.addEdge(x, y, addMissing = True)

    return G, id

def arg_setup():
    parser = argparse.ArgumentParser(description='Betweenness separator.')
    parser.add_argument('-n', '--nrun', dest='nrun', type = int, help='The number of runs', default = 1)
    parser.add_argument('-b', '--budget', dest='budget', type = int, required = True,
                    help='Number of nodes to delete')
    parser.add_argument('-D', '--dist', dest='dist', type = int, default = 3,
                    help='Distance parameter')
    parser.add_argument('file', help = 'The file that contains graph information')
    return parser

def main(): 

    # Let's run in t threads
    num_t = 1
    nk.setNumberOfThreads(num_t)

    parser = arg_setup()
    args = parser.parse_args()

    with open(args.file, 'r') as fin:
        G, n = my_read_graph_mapped(fin)
    
    num_sample = 25*(round(math.log(n,2),0)**2)
    BG = nk.Graph(G)
    print('file', args.file, ', nrun',args.nrun, ', B',args.budget, ', dist', args.dist, ', nsample', num_sample, ', number of threads', num_t)

    best_S = None
    best_obj = None
    best_time = None
    avg = 0
    etime_avg = 0
    # print('---------------------------------------estimate BC----------------------------------')
    # print(num_sample)
    # exit()
    for i in range(args.nrun):
        rd.seed(time.perf_counter()%rd.randint(1,10000))
        S, etime = CR_greedy(nk.Graph(G), n, budget = args.budget, dist = args.dist, nsample = num_sample)
        s_clock = time.perf_counter()
        obj = objective_function(nk.Graph(G), S, args.dist)
        obj = objective_function(nk.Graph(G), S, args.dist)
        e_clock = time.perf_counter()
        eprint('\n\n---time of OF computation with b = %d, time = %f:       ',(args.budget, e_clock-s_clock))
        eprint('---time of total CR-greedy computation:       ',etime,'\n\n')
        etime += e_clock - s_clock
        #obj = -1
        if best_S == None:
            best_S = S
            best_obj = obj
            best_time = etime
        elif best_obj > obj:
            obj = best_obj
            best_S = S
            best_time = etime
        avg += obj
        etime_avg += etime
    avg /= args.nrun
    etime_avg /= args.nrun
    print('min OFV,','average OFV,','time to find the min,','average running time')
    print(best_obj, round(avg,1), round(best_time,2), round(etime_avg, 2)) 
    eprint(len(best_S), args.budget, args.file)

if __name__ == '__main__':
    sys.setrecursionlimit(int(1e9))
    main()
