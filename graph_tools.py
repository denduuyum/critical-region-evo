from queue import Queue
import networkit as nk
import networkx as nx
import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.use('TkAgg')
from matplotlib.colors import LogNorm
import copy
import sys
import ArticulationPoint as AP
#print('Hi from graph tool')
#exit()

# It creates only ONE bfs tree (assuming G is a connected graph)
# Returns:
#     - d    : depth of the bfs tree
#     - T    : edges of the bfs tree
#     - par  : parents nodes of the vertices in the bfs tree
def bfs(G, s, lvl_info = None, dist = None):

    Q = Queue()
    Q.put((s, 0))
    T = []
    visited = {}
    par = {}

    visited[s] = True
    d = 0
    while Q.empty() == False:
        v, d = Q.get()
        if dist != None and d > dist:
            break

        if lvl_info != None:
            if len(lvl_info) <= d:
                lvl_info.append([v])
            else:
                lvl_info[d].append(v)

        for u in G.iterNeighbors(v):
            if visited.get(u) == None:
                visited[u] = True
                par[u] = v
                T.append((v, u))
                Q.put((u, d + 1))
    return d, T, par


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


'''
 G - An undirected graph
 u - source vertex
 d - depth info
 vis - array for marking visitation
 P - array for marking parent vertex
 lvl - furtest vertex info
'''
def dfs(G, u, d, vis, P, lvl):
    vis[u] = True
    if lvl[0] < d:
        lvl[0] = d
        lvl[1] = u
    for v in G.neighbors(u):
        if vis[v] == False:
            P[v] = u 
            dfs(G, v, d+1, vis, P, lvl)

# - P    : Array storing parents of the DFS/BFS tree
# - s    : source
# - t    : sink
def get_path_nodes(P, s, t):
    path = [t]
    while t != s:
        t = P[t]
        path.append(t)
    return path

def dfs(G, u, d, vis):
    vis[u] = d
    for v in G[u]:
        if vis[v] == -1:
            dfs(G, v, d + 1, vis)

def get_edge_number(G, n):
    total_degree = 0
    for i in range(n):
        if G.hasNode(i) == True:
            total_degree = total_degree + len(G.neighbors(i))

    return int(total_degree/4)


# Function to read graph from the source            
def read_graph(fin):
    n = int(fin.readline())
    G = nk.Graph(n)
    
    # print("Number of nodes: ", n)    
    while True:
        try:
            line = fin.readline()
        except:
            break

        line = line.split()
        if len(line) == 0:
            break

        x = int(line[0][:-1])
        arr = [int(y) for y in line[1:]]        
        for y in arr:
            G.addEdge(x, y, addMissing = True)

    G.removeMultiEdges()
    return G, n

def read_graph_mapped(fin):
    n = int(fin.readline())
    G = nk.Graph(n)
    d = {}
    id = 0
    
    # print("Number of nodes: ", n)
    while True:
        try:
            line = fin.readline()
        except:
            break

        line = line.split()
        if len(line) == 0:
            break

        x = int(line[0][:-1])
        if d.get(x) == None:
            d[x] = id
            id += 1
        x = d[x]
        arr = [int(y) for y in line[1:]]
        for y in arr:
            if d.get(y) == None:
                d[y] = id
                id += 1
            yy = d[y]
            G.addEdge(x, yy, addMissing = True)

    return G, id


def nx_read_graph_mapped(fin):
    n = int(fin.readline())
    G = nx.Graph()
    d = {}
    id = 0
    
    # print("Number of nodes: ", n)
    while True:
        try:
            line = fin.readline()
        except:
            break

        line = line.split()
        if len(line) == 0:
            break

        x = int(line[0][:-1])
        if d.get(x) == None:
            d[x] = id
            id += 1
        x = d[x]
        arr = [int(y) for y in line[1:]]
        for y in arr:
            if d.get(y) == None:
                d[y] = id
                id += 1
            yy = d[y]
            G.add_edge(x, yy)

    return G, id

            
def write_graph(G, n, fout):
    fout.write(str(G.numberOfNodes()) + '\n')
    for i in range(n):
        if G.hasNode(i) == True:
            fout.write(str(i) + ': ')
            for u in G.neighbors(i):
                fout.write(str(u) + ' ')
            fout.write('\n')

def largest_component_size(G, cc = None):
    if G.numberOfNodes() == 0:
        return 0
    if cc == None:
        cc = nk.components.ConnectedComponents(G)
    cc.run()
    H = [v for k, v in cc.getComponentSizes().items()]
    if len(H) == 0:
        return 0
    return max(H)

def print_graph(G, del_nodes = None):
    G = nk.Graph(G)
    cc = nk.components.ConnectedComponents(G)

    if del_nodes != None:
        for u in del_nodes:
            G.removeNode(u)
            
    cc.run()
    ccGraphSize = [v for k, v in cc.getComponentSizes().items()]
    ccGraphSize = sorted(ccGraphSize, reverse = True)
    print(ccGraphSize)

def eprint_graph(G, del_nodes = None):
    if del_nodes != None:
        G = nk.Graph(G)
        for u in del_nodes:
            G.removeNode(u)

    cc = nk.components.ConnectedComponents(G)
            
    cc.run()
    ccGraphSize = [v for k, v in cc.getComponentSizes().items()]
    ccGraphSize = sorted(ccGraphSize, reverse = True)
    # eprint(ccGraphSize[:20])
    eprint(ccGraphSize)


def size_of_gcc_after_removing(G, del_nodes = None):
    G = nk.Graph(G)
    cc = nk.components.ConnectedComponents(G)

    if del_nodes != None:
        for u in del_nodes:
            G.removeNode(u)
            
    cc.run()
    ccGraphSize = [v for k, v in cc.getComponentSizes().items()]
    ccGraphSize = sorted(ccGraphSize, reverse = True)
    if len(ccGraphSize) == 0:
        return 0
    return ccGraphSize[0]
    

def remove_nodes(G, nodes):
    for u in nodes:
        if G.hasNode(u):
            G.removeNode(u)

def random_path(G, s, d, DEPTH):

    if d == DEPTH:
        return [s]
    if G.degree(s) == 0:
        return [s]
    u = nk.graphtools.randomNeighbor(G, s)
    return random_path(G, u, d + 1, DEPTH) + [u]

def collective_bc(G, bc_data):
    bc2 = copy.copy(bc_data)
    for u in G.iterNodes():
        for v in G.neighbors(u):
            if u != v:
                bc2[u] += bc_data[v]
    return bc2
    
    
def remove_nodes_ret_edges(G, nodes):
    edges = []
    for u in nodes:
        if G.hasNode(u):
            for v in G.neighbors(u):
                edges.append((u, v))
        G.removeNode(u)
    return reversed(edges)

def remove_top_k_nodes(G, rank, k):
    n = len(rank)
    p = [i for i in range(n)]
    p = sorted(p, key = lambda x: rank[x], reverse = True)
    nodes = []
    for i in range(int(k)):
        nodes.append(p[i])
    remove_nodes(G, nodes)
    return nodes

def save_graph(G, title, bc = None, DIR = 'analysis'):
    cc = nk.components.ConnectedComponents(G)
    cc.run()

    for i, c in enumerate(cc.getComponents()):
        l = len(c)
        fname = DIR + '/' + title + '_' + str(i) + '_' + str(l)
        with open(fname, 'w') as f:
            for u in c:
                if bc != None:
                    adj = str(u) + ' ' + str(bc[u]) + ': '
                else:
                    adj = str(u) + ': '
                    
                for v in G.neighbors(u):
                    adj += str(v) + ' '

                f.write(adj + '\n')

def save_cut(nodes, title, DIR = 'analysis'):
    with open(DIR + '/' + title, 'w') as f:
        print(len(nodes))
        for u in nodes:
            f.write(str(u) + '\n')

def remove_self_edges(G):
    for u in G.iterNodes():
        while G.hasEdge(u, u):
            G.removeEdge(u, u)
            
def clean_graph(G):
    remove_self_edges(G)
    # Remove duplicate edges
    G.removeMultiEdges()

def cumulative_bc(G, bc):
    new_bc = copy.copy(bc)
    for u in G.iterNodes():
        for v in G.neighbors(u):
            new_bc[u] += bc[v]
        new_bc[u] /= G.degree(u) + 1
    return new_bc


def draw_graph(G, S = [], S2 = []):
    # Convert the graph into networkx graph
    Gx = nx.Graph()
    for u in G.iterNodes():
        Gx.add_node(u)

    for e in G.iterEdges():
        Gx.add_edge(e[0], e[1])

    sizes = []
    cols = []
    for u in Gx:
        if u in S:
            cols.append('blue')
            sizes.append(5)
        elif u in S2:
            cols.append('red')
            sizes.append(20)
        else:
            cols.append('green')
            sizes.append(5)
        
            
    nx.draw(Gx, node_color=cols, node_size=sizes)
    # nx.draw(Gx, node_color=cols, node_size=sizes, with_labels = True)
    plt.show()



def draw_highlight_top_bc(G, bc):
    # Convert the graph into networkx graph
    Gx = nx.Graph()
    for u in G.iterNodes():
        Gx.add_node(u)
    for e in G.iterEdges():
        Gx.add_edge(e[0], e[1])

    sizes = []
    cols = []
    mx = max([bc[u] for u in Gx])
    mn = min([bc[u] for u in Gx])
    l = mx - mn
    lg_norm = LogNorm()
    log_data = lg_norm([bc[u] for u in Gx])
    i = 0
    for u in Gx:
        cols.append((1, 1-log_data[i], 1-log_data[i], 1))
        # cols.append((1, 0, 0, 1))
        # cols.append((1, (bc[u] - mn) / l, (bc[u] - mn) / l, 1))
        i += 1
        sizes.append(50)
    print(cols)
    nx.draw(Gx, node_size=sizes, node_color = cols)
    # nx.draw(Gx, node_size=sizes, s = 500, cmap=log_data)
    plt.show()

def avg(x, k = 3):
    y = []
    for i in range(k, len(x)):
        y.append(sum(x[i - k: i])/k)
    return y
    
def draw_2graph(x1, x2, y1, y2, z1, z2, fname = "filename.svg"):
    plt.subplot('121')
    xx1 = [i for i in range(1, len(x1) + 1)]
    yy1 = [i for i in range(1, len(y1) + 1)]
    zz1 = [i for i in range(1, len(z1) + 1)]
    plt.plot(xx1, x1, 'o-', yy1, y1, 'x-', zz1, z1, '+-')
    plt.subplot('122')
    xx2 = [i for i in range(1, len(x2) + 1)]
    yy2 = [i for i in range(1, len(y2) + 1)]
    zz2 = [i for i in range(1, len(z2) + 1)]
    plt.plot(xx2, x2, 'o-', yy2, y2, 'x-', zz2, z2, '+-')
    plt.savefig(fname)
    

class DrawGraph:

    def __init__(self, G):
        self.G = nx.Graph()
        for u in G.iterNodes():
            self.G.add_node(u)

        for e in G.iterEdges():
            self.G.add_edge(e[0], e[1])

        self.pos = nx.spring_layout(self.G)

    def remove_nodes(self, S):
        for u in S:
            if self.G.has_node(u):
                self.G.remove_node(u)

    def draw(self, S):
        # Convert the graph into networkx graph
        sizes = []
        cols = []
        pos = {}
        for u in self.G:
            pos[u] = self.pos[u]
            if u in S:
                cols.append('red')
                sizes.append(2)
            else:
                cols.append('black')
                sizes.append(1)
            
        nx.draw(self.G, pos = pos, node_color=cols, node_size=sizes)
        plt.show()

def cumulative_bc(G, bc):
    new_bc = {}
    for u in G.iterNodes():
        s = 0
        for v in G.iterNeighbors(u):
            s += bc[v]

        new_bc[u] = bc[u] * s
        

    return new_bc


def get_core(G):
    deleted_aps = []
    while True:
        ap = AP.ArticulationPoint()
        aps = ap.AP(nk.components.ConnectedComponents.extractLargestConnectedComponent(G))
        deleted_aps += aps
        if len(aps) == 0:
            break
        for u in aps:
            G.removeNode(u)

    return deleted_aps
    
def get_ap(G):
    ap = AP.ArticulationPoint()
    aps = ap.AP(nk.components.ConnectedComponents.extractLargestConnectedComponent(G))
    for u in aps:
        G.removeNode(u)

    return aps
    
