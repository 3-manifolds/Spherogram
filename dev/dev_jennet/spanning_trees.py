from __future__ import print_function
"""
Graph functions to be used by the Jones polynomial function in links.py
Potentially could be implemented in sage: 
http://www.sagemath.org/doc/reference/graphs/sage/graphs/graph.html
"""

import sage.graphs.graph as graph

def spanning_trees(G):
    """
    Function 'span' from Read paper
    """

    if G.is_connected():
        part_G = graph.Graph([])
        part_G.add_vertices(G.vertices())
        Bridge = BridgeFinder()
        part_G.add_edges(Bridge.fast_find_bridges(G))
        return rec(part_G,G)
    else:
        return []

def rec(part_G,G):
    """
    Recursive function 'rec' from Read paper
    """
    trees = []

    if not G.is_connected():
        return []

    if len(G.edges()) == len(part_G.edges()):
        trees +=[part_G.copy()]
    else:

        # let e be an edge not in part_G
        X = G.edges()
  
        for y in part_G.edges():
            X.remove(y)
        e = X[0]

        # add e to part_G
        part_G.add_edge(e)
        B = get_B(part_G,G)

        # delete all edges in B from G
        G.delete_edges(B)

        # rec
        trees += rec(part_G, G)

        # add all edges in B to G
        G.add_edges(B)

        # delete e from part_G and from G
        G.delete_edge(e)
        part_G.delete_edge(e)

        # let C be the set of bridges which are not tree edges
        Bridge = BridgeFinder()
        C = Bridge.fast_find_bridges(G)
 
        for x in part_G.edges():
            if x in C:
                C.remove(x)
 
        # add all edges in C to part_G
        part_G.add_edges(C)

        # rec
        trees += rec(part_G, G)        
        part_G.delete_edges(C)
        G.add_edge(e)

    return trees

def get_B(part_G,G):
    """ 
    Returns the set of edges not in part_G joining vertices already connected in part_G
    """
    B = []
    comps = part_G.connected_components()
    vc_dict = dict()
    for i in range(len(comps)):
        for v in comps[i]:
            vc_dict[v] = i
    
    X = G.edges()
    for y in part_G.edges(): 
        X.remove(y)
    
    for e in X:
        if vc_dict[e[0]] == vc_dict[e[1]]:
            B.append(e)
  
    return B

class BridgeFinder():
    """
    Used to find all bridges of a graph G in an efficient way.
    Follows the solution to 4.1.36 in Sedgewick's _Algorithms_ (4th ed.)
    """

    def __init__(self):
        self.cnt = 0    
        self.low = dict()
        self.pre = dict()
        self.bridges = []

    def fast_find_bridges(self,G):
        verts = G.vertices()
        self.low = dict({v:-1 for v in verts})
        self.pre = dict({v:-1 for v in verts})

        for v in verts:
            if self.pre[v] == -1:
                self.dfs(G, v, v)

        return self.bridges

    def dfs(self,G, u, v):
        self.cnt += 1
        self.pre[v] = self.cnt
        self.low[v] = self.pre[v];
        for w in G.neighbor_iterator(v):
            if (self.pre[w] == -1):
                self.dfs(G, v, w)
                self.low[v] = min(self.low[v], self.low[w])
                if (self.low[w] == self.pre[w]):
                    #print(repr(v) + "-" + repr(w) + " is a bridge")
                    self.bridges.append( (min(v,w), max(v,w), None) )
            # update low number - ignore reverse of edge leading to v
            elif (w != u):
                self.low[v] = min(self.low[v], self.pre[w])
