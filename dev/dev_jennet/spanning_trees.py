from __future__ import print_function
"""
Graph functions to be used by the Jones polynomial function in links.py
Potentially could be implemented in sage: 
http://www.sagemath.org/doc/reference/graphs/sage/graphs/graph.html
"""

import sage.graphs.graph as graph
import bridge_finding

def spanning_trees(G):
    """
    Function 'span' from Read paper
    """

    if G.is_connected():
        part_G = graph.Graph([])
        part_G.add_vertices(G.vertices())
        part_G.add_edges(bridge_finding.find_bridges(G))
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
        C = bridge_finding.find_bridges(G)
 
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

