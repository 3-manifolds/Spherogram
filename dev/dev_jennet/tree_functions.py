import sage.graphs.graph as graph
from sage.graphs.graph import *

def subgraphs(G):
    "Returns a list of all graphs gotten by taking a subset of the edges of G."
    if not G.edges():
        return [G,]
    e = G.edges()[0]
    H = G.copy()
    H.delete_edge(e)
    without_e = subgraphs(H)
    with_e = [x.copy() for x in without_e]
    for x in with_e:
        x.add_edge(e)
    return without_e + with_e

def orient_tree(T, root):
    "Returns a digraph by orienting T outwards from root."
    if not T.is_tree():
        raise Exception("Input must be a tree.")
    if not T.edges():
        verts = T.vertices()
        v = verts[0]
        return graph.DiGraph({v:list()})
    H = T.copy()
    H.delete_vertex(root)
    root_edges = list()
    for e in T.edges_incident(root):
        if e[0] == root:
            root_edges.append(e)
        else:
            root_edges.append((e[1],e[0],e[2]))
    pieces = list() #Orient the subtrees gotten by deleting root.
    for e in root_edges:
        comp = H.subgraph(H.connected_component_containing_vertex(e[1]))
        pieces.append(orient_tree(comp,e[1]))
    all_edges = root_edges
    for P in pieces:
        all_edges = all_edges+P.edges()
    return graph.DiGraph(all_edges)

def cut(G,T,e):
    """
    Input:
    --A graph G.
    --A spanning tree T for G
    --And edge e of T
    
    Cutting T along e separates T into two components.
    Returns: The list of edges in G - e connecting the two different components of T-e."""
    if not e in T.edges():
        raise Exception("e must be an edge of T.")
    H = G.copy()
    S = T.copy()
    S.delete_edge(e)
    (C1, C2) = S.connected_components()
    answer = list()
    for f in G.edges():
        if (f[0] in C1 and f[1] in C2) or (f[0] in C2 and f[1] in C1):
            if f != e:
                answer.append(f)
    return answer
    

def is_internally_active(G, T, e):
    """
    Input:
    --A graph G.
    --A spanning tree T for G
    --And edge e of G

    Returns: True is e is in T and e is internally active for T, False otherwise. Uses the ordering on G.edges()."""
    if not e in T.edges():
        return False
    edges = G.edges()
    for f in cut(G,T,e):
        if edges.index(f)<edges.index(e):
            return False
    return True

def cyc(G,T,e):
    """
    Input:
    --A graph G.
    --A spanning tree T for G
    --And edge e of G not in T

    Adjoining e to T creates a cycle.
    Returns: this cycle."""
    if not e in G.edges():
        raise Exception("e must be an edge of G.")
    if e in T.edges():
        raise Exception("e must not be an edge of T.")
    #First thing: catch exceptional case that e is a multiple for an edge in T (giving a 2-cycle).
    try:
        l = T.edge_label(e[0],e[1])
        if isinstance(l,list):
            l = l[0] #For multigraphs, edge_label returns a list. In this case, it's a list with one element...
        if (e[0],e[1], l) in T.edges():
            return [(e[0],e[1],l),e]
        return [(e[1],e[0],l),e]
    except:
        pass
    #Now the typical case.
    S = graph.Graph(T.edges()) #Hack because otherwise sage thinks this is a multigraph, and cycle_basis is not implemented for such.
    S.add_edge(e)
    cb = S.cycle_basis()[0]
    answer = list()
    for i in range(len(cb)):
        l = S.edge_label(cb[i],cb[(i+1)%len(cb)])
        if (cb[i],cb[(i+1)%len(cb)],l) in S.edges():
            answer.append((cb[i],cb[(i+1)%len(cb)],l))
        else:
            answer.append((cb[(i+1)%len(cb)],cb[i],l))
    return answer

def is_externally_active(G,T,e):
    """
    Input:
    --A graph G.
    --A spanning tree T for G
    --And edge e of G

    Returns: True is e is not in T and e is externally active for T, False otherwise. Uses the ordering on G.edges()."""    
    if e in T.edges():
        return False
    edges = G.edges()
    for f in cyc(G,T,e):
        if edges.index(f)<edges.index(e):
            return False
    return True
