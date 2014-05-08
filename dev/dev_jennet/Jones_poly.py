from __future__ import print_function
"""
Functions needed to calculate the Jones polynomial of K. Still needs work ...
"""
import spherogram
#import spanning_trees
import tree_functions
from sage.symbolic.ring import var

def _crossing_to_edge(K,G,c):
    "Returns the edge of G corresponding to crossing c, where G is either the black graph or white graph."
    edge = list()
    verts = G.vertices()
    black_faces = list() #Strip off the corners (0, 1, 2 or 3) from the crossings
    for v in verts:
        black_faces.append([x[0] for x in v])
    for i in range(len(black_faces)):
        if c in black_faces[i]:
            edge.append(verts[i])
    edge.append(c)
    if len(edge) != 3:
        raise Exception("Did not find two faces incident to c="+repr(c))
    return tuple(edge)

def _dual_spanning_tree(K,T):
    "Returns the spanning tree for the white graph dual to the spanning tree T for the black graph. Here, dual means the edges do not intersect, i.e., no common crossings."
    edges = list()
    G = K.white_graph()
    crossings_to_ignore = [e[2] for e in T.edges()]
    crossings_to_use = [c for c in K.crossings if not c in crossings_to_ignore]
    return Graph([K._crossing_to_edge(G,c) for c in crossings_to_use])

def _edge_sign(K, edge): #Giving weird answers; maybe not working right?
    "Returns the sign (+/- 1) associated to given edge in the black graph."
    crossing = edge[2]
    if set(((crossing,0),(crossing,1))).issubset(set(edge[0])) or set(((crossing,0),(crossing,1))).issubset(set(edge[1])):
        return +1
    return -1

def Kauffman_states(K):
    "Returns the set of Kauffman states for the Alexander polynomial, corresponding to the spanning trees in the black graph. Returns a list of dictionaries, with keys crossings and values faces in the knot projection."
    G = K.black_graph()
    trees = G.spanning_trees()
    marked_edge = ((K.crossings[0],0), K.crossings[0].adjacent[0]) #We (arbitrarily) mark an edge in K.
    states = list()
    for T in trees:
        #Find the root of T. 
        for v in T.vertices():
            if marked_edge[0] in v:
                root = v
        #Orient T out from the root.
        oT = orient_tree(T,root)
        #Find the planar dual spanning tree for the white graph
        dT = _dual_spanning_tree(K,T)
        #Find the root of dT. 
        for v in dT.vertices():
            if marked_edge[0] in v:
                droot = v
        odT = orient_tree(dT,droot)
        state = dict()
        for e in oT.edges():
            state[e[2]] = e[1]
        for e in odT.edges():
            state[e[2]] = e[1]
        states.append(state)
    return states

def _Jones_contrib_edge(K, G, T, e, A):
    "Returns the contribution to the Jones polynomial of the specified tree T and edge e."
    #Need to also take crossing sign into account -- A -> 1/A in negative case.
    s = _edge_sign(K,e)
    if tree_functions.is_internally_active(G,T,e):
        return -A**(-3*s)
    if (e in T.edges()) and (not tree_functions.is_internally_active(G,T,e)):
        return A**s
    if tree_functions.is_externally_active(G,T,e):
        return -A**(3*s)
    if (not e in T.edges()) and (not tree_functions.is_externally_active(G,T,e)):
        return A**(-1*s)
    
def _Jones_contrib(K, G, T, A):
    "Returns the contribution to the Jones polynomial of the tree T. G should be self.black_graph()."
    answer = 1
    # 2 loops, edges in T and edges not in T
    for e in G.edges():
        answer = answer*_Jones_contrib_edge(K,G,T,e,A)
    return answer

def Jones_poly(K,variable=None):
    if not variable:
        variable = var('t')
    answer = 0
    A = var('A')
    G = K.black_graph()
    for T in G.spanning_trees():
        answer = answer + _Jones_contrib(K,G,T,A)
    answer = answer * (-A)**(3*K.writhe()) #Switched this to a -. Sign of exponent seems to be off from standard conventions.
    answer = answer.expand()
    #Python doesn't deal well with rational powers, so renormalizing (divide exponents by 4) is a pain. (Sage would do this fine.)
    ans = 0
    for [coeff, exp] in answer.coeffs():
        ans = ans + coeff*(variable**(-exp/4))
    return ans
