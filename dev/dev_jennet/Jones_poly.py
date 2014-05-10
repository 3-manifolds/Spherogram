from __future__ import print_function
"""
Functions needed to calculate the Jones polynomial of K. Still needs work ...
"""
import spherogram
#import spanning_trees
import tree_functions
from sage.symbolic.ring import var
import sage.graphs.graph as graph


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
    #Label the edges...
    newG = graph.Graph()
    newG.allow_multiple_edges(True)
    for i in range(len(G.edges())):
        e = G.edges()[i]
        newG.add_edge(e[0],e[1],label=i)
    for T in newG.spanning_trees():
        answer = answer + _Jones_contrib(K,newG,T,A)
    answer = answer * (-A)**(-3*K.writhe())
    answer = answer.expand()
    #Python doesn't deal well with rational powers, so renormalizing (divide exponents by 4) is a pain. (Sage would do this fine.)
    ans = 0
    for [coeff, exp] in answer.coeffs():
        ans = ans + coeff*(variable**(-exp/4))
    return ans
