from __future__ import print_function
"""
Functions needed to calculate the Jones polynomial of K. Still needs work ...
"""
from sage.symbolic.ring import var
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
import sage.graphs.graph as graph
from sage.rings.rational_field import QQ

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

def _edge_sign(K, edge):
    "Returns the sign (+/- 1) associated to given edge in the black graph."
    crossing = edge[2]
    if set(((crossing,0),(crossing,1))).issubset(set(edge[0])) or set(((crossing,0),(crossing,1))).issubset(set(edge[1])):
        return +1
    return -1

def _old_Jones_contrib_edge(K, G, T, e, A):
    "Returns the contribution to the Jones polynomial of the specified tree T and edge e."
    #Need to also take crossing sign into account -- A -> 1/A in negative case.
    s = _edge_sign(K,e)
    if is_internally_active(G,T,e):
        return -A**(-3*s)
    if (e in T.edges()) and (not is_internally_active(G,T,e)):
        return A**s
    if is_externally_active(G,T,e):
        return -A**(3*s)
    if (not e in T.edges()) and (not is_externally_active(G,T,e)):
        return A**(-1*s)
    
def _old_Jones_contrib(K, G, T, A):
    "Returns the contribution to the Jones polynomial of the tree T. G should be self.black_graph()."
    answer = 1
    # 2 loops, edges in T and edges not in T
    for e in G.edges():
        answer = answer*_old_Jones_contrib_edge(K,G,T,e,A)
    return answer

def old_Jones_poly(K,variable=None):
    if not variable:
        L = LaurentPolynomialRing(QQ,'q')
        variable = L.gen()
    answer = 0
    A = var('A')
    G = K.black_graph()
    for T in spanning_trees(G):
        answer = answer + _old_Jones_contrib(K,G,T,A)
    answer = answer * (-A)**(3*K.writhe())
    answer = answer.expand()
    #Python doesn't deal well with rational powers, so renormalizing (divide exponents by 4) is a pain. (Sage would do this fine.)
    ans = 0
    for [coeff, exp] in answer.coefficients():
        ans = ans + coeff*(variable**(exp/4))
    return ans

def _Jones_contrib_edge(K, G, T, e, A):
    "Returns the contribution to the Jones polynomial of the specified tree T and edge e."
    #Need to also take crossing sign into account -- A -> 1/A in negative case.
    s = e[2]['sign']
    if is_internally_active(G,T,e):
        return -A**(-3*s)
    if (e in T.edges()) and (not is_internally_active(G,T,e)):
        return A**s
    if is_externally_active(G,T,e):
        return -A**(3*s)
    if (not e in T.edges()) and (not is_externally_active(G,T,e)):
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
        L = LaurentPolynomialRing(QQ,'q')
        variable = L.gen()
    answer = 0
    L_A = LaurentPolynomialRing(QQ,'A')
    A = L_A.gen()
    G = K.white_graph()
    writhe = K.writhe()
    for T in spanning_trees(G):
        answer = answer + _Jones_contrib(K,G,T,A)
    answer = answer * (-A)**(3*writhe)
    #Python doesn't deal well with rational powers, so renormalizing (divide exponents by 4) is a pain. (Sage would do this fine.)
    ans = 0
    for i in range(len(answer.coefficients())):
        coeff = answer.coefficients()[i]
        exp = answer.exponents()[i]
        ans = ans + coeff*(variable**(exp/4))
    return ans


    
def spanning_trees(G):
    """
    Returns a list of all spanning trees.

    If the graph is disconnected, returns the empty list.

    Uses the Read-Tarjan backtracking algorithm [RT75]_.

    EXAMPLES::

        sage: G = Graph([(1,2),(1,2),(1,3),(1,3),(2,3),(1,4)],multiedges=True)
        sage: len(G.spanning_trees())
        8
        sage: G.spanning_trees_count()
        8
        sage: G = Graph([(1,2),(2,3),(3,1),(3,4),(4,5),(4,5),(4,6)],multiedges=True)
        sage: len(G.spanning_trees())
        6
        sage: G.spanning_trees_count()
        6

    .. SEEALSO::

        - :meth:`~sage.graphs.generic_graph.GenericGraph.spanning_trees_count`
          -- counts the number of spanning trees.

        - :meth:`~sage.graphs.graph.Graph.random_spanning_tree`
          -- returns a random spanning tree.

    TESTS:

    Works with looped graphs::

        sage: g = Graph({i:[i,(i+1)%6] for i in range(6)})
        sage: g.spanning_trees()
        [Graph on 6 vertices,
         Graph on 6 vertices,
         Graph on 6 vertices,
         Graph on 6 vertices,
         Graph on 6 vertices,
         Graph on 6 vertices]

    REFERENCES:

    .. [RT75] Read, R. C. and Tarjan, R. E.
      Bounds on Backtrack Algoritms for Listing Cycles, Paths, and Spanning Trees
      Networks, Volume 5 (1975), numer 3, pages 237-252.
    """

    def _recursive_spanning_trees(G,forest):
        """
        Returns all the spanning trees of G containing forest
        """
        if not G.is_connected():
            return []

        if G.size() == forest.size():
            return [forest.copy()]
        else:
            # Pick an edge e from G-forest
            for e in G.edge_iterator():
                if not forest.has_edge(e):
                    break

            # 1) Recursive call with e removed from G
            G.delete_edge(e)
            trees = _recursive_spanning_trees(G,forest)
            G.add_edge(e)

            # 2) Recursive call with e include in forest
            #
            # e=xy links the CC (connected component) of forest containing x
            # with the CC containing y. Any other edge which does that
            # cannot be added to forest anymore, and B is the list of them
            c1 = forest.connected_component_containing_vertex(e[0])
            c2 = forest.connected_component_containing_vertex(e[1])
            G.delete_edge(e)
            B = G.edge_boundary(c1,c2,sort=False)
            G.add_edge(e)

            # Actual call
            forest.add_edge(e)
            G.delete_edges(B)
            trees.extend(_recursive_spanning_trees(G,forest))
            G.add_edges(B)
            forest.delete_edge(e)

            return trees

    if G.is_connected() and len(G):
        forest = graph.Graph([])
        forest.add_vertices(G.vertices())
        forest.add_edges(G.bridges())
        return _recursive_spanning_trees(graph.Graph(G,immutable=False,loops=False), forest)
    else:
        return []
