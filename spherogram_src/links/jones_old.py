"""
Functions needed to calculate the Jones polynomial of K. Still needs work ...
"""
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
import sage.graphs.graph as graph
from sage.rings.rational_field import QQ


def edge_index(edge_datum):
    return edge_datum[2]['edge_index']


def cut(G, T, e):
    """
    Input:
    --A graph G.
    --A spanning tree T for G
    --And edge e of T

    Cutting T along e separates T into two components.
    Returns: The list of edges in G - e connecting the two different components of T-e."""
    if not T.has_edge(*e):
        raise ValueError("e must be an edge of T.")
    S = T.copy()
    S.delete_edge(e)
    (C1, C2) = S.connected_components()
    answer = list()
    for f in G.edges(sort=True, key=edge_index):
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

    Returns: True if e is in T and e is internally active for T, False otherwise. Uses the ordering on G.edges()."""
    if not T.has_edge(*e):
        return False
    for f in cut(G, T, e):
        if edge_index(f) < edge_index(e):
            return False
    return True


def cyc(G, T, e):
    """
    Input:
    --A graph G.
    --A spanning tree T for G
    --And edge e of G not in T

    Adjoining e to T creates a cycle.
    Returns: this cycle."""
    if not G.has_edge(*e):
        raise ValueError("e must be an edge of G.")
    if T.has_edge(*e):
        raise ValueError("e must not be an edge of T.")
    # First thing: catch exceptional case that e is a multiple for an edge in T (giving a 2-cycle).
    try:
        l = T.edge_label(e[0], e[1])
        if isinstance(l, list):
            l = l[0]  # For multigraphs, edge_label returns a list. In this case, it's a list with one element...
        if (e[0], e[1], l) in T.edges(sort=True, key=edge_index):
            return [(e[0], e[1], l), e]
        return [(e[1], e[0], l), e]
    except Exception:
        pass

    # Now the typical case.  First, need to turn T into a Graph which
    # doesn't allow multiedges and also make a copy since we will
    # modify it.
    S = graph.Graph(T.edges(sort=True, key=edge_index))
    S.add_edge(e)
    cb = S.cycle_basis()[0]
    answer = list()
    for i in range(len(cb)):
        l = S.edge_label(cb[i], cb[(i + 1) % len(cb)])
        if S.has_edge(cb[i], cb[(i + 1) % len(cb)], l):
            answer.append((cb[i], cb[(i + 1) % len(cb)], l))
        else:
            answer.append((cb[(i + 1) % len(cb)], cb[i], l))
    return answer


def is_externally_active(G, T, e):
    """
    Input:
    --A graph G.
    --A spanning tree T for G
    --And edge e of G

    Returns: True is e is not in T and e is externally active for T, False otherwise. Uses the ordering on G.edges()."""
    if T.has_edge(*e):
        return False
    for f in cyc(G, T, e):
        if edge_index(f) < edge_index(e):
            return False
    return True


def _edge_sign(K, edge):
    "Returns the sign (+/- 1) associated to given edge in the black graph."
    crossing = edge[2]
    if set(((crossing, 0), (crossing, 1))).issubset(set(edge[0])) or set(((crossing, 0), (crossing, 1))).issubset(set(edge[1])):
        return +1
    return -1


def _Jones_contrib_edge(K, G, T, e, A):
    "Returns the contribution to the Jones polynomial of the specified tree T and edge e."
    # Need to also take crossing sign into account -- A -> 1/A in negative case.
    s = e[2]['sign']
    if is_internally_active(G, T, e):
        return -A**(-3 * s)
    if T.has_edge(*e) and (not is_internally_active(G, T, e)):
        return A**s
    if is_externally_active(G, T, e):
        return -A**(3 * s)
    if (not T.has_edge(*e)) and (not is_externally_active(G, T, e)):
        return A**(-1 * s)


def _Jones_contrib(K, G, T, A):
    "Returns the contribution to the Jones polynomial of the tree T. G should be self.black_graph()."
    answer = 1
    # 2 loops, edges in T and edges not in T
    for e in G.edges(sort=True, key=edge_index):
        answer = answer * _Jones_contrib_edge(K, G, T, e, A)
    return answer


def Jones_poly(K, variable=None, new_convention=False):
    """
    The old convention should really have powers of q^(1/2) for links
    with an odd number of components, but it just multiplies the
    answer by q^(1/2) to get rid of them.  Moroever, the choice of
    value for the unlink is a little screwy, essentially::

      (-q^(1/2) - q^(-1/2))^(n - 1).

    In the new convention, powers of q^(1/2) never appear, i.e. the
    new q is the old q^(1/2) and moreover the value for an n-component
    unlink is (q + 1/q)^(n - 1).  This should match Bar-Natan's paper
    on Khovanov homology.
    """
    if not variable:
        L = LaurentPolynomialRing(QQ, 'q')
        variable = L.gen()
    answer = 0
    L_A = LaurentPolynomialRing(QQ, 'A')
    A = L_A.gen()
    G = K.white_graph()
    for i, labels in enumerate(G.edge_labels()):
        labels['edge_index'] = i
    writhe = K.writhe()
    for T in spanning_trees(G):
        answer = answer + _Jones_contrib(K, G, T, A)
    answer = answer * (-A)**(3 * writhe)
    ans = 0
    for i in range(len(answer.coefficients())):
        coeff = answer.coefficients()[i]
        exp = answer.exponents()[i]
        if new_convention:
            # Now do the substitution A = i q^(1/2) so A^2 = -q
            assert exp % 2 == 0
            ans = ans + coeff * ((-variable)**(exp // 2))
        else:
            ans = ans + coeff * (variable**(exp // 4))
    return ans


def spanning_trees(G):
    """
    NOTE: This code was essentially merged into SageMath proper
    in 2014.  However, because of sorting-related changes needed to
    support Python 3, the *other* code in this file will not work with
    SageMath's version of Graph.spanning_trees.  Hence we retain our
    original version, somewhat modified to work around the Python 3
    issues; specifically, each edge must have an "edge_index" label
    which uniquely identifies it and is sortable.

    Returns a list of all spanning trees.

    If the graph is disconnected, returns the empty list.

    Uses the Read-Tarjan backtracking algorithm [RT75]_.

    EXAMPLES::

        sage: G = Graph([(1,2),(1,2),(1,3),(1,3),(2,3),(1,4)],multiedges=True)
        sage: len(spanning_trees(G))
        8
        sage: G.spanning_trees_count()
        8
        sage: G = Graph([(1,2),(2,3),(3,1),(3,4),(4,5),(4,5),(4,6)],multiedges=True)
        sage: len(spanning_trees(G))
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
        sage: spanning_trees(G)
        [Graph on 6 vertices,
         Graph on 6 vertices,
         Graph on 6 vertices,
         Graph on 6 vertices,
         Graph on 6 vertices,
         Graph on 6 vertices]

    REFERENCES:

    .. [RT75] Read, R. C. and Tarjan, R. E.
      Bounds on Backtrack Algorithms for Listing Cycles, Paths, and Spanning Trees
      Networks, Volume 5 (1975), number 3, pages 237-252.
    """

    def _recursive_spanning_trees(G, forest):
        """
        Returns all the spanning trees of G containing forest
        """
        if not G.is_connected():
            return []

        if G.size() == forest.size():
            return [forest.copy()]
        else:
            # Pick an edge e from G-forest
            for e in G.edges(sort=True, key=edge_index):
                if not forest.has_edge(e):
                    break

            # 1) Recursive call with e removed from G
            G.delete_edge(e)
            trees = _recursive_spanning_trees(G, forest)
            G.add_edge(e)

            # 2) Recursive call with e include in forest
            #
            # e=xy links the CC (connected component) of forest containing x
            # with the CC containing y. Any other edge which does that
            # cannot be added to forest anymore, and B is the list of them
            c1 = forest.connected_component_containing_vertex(e[0])
            c2 = forest.connected_component_containing_vertex(e[1])
            G.delete_edge(e)
            B = G.edge_boundary(c1, c2, sort=False)
            G.add_edge(e)

            # Actual call
            forest.add_edge(e)
            G.delete_edges(B)
            trees.extend(_recursive_spanning_trees(G, forest))
            G.add_edges(B)
            forest.delete_edge(e)

            return trees

    if G.is_connected() and len(G):
        forest = graph.Graph()
        forest.add_vertices(G.vertices(sort=True))
        forest.add_edges(G.bridges())
        return _recursive_spanning_trees(graph.Graph(G, immutable=False, loops=False),
                                         forest)
    else:
        return []
