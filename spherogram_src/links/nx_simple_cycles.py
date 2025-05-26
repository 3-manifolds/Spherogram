"""
========================
Cycle finding algorithms
========================

This code copied out of NetworkX 3.4.2 as a backport for older versions
of NetworkX.  It works at least down to NetworkX 2.6.

The needed feature appeared in NetworkX 3.1 from March 2023.
"""

from collections import defaultdict
from itertools import combinations, product
import networkx as nx

__all__ = [
    "simple_cycles",
]


def simple_cycles(G, length_bound=None):
    """Find simple cycles (elementary circuits) of a graph.

    A "simple cycle", or "elementary circuit", is a closed path where
    no node appears twice.  In a directed graph, two simple cycles are distinct
    if they are not cyclic permutations of each other.  In an undirected graph,
    two simple cycles are distinct if they are not cyclic permutations of each
    other nor of the other's reversal.

    Optionally, the cycles are bounded in length.  In the unbounded case, we use
    a nonrecursive, iterator/generator version of Johnson's algorithm [1]_.  In
    the bounded case, we use a version of the algorithm of Gupta and
    Suzumura [2]_. There may be better algorithms for some cases [3]_ [4]_ [5]_.

    The algorithms of Johnson, and Gupta and Suzumura, are enhanced by some
    well-known preprocessing techniques.  When `G` is directed, we restrict our
    attention to strongly connected components of `G`, generate all simple cycles
    containing a certain node, remove that node, and further decompose the
    remainder into strongly connected components.  When `G` is undirected, we
    restrict our attention to biconnected components, generate all simple cycles
    containing a particular edge, remove that edge, and further decompose the
    remainder into biconnected components.

    Note that multigraphs are supported by this function -- and in undirected
    multigraphs, a pair of parallel edges is considered a cycle of length 2.
    Likewise, self-loops are considered to be cycles of length 1.  We define
    cycles as sequences of nodes; so the presence of loops and parallel edges
    does not change the number of simple cycles in a graph.

    Parameters
    ----------
    G : NetworkX Graph
       A networkx graph. Undirected, directed, and multigraphs are all supported.

    length_bound : int or None, optional (default=None)
       If `length_bound` is an int, generate all simple cycles of `G` with length at
       most `length_bound`.  Otherwise, generate all simple cycles of `G`.

    Yields
    ------
    list of nodes
       Each cycle is represented by a list of nodes along the cycle.

    Examples
    --------
    >>> G = nx.DiGraph([(0, 0), (0, 1), (0, 2), (1, 2), (2, 0), (2, 1), (2, 2)])
    >>> sorted(nx.simple_cycles(G))
    [[0], [0, 1, 2], [0, 2], [1, 2], [2]]

    To filter the cycles so that they don't include certain nodes or edges,
    copy your graph and eliminate those nodes or edges before calling.
    For example, to exclude self-loops from the above example:

    >>> H = G.copy()
    >>> H.remove_edges_from(nx.selfloop_edges(G))
    >>> sorted(nx.simple_cycles(H))
    [[0, 1, 2], [0, 2], [1, 2]]

    Notes
    -----
    When `length_bound` is None, the time complexity is $O((n+e)(c+1))$ for $n$
    nodes, $e$ edges and $c$ simple circuits.  Otherwise, when ``length_bound > 1``,
    the time complexity is $O((c+n)(k-1)d^k)$ where $d$ is the average degree of
    the nodes of `G` and $k$ = `length_bound`.

    Raises
    ------
    ValueError
        when ``length_bound < 0``.

    References
    ----------
    .. [1] Finding all the elementary circuits of a directed graph.
       D. B. Johnson, SIAM Journal on Computing 4, no. 1, 77-84, 1975.
       https://doi.org/10.1137/0204007
    .. [2] Finding All Bounded-Length Simple Cycles in a Directed Graph
       A. Gupta and T. Suzumura https://arxiv.org/abs/2105.10094
    .. [3] Enumerating the cycles of a digraph: a new preprocessing strategy.
       G. Loizou and P. Thanish, Information Sciences, v. 27, 163-182, 1982.
    .. [4] A search strategy for the elementary cycles of a directed graph.
       J.L. Szwarcfiter and P.E. Lauer, BIT NUMERICAL MATHEMATICS,
       v. 16, no. 2, 192-204, 1976.
    .. [5] Optimal Listing of Cycles and st-Paths in Undirected Graphs
        R. Ferreira and R. Grossi and A. Marino and N. Pisanti and R. Rizzi and
        G. Sacomoto https://arxiv.org/abs/1205.2766

    See Also
    --------
    cycle_basis
    chordless_cycles
    """

    if length_bound is not None:
        if length_bound == 0:
            return
        elif length_bound < 0:
            raise ValueError("length bound must be non-negative")

    directed = G.is_directed()
    yield from ([v] for v, Gv in G.adj.items() if v in Gv)

    if length_bound is not None and length_bound == 1:
        return

    if G.is_multigraph() and not directed:
        visited = set()
        for u, Gu in G.adj.items():
            multiplicity = ((v, len(Guv)) for v, Guv in Gu.items() if v in visited)
            yield from ([u, v] for v, m in multiplicity if m > 1)
            visited.add(u)

    # explicitly filter out loops; implicitly filter out parallel edges
    if directed:
        G = nx.DiGraph((u, v) for u, Gu in G.adj.items() for v in Gu if v != u)
    else:
        G = nx.Graph((u, v) for u, Gu in G.adj.items() for v in Gu if v != u)

    # this case is not strictly necessary but improves performance
    if length_bound is not None and length_bound == 2:
        if directed:
            visited = set()
            for u, Gu in G.adj.items():
                yield from (
                    [v, u] for v in visited.intersection(Gu) if G.has_edge(v, u)
                )
                visited.add(u)
        return

    if directed:
        yield from _directed_cycle_search(G, length_bound)
    else:
        yield from _undirected_cycle_search(G, length_bound)


def _directed_cycle_search(G, length_bound):
    """A dispatch function for `simple_cycles` for directed graphs.

    We generate all cycles of G through binary partition.

        1. Pick a node v in G which belongs to at least one cycle
            a. Generate all cycles of G which contain the node v.
            b. Recursively generate all cycles of G \\ v.

    This is accomplished through the following:

        1. Compute the strongly connected components SCC of G.
        2. Select and remove a biconnected component C from BCC.  Select a
           non-tree edge (u, v) of a depth-first search of G[C].
        3. For each simple cycle P containing v in G[C], yield P.
        4. Add the biconnected components of G[C \\ v] to BCC.

    If the parameter length_bound is not None, then step 3 will be limited to
    simple cycles of length at most length_bound.

    Parameters
    ----------
    G : NetworkX DiGraph
       A directed graph

    length_bound : int or None
       If length_bound is an int, generate all simple cycles of G with length at most length_bound.
       Otherwise, generate all simple cycles of G.

    Yields
    ------
    list of nodes
       Each cycle is represented by a list of nodes along the cycle.
    """

    scc = nx.strongly_connected_components
    components = [c for c in scc(G) if len(c) >= 2]
    while components:
        c = components.pop()
        Gc = G.subgraph(c)
        v = next(iter(c))
        if length_bound is None:
            yield from _johnson_cycle_search(Gc, [v])
        else:
            yield from _bounded_cycle_search(Gc, [v], length_bound)
        # delete v after searching G, to make sure we can find v
        G.remove_node(v)
        components.extend(c for c in scc(Gc) if len(c) >= 2)


def _undirected_cycle_search(G, length_bound):
    """A dispatch function for `simple_cycles` for undirected graphs.

    We generate all cycles of G through binary partition.

        1. Pick an edge (u, v) in G which belongs to at least one cycle
            a. Generate all cycles of G which contain the edge (u, v)
            b. Recursively generate all cycles of G \\ (u, v)

    This is accomplished through the following:

        1. Compute the biconnected components BCC of G.
        2. Select and remove a biconnected component C from BCC.  Select a
           non-tree edge (u, v) of a depth-first search of G[C].
        3. For each (v -> u) path P remaining in G[C] \\ (u, v), yield P.
        4. Add the biconnected components of G[C] \\ (u, v) to BCC.

    If the parameter length_bound is not None, then step 3 will be limited to simple paths
    of length at most length_bound.

    Parameters
    ----------
    G : NetworkX Graph
       An undirected graph

    length_bound : int or None
       If length_bound is an int, generate all simple cycles of G with length at most length_bound.
       Otherwise, generate all simple cycles of G.

    Yields
    ------
    list of nodes
       Each cycle is represented by a list of nodes along the cycle.
    """

    bcc = nx.biconnected_components
    components = [c for c in bcc(G) if len(c) >= 3]
    while components:
        c = components.pop()
        Gc = G.subgraph(c)
        uv = list(next(iter(Gc.edges)))
        G.remove_edge(*uv)
        # delete (u, v) before searching G, to avoid fake 3-cycles [u, v, u]
        if length_bound is None:
            yield from _johnson_cycle_search(Gc, uv)
        else:
            yield from _bounded_cycle_search(Gc, uv, length_bound)
        components.extend(c for c in bcc(Gc) if len(c) >= 3)


class _NeighborhoodCache(dict):
    """Very lightweight graph wrapper which caches neighborhoods as list.

    This dict subclass uses the __missing__ functionality to query graphs for
    their neighborhoods, and store the result as a list.  This is used to avoid
    the performance penalty incurred by subgraph views.
    """

    def __init__(self, G):
        self.G = G

    def __missing__(self, v):
        Gv = self[v] = list(self.G[v])
        return Gv


def _johnson_cycle_search(G, path):
    """The main loop of the cycle-enumeration algorithm of Johnson.

    Parameters
    ----------
    G : NetworkX Graph or DiGraph
       A graph

    path : list
       A cycle prefix.  All cycles generated will begin with this prefix.

    Yields
    ------
    list of nodes
       Each cycle is represented by a list of nodes along the cycle.

    References
    ----------
    .. [1] Finding all the elementary circuits of a directed graph.
       D. B. Johnson, SIAM Journal on Computing 4, no. 1, 77-84, 1975.
       https://doi.org/10.1137/0204007
    """
    G = _NeighborhoodCache(G)
    blocked = set(path)
    B = defaultdict(set)  # graph portions that yield no elementary circuit
    start = path[0]
    stack = [iter(G[path[-1]])]
    closed = [False]
    while stack:
        nbrs = stack[-1]
        for w in nbrs:
            if w == start:
                yield path[:]
                closed[-1] = True
            elif w not in blocked:
                path.append(w)
                closed.append(False)
                stack.append(iter(G[w]))
                blocked.add(w)
                break
        else:  # no more nbrs
            stack.pop()
            v = path.pop()
            if closed.pop():
                if closed:
                    closed[-1] = True
                unblock_stack = {v}
                while unblock_stack:
                    u = unblock_stack.pop()
                    if u in blocked:
                        blocked.remove(u)
                        unblock_stack.update(B[u])
                        B[u].clear()
            else:
                for w in G[v]:
                    B[w].add(v)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
