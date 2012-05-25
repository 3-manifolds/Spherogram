from __future__ import print_function
"""
Python implementation of graphs, reduced graphs, directed graphs, fat
graphs (graphs with ordered adjacency lists) and fat directed graphs.

Vertices are arbitrary hashable python objects. Graph methods must
not change vertices.

To instantiate a graph, provide a list of edges, as pairs of vertices,
and an optional list of additional vertices, which may be isolated.
Vertices are saved as a set, so redundancies in the vertex list will
be ignored.  Edges are saved as a list, so multiple edges are allowed.

If G is a graph and v a vertex of G then G[v] is a list of distinct
vertices adjacent to v.  If G is a digraph, then G[v] is the list of
distinct heads of directed edges with tail v.

When called as a function, an edge is the non-trivial involution of
its endpoints; calling it with one endpoint returns the other one.

One can iterate over edges as:

for x, y in G.edges:
  ...

where x and y will be the two endpoints of the edge (tail and head
in the case of directed edges).
"""

try:
    import sage.graphs.graph
    _within_sage = True
except ImportError:
    from .planarity import planar
    
from collections import deque
import operator

class Edge:
    """
    A basic edge.
    """
    def __init__(self, object1, object2):
        self.ends = [object1, object2]

    def __contains__(self, object):
        """
        The endpoints of an edge are "in" the edge.
        """
        return object in self.ends
    
    def __call__(self, end):
        """
        Calling an edge with one endpoint returns the other one.
        """
        try:
            return [v for v in self if v != end].pop()
        except IndexError:
            raise ValueError('Vertex is not an endpoint')

    def __repr__(self):
        return '%s --- %s'%tuple(self.ends)

    def __iter__(self):
        return self.ends.__iter__()

    def is_loop(self):
        return self.ends[0] == self.ends[1]
        
class DirectedEdge(Edge):
    """
    An edge with a head and a tail.
    """
    def head(self):
        return self.ends[1]

    def tail(self):
        return self.ends[0]

    def __repr__(self):
        return '%s --> %s'%tuple(self.ends)

class EdgesBFO:
    """
    Iterator for non-loop edges of a graph in the complement of
    a forbidden set, ordered by distance from the source.
    
    Yields triples (e, v, f) where e and f are edges containing v
    and e precedes f in the depth-first ordering.
    """
    def __init__(self, graph, source, forbidden=set()):
        self.graph = graph
        self.initial = initial = graph.incident(source) - forbidden
        self.seen = forbidden | initial
        self.fifo = deque([ (None, source, e) for e in initial ])

    def __iter__(self):
        return self
    
    def __next__(self):
        if self.fifo:
            parent, vertex, child = self.fifo.popleft()
            new_edges = self.graph.incident(child(vertex)) - self.seen
            self.seen |= new_edges
            self.fifo.extend([(child, child(vertex), edge)
                              for edge in new_edges])
            return parent, vertex, child
        else:
            raise StopIteration

    def next(self):  #For Python 2 compatibility
        return self.__next__()

class Graph:
    """
    A set of vertices and a set of edges joining pairs of vertices.
    Vertices are arbitrary hashable objects.
    """
    edge_class = Edge

    def __init__(self, pairs=[], singles=[]):
        self.vertices = set()
        self.edges = set()
        self.Edge = self.__class__.edge_class
        for x, y in pairs:
            self.add_edge(x,y)
        for vertex in singles:
            self.add_vertex(vertex)

    def __repr__(self):
        V = 'Vertices:\n  ' + '\n  '.join([str(v) for v in self.vertices])
        E = 'Edges:\n  ' + '\n  '.join([str(e) for e in self.edges])
        return '%s\n%s'%(V,E)

    def __getitem__(self, vertex):
        """
        Return a list of adjacent vertices.
        """
        return [e(vertex) for e in self.incident(vertex)]

    def add_edge(self, x, y):
        edge = self.Edge(x,y)
        self.edges.add(edge)
        self.vertices.update(edge)

    def add_vertex(self, hashable):
        self.vertices.add(hashable)
        
    def incident(self, vertex):
        """
        Return the set of non-loops incident to the vertex.
        """
        return set([ e for e in self.edges
                     if vertex in e and not e.is_loop() ])

    def valence(self, vertex):
        valence = 0
        for e in self.edges:
            if vertex in e:
                valence += 1
            if e(vertex) == vertex:
                valence += 1
        return valence

    def components(self, deleted_vertices=[]):
        """
        Return the vertex sets of the connected components of the
        graph obtained by removing the deleted_vertices and any edges
        incident to them.

        >>> G = Graph([(0,1),(1,2),(2,0),(2,3),(3,4),(4,2)])
        >>> G.components()
        [set([0, 1, 2, 3, 4])]
        >>> G.components(deleted_vertices=[2])
        [set([3, 4]), set([0, 1])]
        >>> G.components(deleted_vertices=[0])
        [set([1, 2, 3, 4])]
        """
        forbidden = set()
        for vertex in deleted_vertices:
            forbidden |= self.incident(vertex)
        vertices, result = list(self.vertices - set(deleted_vertices)), []
        while vertices:
            component, start = set(), vertices.pop()
            component.add(start)
            for parent, vertex, child in EdgesBFO(self, start, forbidden):
                new_vertex = child(vertex)
                component.add(new_vertex)
                if new_vertex in vertices:
                    vertices.remove(new_vertex)
            result.append(component)
        return result

    def is_connected(self, deleted_vertices=[]):
        """
        Determine whether the graph obtained by removing the
        deleted_vertices and incident edges is connected.

        >>> G = Graph([(0,1),(1,2),(2,0),(2,3),(3,4),(4,2)])
        >>> G.is_connected(deleted_vertices=[2])
        False
        """
        return len(self.components(deleted_vertices)) <= 1

    def one_min_cut(self, source, sink, capacity=None):
        """
        Find one minimal cut which separates source from sink.

        Returns the cut set of vertices on the source side, the
        set of edges that cross the cut, a maximal family of
        weighted edge-disjoint paths from source to sink, and
        the set of edges with non-zero residual.
        """
        if sink == source:
            return None
        if capacity is None:
            residual = dict.fromkeys(self.edges, 1)
        else:
            residual = dict.copy(capacity)
        full_edges = set([e for e in self.edges if residual[e] == 0])
        children = {}
        for vertex in self.vertices:
            children[vertex] = set()
        cut_set = set()
        path_list = []
        while True:
            # Try to find a new path from source to sink
            parents, cut_set, reached_sink = {}, set([source]), False
            for parent, vertex, child in EdgesBFO(self, source, full_edges):
                 parents[child] = (parent, vertex)
                 cut_set.add(child(vertex))
                 if child(vertex) == sink:
                     reached_sink = True
                     break
            # If we did not get to the sink, we visited every vertex
            # reachable from the source, thereby providing the cut
            # set.
            if not reached_sink:
                break
            # If we got to the sink, do the bookkeeping and continue.
            path = deque()
            weight = residual[child]
            while True:
                path.appendleft( (vertex, child) )
                children[vertex].add(child)
                if vertex == source:
                    break
                child = parent
                parent, vertex = parents[child]
                weight = min(weight, residual[child])
            for vertex, edge in path:
                residual[edge] -= weight
                if residual[edge] == 0:
                    full_edges.add(edge)
            path_list.append( (weight, path) )
        # Find the cut edges.
        cut_edges = set()
        for vertex in cut_set:
            cut_edges |= set([edge for edge in self.edges
                              if vertex in edge
                              and edge(vertex) not in cut_set])
        unsaturated = [ e for e in self.edges if residual[e] > 0 ]
        return {'set': cut_set, 'edges': cut_edges, 'paths': path_list,
                'unsaturated': unsaturated}

    def is_planar(self):
        """
        Return the planarity.
        """
        if _within_sage:
            S = self.sage()
            return S.is_planar()
        return planar(self)
    
    def merge(self, V1, V2):
        """
        Merge two vertices and remove all edges between them.  The
        vertices must support the __or__ method (e.g. all vertices
        might be frozensets).  The two vertices V1 and V2 are replaced
        by V1|V2.

        >>> G = Graph([(0,1),(1,2),(2,0)]).mergeable()
        >>> F = lambda x: frozenset([x])
        >>> G.merge(F(1),F(2))
        >>> G
        Vertices:
          frozenset([1, 2])
          frozenset([0])
        Edges:
          frozenset([1, 2]) --- frozenset([0])
          frozenset([0]) --- frozenset([1, 2])
        """
        new_vertex = V1|V2
        if new_vertex in self.vertices:
            raise ValueError('Merged vertex already exists!')
        self.edges -= set([e for e in self.edges if V1 in e and V2 in e])
        self.vertices.remove(V1)
        self.vertices.remove(V2)
        old_vertices = (V1, V2)
        self.vertices.add(new_vertex)
        for edge in [e for e in self.edges if V1 in e or V2 in e]:
            self.edges.remove(edge)
            x, y = edge
            if x in old_vertices:
                self.edges.add( self.Edge(new_vertex, y) )
            if y in old_vertices:
                self.edges.add( self.Edge(x, new_vertex) )

    def mergeable(self):
        """
        Return a new graph whose vertices are singleton frozensets,
        each containing a vertex of this graph, and with edges
        corresponding to those of this graph.  The new graph
        will support merging vertices.
        """
        V = [frozenset([v]) for v in self.vertices]
        E = [self.Edge(frozenset([x]), frozenset([y])) for x, y in self.edges]
        return self.__class__(E,V)

class ReducedGraph(Graph):
    """
    A graph with at most one edge between any two vertices,
    but having edges with multiplicities.
    """
    def __init__(self, pairs=[], singles=[]):
        self.multiplicities = {}
        Graph.__init__(self, pairs, singles) 

    def find_edge(self, vertex1, vertex2):
        vertices = set([vertex1, vertex2])
        for edge in self.edges:
            if set(edge.ends) == vertices:
                return edge

    def add_edge(self, x, y, multiplicity=1):
        edge = self.find_edge(x, y)
        if edge:
            self.multiplicities[edge] += multiplicity
        else:
            new_edge = self.Edge(x,y)
            self.vertices.update([x,y])
            self.edges.add(new_edge)
            self.multiplicities[new_edge] = multiplicity

    def multi_valence(self, vertex):
        """
        Return the valence of a vertex, counting edge multiplicities.
        """
        valence = 0
        for e in self.edges:
            if vertex in e:
                valence += self.multiplicities[e]
            if e(vertex) == vertex:
                valence += self.multiplicities[e]
        return valence

    def one_min_cut(self, source, sink):
        cut = Graph.one_min_cut(self, source, sink,
                            self.multiplicities.copy())
        cut['size'] = sum([ self.multiplicities[e] for e in cut['edges'] ])
        return cut

    def cut_pairs(self):
        """
        Return a list of cut_pairs.  The graph is assumed to be
        connected and to have no cut vertices.
        """
        pairs = []
        majors = [v for v in self.vertices if self.valence(v) > 2]
        if len(majors) == 2:
            v, V = majors
            if self.valence(v) == 3:
                return []
            edge = self.find_edge(v,V)
            if not edge or self.multiplicities[edge] < 2:
                return []
            return majors
        major_set = set(majors)
        for n in xrange(1,len(majors)):
            for v in majors[:n]:
                pair = (v, majors[n])
                components = self.components(deleted_vertices=pair)
                if len(components) > 2:
                    pairs.append(pair)
                elif len(components) == 2:
                    M0 = len(major_set & components[0])
                    M1 = len(major_set & components[1])
                    edge = self.find_edge(*pair)
                    if edge:
                        if M0 or M1:
                            pairs.append(pair)
                            continue
                    else:
                        if M0 and M1:
                            pairs.append(pair)
        return pairs
    
class Digraph(Graph):
    edge_class = DirectedEdge

    def incident(self, vertex):
        """
        Return the set of non-loops which begin at the vertex.
        """
        return set([ e for e in self.edges
                     if e.tail() == vertex and e.head() != vertex ])

    def components(self):
        """
        Return the vertex sets of the strongly connected components.

        >>> G = Digraph([(0,1),(0,2),(1,2),(2,3),(3,1)])
        >>> G.components()
        [frozenset([1, 2, 3]), frozenset([0])]
        >>> G = Digraph([(0,1),(0,2),(1,2),(2,3),(1,3)])
        >>> G.components()
        [frozenset([3]), frozenset([2]), frozenset([1]), frozenset([0])]
        """
        return StrongConnector(self).components

    def is_connected(self):
        return len(self.components()) <= 1

    def component_DAG(self):
        """
        Return the acyclic digraph whose vertices are the strong
        components of this digraph.  Two components are joined by an
        edge if this digraph has an edge from one component to the
        other.

        >>> G = Digraph([(0,1),(0,2),(1,2),(2,3),(3,1)])
        >>> G.component_DAG()
        Vertices:
          frozenset([1, 2, 3])
          frozenset([0])
        Edges:
          frozenset([0]) --> frozenset([1, 2, 3])
        """
        return StrongConnector(self).DAG()
        
class StrongConnector:
    """
    Finds strong components of a digraph using Tarjan's algorithm;
    see http://en.wikipedia.org/wiki/
    Tarjan%27s_strongly_connected_components_algorithm
    """
    def __init__(self, digraph):
        self.digraph = digraph
        self.seen = []
        self.unclassified = []
        self.components = []
        self.root = {}
        self.which_component = {}
        self.links = set()
        for vertex in self.digraph.vertices:
            if vertex not in self.seen:
                self.search(vertex)
                
    def search(self, vertex):
        self.root[vertex] = len(self.seen)
        self.seen.append(vertex)
        self.unclassified.append(vertex)
        for child in self.digraph[vertex]:
            if child not in self.seen:
                self.search(child)
                self.root[vertex] = min(self.root[child],
                                        self.root[vertex])
            elif child in self.unclassified:
                self.root[vertex] = min(self.seen.index(child),
                                        self.root[vertex])
            else:
                self.links.add((vertex, child))
        if self.root[vertex] == self.seen.index(vertex):
            component = []
            while True:
                child = self.unclassified.pop()
                component.append(child)
                if child == vertex:
                    break
            if self.unclassified:
                self.links.add((self.unclassified[-1], vertex))
            component = frozenset(component)
            self.components.append(component)
            for child in component:
                self.which_component[child] = component

    def DAG(self):
        """
        Return the acyclic directed graph whose vertices are the
        strong components of the underlying digraph.  There is an edge
        joining two components if and only if there is an edge of the
        underlying digraph having an endpoint in each component.
        Using self.links, rather than self.digraph.edges, is a slight
        optimization.
        """
        edges = set()
        for tail, head in self.links:
            dag_tail= self.which_component[tail]
            dag_head = self.which_component[head]
            if dag_head != dag_tail:
                edges.add( (dag_tail, dag_head) )
        return Digraph(edges, self.components)

class Poset(set):
    """
    A partially ordered set, generated from an acyclic directed graph.
    Instantiate with a Digraph.  A ValueError exception is raised if the
    Digraph contains a cycle.
    """
    def __init__(self, digraph):
        self.elements = set(digraph.vertices)
        self.larger = {}
        self.smaller = {}
        self.successors = {}
        self.closed = set()
        for vertex in self:
            self.larger[vertex] = set()
            self.smaller[vertex] = set()
            self.successors[vertex] = set(digraph[vertex])
        seen = []
        for vertex in self:
            if vertex not in seen:
                self.search(vertex, seen, digraph)

    def __iter__(self):
        return self.elements.__iter__()

    def __len__(self):
        return len(self.elements)
    
    def search(self, vertex, seen, digraph):
        seen.append(vertex)
        for child in digraph[vertex]:
            if child in self.smaller[vertex]:
                raise ValueError('Digraph is not acyclic.')
            self.smaller[child].add(vertex)
            self.smaller[child] |= self.smaller[vertex]
            self.search(child, seen, digraph)
            self.larger[vertex].add(child)
            self.larger[vertex] |= self.larger[child]

    def compare(self, x, y):
        if x == y:
            return 0
        elif x in self.smaller[y]:
            return 1
        elif y in self.smaller[x]:
            return -1
        else:
            return None

    def incomparable(self, x):
        """
        Return the elements which are not comparable to x.
        """
        return self.elements - self.smaller[x] - self.larger[x] - set([x])

    def smallest(self):
        """
        Return the subset of minimal elements.
        """
        return frozenset( [ x for x in self if not self.smaller[x] ] )

    def largest(self):
        """
        Return the subset of maximal elements.
        """
        return frozenset( [ x for x in self if not self.larger[x] ] )

    def closure(self, A):
        """
        Return the smallest set X containing A which is is closed
        under < , i.e. such that (x in X and y < x) => y in X.
        """
        result = frozenset(A)
        for a in A:
            result |= self.smaller[a]
        if len(result) == len(A):
            return result
        else:
            return self.closure(result)

    def closed_subsets(self, start=None):
        """
        Generator for all transitively closed subsets.
        """
        if start is None:
            if self.closed:
                for subset in self.closed:
                    yield subset
                return
            else:
                start = self.smallest()
        complement = self.elements - start
        if start not in self.closed:
            self.closed.add(start)
            yield start
        for element in complement:
            print( 'adding ', element)
            extended = self.closure(start | set([element]))
            for subset in self.closed_subsets(extended):
                yield subset

    def closed_subsets(self, start=None):
        """
        Generator for all transitively closed subsets.  The subsets
        are computed once, then cached for use in subsequent calls.
        """
        if start is None:
            if self.closed:
                for subset in self.closed:
                    yield subset
                return
            else:
                start = self.smallest()
        if start not in self.closed:
            self.closed.add(start)
            yield start
        children = set()
        for element in start:
            children.update(self.successors[element] - start)
        for child in children:
            extended = self.closure(start | set([child]))
            for subset in self.closed_subsets(extended):
                yield subset

class FatGraph(Graph):

    def __init__(self, pairs, singles=[]):
        Graph.__init__(self, pairs, singles)
        self.adjacent = {}
        for v in self.vertices:
            self.adjacent[v] = [e(v) for e in self.edges if v in e]

    def __getitem__(self, vertex):
        return self.adjacent[vertex]

    def __repr__(self):
        return '\n'.join(['%s: %s'%(v, self.adjacent[v])
                          for v in self.vertices])

class FatDigraph(FatGraph):
    edge_class = DirectedEdge


if _within_sage:
    def _to_sage(self):
        S = sage.graphs.graph.Graph(loops=True, multiedges=True)
        S.add_vertices(self.vertices)
        for e in self.edges:
            S.add_edge(e.ends[0], e.ends[1], repr(e))
        return S
    Graph.sage = _to_sage

    def _to_sage_digraph(self):
        S = sage.graphs.graph.DiGraph(loops=True, multiedges=True)
        S.add_vertices(self.vertices)
        for e in self.edges:
            S.add_edge(e.ends[0], e.ends[1], repr(e))
        return S
    Digraph.sage = _to_sage_digraph


