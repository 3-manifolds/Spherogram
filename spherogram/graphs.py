"""
Python implementation of graphs, reduced graphs, directed graphs, fat
graphs (graphs with ordered adjacency lists) and fat directed graphs.

Vertices are arbitrary python objects.  To instantiate a graph,
provide a list of edges, as pairs of vertices, and an optional list of
additional vertices, which may be isolated.  Vertices are saved as a
set, so redundancies in the vertex list will be ignored.  Edges are
saved as a list, so multiple edges are allowed.

If G is a graph and v a vertex of G then G[v] is a list of vertices
adjacent to v.
"""

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
            raise ValueError, 'Vertex is not an endpoint'

    def __repr__(self):
        return '%s --- %s'%tuple(self.ends)

    def __iter__(self):
        return self.ends.__iter__() 
        
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
    and e is the parent of f with respect to the depth-first ordering.
    """
    def __init__(self, graph, source, forbidden=set()):
        self.graph = graph
        self.initial = initial = graph.incident(source) - forbidden
        self.seen = forbidden | initial
        self.fifo = deque([ (None, source, e) for e in initial ])

    def __iter__(self):
        return self
    
    def next(self):
        if self.fifo:
            parent, vertex, child = self.fifo.popleft()
            new_edges = self.graph.incident(child(vertex)) - self.seen
            self.seen |= new_edges
            self.fifo.extend([(child, child(vertex), edge)
                              for edge in new_edges])
            return parent, vertex, child
        else:
            raise StopIteration

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
        return [e(vertex) for e in self.edges if vertex in e]

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
                     if vertex in e and e(vertex) != vertex ])

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
        set of edges that cross the cut, and a maximal family of
        edge-disjoint paths from source to sink.
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
            path = []
            weight = residual[child]
            while True:
                path.insert(0, (vertex, child))
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
            cut_edges |= set([edge for edge in children[vertex]
                              if edge(vertex) not in cut_set])
        return {'set': cut_set, 'edges': cut_edges, 'paths': path_list}

    def is_planar(self):
        """
        Return the planarity.
        """
        return planar(self)

    def merge(self, V1, V2):
        """
        Merge two vertices and remove all edges between them.
        The vertex objects must support the __or__ method
        (e.g. all vertices might be frozensets).
        The two vertices V1 and V2 are replaced by V1|V2.
        >>> G = Graph([(0,1),(1,2),(2,0)]).mergeable()
        >>> F = lambda x: frozenset([x])
        >>> G.merge(F(1),F(2))
        >>> G
        V: frozenset([1, 2]), frozenset([0])
        E: frozenset([0]) --- frozenset([1, 2]), frozenset([1, 2]) --- frozenset([0])
        """
        new_vertex = V1|V2
        if new_vertex in self.vertices:
            raise ValueError, 'Merged vertex already exists!'
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

    def valence(self, vertex):
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
        """
        print 'Not written'

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
