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
        vertices = 'V: ' + ', '.join([str(v) for v in self.vertices])
        edges = 'E: ' + ', '.join([str(e) for e in self.edges])
        return '%s\n%s'%(vertices, edges)

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
    
    def edges_bfo(self, source, forbidden=set()):
        """
        Generator of embedded edges of the graph that are not in the
        forbidden set, ordered by distance from the source.

        Yields triples (e, v, f) where e and f are edges containing v
        and e is the parent of f with respect to the depth-first ordering.
        """
        initial = self.incident(source) - forbidden
        seen = forbidden | initial
        fifo = deque([ (None, source, e) for e in initial ])
        while fifo:
            parent, vertex, child = fifo.popleft()
            neighbors = self.incident(child(vertex)) - seen
            seen |= neighbors
            fifo.extend([(child, child(vertex), neighbor)
                         for neighbor in neighbors])
            yield parent, vertex, child

    def min_cut(self, source, sink, capacity=None):
        """
        Find a minimal cut which separates source from sink.

        Returns the cut set of vertices on the source side, and the
        set of edges that cross the cut.
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
        while True:
            # Try to find a new path from source to sink
            parents, cut_set, failed = {}, set([source]), True
            for parent, vertex, child in self.edges_bfo(source, full_edges):
                 parents[child] = (parent, vertex)
                 cut_set.add(child(vertex))
                 if child(vertex) == sink:
                     failed = False
                     break
            # If we fail, we have visited every vertex reachable
            # from the source, thereby providing the cut set.
            if failed:
                break
            # On success, adjust the flow, record the children, and go on.
            path = []
            flow = residual[child]
            while True:
                path.insert(0, child)
                children[vertex].add(child)
                if vertex == source:
                    break
                child = parent
                parent, vertex = parents[child]
                flow = min(flow, residual[child])
            for edge in path:
                residual[edge] -= flow
                if residual[edge] == 0:
                    full_edges.add(edge)
            #print path
        #print cut_set
        # Find the cut edges.
        cut_edges = set()
        for vertex in cut_set:
            cut_edges |= set([edge for edge in children[vertex]
                              if edge(vertex) not in cut_set])
        return cut_set, cut_edges

    def is_planar(self):
        return planar(self)
    
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

    def min_cut(self, source, sink):
        cut_set, cut_edges = Graph.min_cut(self, source, sink,
                                           self.multiplicities.copy())
        cut_size = sum([self.multiplicities[e] for e in cut_edges])
        return cut_set, cut_edges, cut_size
        
class Digraph(Graph):
    edge_class = DirectedEdge

    def incident(self, vertex):
        """
        Return the set of non-loops which begin at the vertex.
        """
        return set([ e for e in self.edges
                     if e.tail() == vertex and e.head() != vertex ])

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
