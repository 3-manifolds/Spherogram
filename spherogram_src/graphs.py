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

For most graphs G, G(v) returns the set of edges incident to G.
In the case of a FatGraph, the return value is the ordered list
of incident edges.  (Loops will appear twice).

When called as a function, an edge is the non-trivial involution of
its endpoints; calling it with one endpoint returns the other one.

One can iterate over edges as:

for x, y in G.edges:
  ...

where x and y will be the two endpoints of the edge (tail and head
in the case of directed edges).
"""

try:
    import sage.all
    import sage.graphs.graph
    _within_sage = True
except ImportError:
    from .planarity import planar
    from spherogram.planarity import planar
    _within_sage = False
    
from collections import deque
import operator

class BaseEdge(tuple):
    """
    Base class for edges: a 2-tuple of vertices with extra methods.
    Calling a BaseEdge with one of its vertices returns the other one.
    """

    def __new__(cls, x, y):
        return tuple.__new__(cls, (x,y))
        
    def __call__(self, end):
        """
        Calling an edge with one endpoint returns the other one.
        """
        if end is self[0]:
            return self[1]
        elif end is self[1]:
            return self[0]
        else:
            raise ValueError('Vertex is not an endpoint')

    def __hash__(self):
        return id(self)
    
    def inciden_to(self):
        return list(self)
    
    def is_loop(self):
        return self[0] is self[1]
    
class Edge(BaseEdge):
    """
    An undirected edge.  We allow multiple Edges between the same
    vertices, so the __eq__ operator is overloaded.
    """

    def __repr__(self):
        return '%s --- %s'%self
        
    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

class MultiEdge(BaseEdge):
    """
    An undirected edge.  MultiEdges are equal if they have the
    same vertices.  The multiplicity is initialized to 1.
    """

    def __init__(self, x, y):
        self.multiplicity = 1
        
    def __repr__(self):
        return '%s --%d-- %s'%(self[0], self.multiplicity, self[1])

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return set(self) == set(other)

class DirectedEdge(BaseEdge):
    """
    An Edge with a tail and a head.  The two vertices can be accessed as
    E.tail and E.head
    """
    
    def __repr__(self):
        return '%s --> %s'%self

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    @property
    def head(self):
        return self[1]

    @property
    def tail(self):
        return self[0]

class DirectedMultiEdge(DirectedEdge):
    """
    A DirectedEdge with multiplicity.  DirectedMultiEdges are equal if
    they have the same head and tail.  The multiplicity is initialized to 1.
    """

    def __init__(self, v, w):
        self.multiplicity = 1
        
    def __repr__(self):
        return '%s --%d-> %s'%(self[0], self.multiplicity, self[1])

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return tuple(self) == tuple(other)

class FatEdge(Edge):
    """
    An Edge that knows its place among the edges incident to each
    of its vertices.  Initialize with two pairs (v,n) and (w,m) meaning
    that this edge joins v to w and has index n at v and m at w.
    The parity of the optional integer argument twists determines
    whether the edge is twisted or not.
    """
    def __new__(cls, x, y, twists=0):
        return tuple.__new__(cls, (x[0],y[0]))

    def __init__(self, x, y, twists=0):
        self.slots = [x[1],y[1]]
        self.twisted = True if twists%2 != 0 else False

    def __repr__(self):
        return '%s[%d] -%s- %s[%d]'%(self[0], self.slots[0],
                                     'x' if self.twisted else '-',
                                     self[1], self.slots[1])

    def __hash__(self):
        return id(self)

    def slot(self, vertex):
        try:
            return self.slots[self.index(vertex)]
        except ValueError:
            raise ValueError('Vertex is not an end of this edge.')

    def set_slot(self, vertex, n):
        try:
            self.slots[self.index(vertex)] = n
        except ValueError:
            raise ValueError('Vertex is not an end of this edge.')


class EdgesBFO(object):
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

class Graph(object):
    """
    A set of vertices and a set of edges joining pairs of vertices.
    Vertices are arbitrary hashable objects and should not be
    modified by Graph methods.
    """
    edge_class = Edge

    def __init__(self, pairs=[], singles=[]):
        self.vertices = set()
        self.edges = set()
        self.incidence_dict = {}
        self.Edge = self.__class__.edge_class
        for pair in pairs:
            self.add_edge(*pair)
        for vertex in singles:
            self.add_vertex(vertex)
        self._validate()

    def _validate(self):
        pass
    
    def __repr__(self):
        V = 'Vertices:\n  ' + '\n  '.join([str(v) for v in self.vertices])
        E = 'Edges:\n  ' + '\n  '.join([str(e) for e in self.edges])
        return '%s\n%s'%(V,E)

    def __call__(self, vertex):
        """
        Return the set of incident edges.
        """
        return self.incidence_dict[vertex]
            
    def __getitem__(self, vertex):
        """
        Return a list of adjacent vertices.
        """
        return [e(vertex) for e in self.incidence_dict[vertex]]

    def add_edge(self, *args):
        edge = self.Edge(*args)
        self.edges.add(edge)
        self.vertices.update(edge)
        for v in edge:
            try:
                self.incidence_dict[v].append(edge)
            except KeyError:
                self.incidence_dict[v] = [edge]
        return edge

    def add_vertex(self, hashable):
        self.vertices.add(hashable)
        if not self.incidence_dict.has_key(hashable):
            self.incidence_dict[hashable] = []

    def remove_vertex(self, hashable):
        self.vertices.remove(hashable)
        incident = self.incidence_dict.pop(hashable)
        self.edges -= set(incident)

    def incident(self, vertex):
        """
        Return the set of non-loops incident to the vertex.
        """
        return set(e for e in self.incidence_dict[vertex] if not e.is_loop())

    def valence(self, vertex):
        """
        Return the valence of a vertex.
        """
        return len(self.incidence_dict[vertex])

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
        vertices, result = self.vertices - set(deleted_vertices), []
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

        The edge capacities are supplied as a dictionary, with
        edges as keys and the capacity of an edge as value.  If
        no capacity dict is supplied, every edge is given capacity
        1.
        
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
            flow = residual[child]
            while True:
                path.appendleft( (vertex, child) )
                children[vertex].add(child)
                if vertex == source:
                    break
                child = parent
                parent, vertex = parents[child]
                flow = min(flow, residual[child])
            for vertex, edge in path:
                residual[edge] -= flow
                if residual[edge] == 0:
                    full_edges.add(edge)
            path_list.append( (flow, path) )
        # Find the cut edges.
        cut_edges = set()
        for vertex in cut_set:
            cut_edges |= set([edge for edge in self.edges
                              if vertex in edge
                              and edge(vertex) not in cut_set])
        unsaturated = [ e for e in self.edges if residual[e] > 0 ]
        return {'set': cut_set, 'edges': cut_edges, 'paths': path_list,
                'unsaturated': unsaturated}

    def reduced(self):
        R = ReducedGraph()
        for e in self.edges:
            R.add_edge(*e)
        return R

    def is_planar(self):
        return self.reduced().is_planar()
    
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
    edge_class = MultiEdge
    
    def __init__(self, pairs=[], singles=[]):
        self.Edge = self.__class__.edge_class
        Graph.__init__(self, pairs, singles) 

    def find_edge(self, vertex1, vertex2):
        pair = (vertex1, vertex2)
        for edge in self.edges:
            if edge == pair:
                return edge

    def add_edge(self, x, y):
        edge = self.find_edge(x, y)
        if edge:
            edge.multiplicity += 1
        else:
            edge = self.Edge(x, y)
            self.vertices.update([x,y])
            self.edges.add(edge)
            for v in edge:
                try:
                    self.incidence_dict[v].append(edge)
                except KeyError:
                    self.incidence_dict[v] = [edge]
        return edge

    def multi_valence(self, vertex):
        """
        Return the valence of a vertex, counting edge multiplicities.
        """
        valence = 0
        for e in self.incidence_dict[v]:
                valence += e.multiplicity
        return valence

    def is_planar(self):
        """
        Return the planarity.
        """

        verts_with_loops = set()
        non_loop_edges = set()
        for e in self.edges:
            v, w = e
            if v != w:
                non_loop_edges.add( (v, w) )
            else:
                verts_with_loops.add(v)

        sans_loops = ReducedGraph(non_loop_edges)
        if _within_sage:
            S = sans_loops.sage(loops=False, multiedges=False)
            is_planar = S.is_planar(set_embedding=True)
            embedding = S.get_embedding()
        else:
            is_planar, embedding = planar(sans_loops)

        if is_planar:
            for v in verts_with_loops:
                embedding[v].append(v)
            self._embedding = embedding

        return is_planar

    def embedding(self):
        if self.is_planar():
            return self._embedding

    def one_min_cut(self, source, sink):
        capacity = dict((e, e.multiplicity) for e in self.edges)
        cut = Graph.one_min_cut(self, source, sink, capacity)
        cut['size'] = sum([e.multiplicity for e in cut['edges'] ])
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
            if not edge or edge.multiplicity < 2:
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

class CyclicList(list):
    def __getitem__(self, n):
        return list.__getitem__(self, n%len(self))
    def succ(self, x):
        return self[(self.index(x)+1)%len(self)]
    def pred(self, x):
        return self[(self.index(x)-1)%len(self)]

class FatGraph(Graph):
    """
    A FatGraph is a Graph which maintains a CyclicList of incident
    edges for each vertex.  The edges are FatEdges, which come in two
    flavors: twisted and untwisted.  Since the incident edges are in
    particulary cyclically ordered, a FatGraph has a canonical
    embedding as a spine of a surface with boundary.  However, a
    CyclicList has a first element.  This extra data, namely a choice
    of distinguished edge for each vertex, is used to give a canonical
    mapping from FatGraphs to link or tangle diagrams.

    A FatGraph should be initialized with a sequence of either double
    pairs ((v,n), (w,m)) or triples ((v,n), (w,m), T).  These indicate
    that there is an edge from v to w which has index n at v and m at
    w.  If T is supplied then the parity of T determines the twist.
    """
    edge_class = FatEdge

    def __call__(self, vertex):
        return self.incidence_dict[vertex]

    def add_edge(self, *args):
        incidences = self.incidence_dict
        edge = self.Edge(*args)
        self.edges.add(edge)
        self.vertices.update(edge)
        for v in edge:
            # the values of incidence_dict should be objects that keep
            # themselves sorted.
            if v in incidences:
                incidences[v].append(edge)
                incidences[v].sort(key=lambda e : e.slot(v))
            else:
                incidences[v] = CyclicList([edge])
        return edge

    def _validate(self):
        for v in self.vertices:
            slots = [e.slot(v) for e in self(v)]
            assert slots == range(len(slots))

    def reorder(self, vertex, cyclist):
        for e, n in zip(self(vertex), cyclist):
            e.set_slot(vertex, n)
        self.incidence_dict[vertex].sort(key=lambda e : e.slot(vertex))

    def boundary_cycles(self):
        left  = [(e[0], e, 'L') for e in self.edges]
        right = [(e[0], e, 'R') for e in self.edges]
        sides = left + right
        cycles = []
        while sides:
            cycle = []
            v, e, s = start = sides.pop()
            #print('start:', start)
            while True:
                cycle.append(e)
                # cross the edge
                v = e(v)
                if (
                    ( e.twisted and s=='L') or
                    ( not e.twisted and (s=='L')==(v==e[0]) )
                    ): # go counter-clockwise
                    e = self(v).succ(e)
                    s = 'R' if e.twisted or v == e[0] else 'L'
                else: # go clockwise
                    e = self(v).pred(e)
                    s = 'L' if e.twisted or v == e[0] else 'R'
                if (e[0],e,s) == start:
                    cycles.append(cycle)
                    break
                #print('next', (e[0],e,s))
                sides.remove( (e[0],e,s) )
        return cycles

    def filled_euler(self):
        return len(self.vertices) - len(self.edges) + len(self.boundary_cycles())

class Digraph(Graph):
    edge_class = DirectedEdge

    def incident(self, vertex):
        """
        Return the set of non-loops which *begin* at the vertex.
        """
        return set(e for e in self.incidence_dict[vertex]
                     if e.tail is vertex and not e.head is vertex)

    def incident_to(self, vertex):
        """
        Return the set of non-loops which *end* at the vertex.
        """
        return set(e for e in self.incidence_dict[vertex]
                     if e.head is vertex and not e.tail is vertex)

    def in_valence(self, vertex):
        return len([e for e in self.incidence_dict[vertex] if e.head is vertex])

    def out_valence(self, vertex):
        return len([e for e in self.incidence_dict[vertex] if e.tail is vertex])

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
        
class StrongConnector(object):
    """
    Finds strong components of a digraph using Tarjan's algorithm;
    see http://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
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
        Return the directed acyclic graph whose vertices are the
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
    A partially ordered set, generated from a directed acyclic graph.
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
    
if _within_sage:
    def _to_sage(self, loops=True, multiedges=True):
        S = sage.graphs.graph.Graph(loops=loops, multiedges=multiedges)
        S.add_vertices(self.vertices)
        for e in self.edges:
            v, w = e
            S.add_edge(v, w, repr(e))
        return S
    Graph.sage = _to_sage

    def _to_sage_digraph(self):
        S = sage.graphs.graph.DiGraph(loops=True, multiedges=True)
        S.add_vertices(self.vertices)
        for e in self.edges:
            v, w = e
            S.add_edge(v, w, repr(e))
        return S
    Digraph.sage = _to_sage_digraph
