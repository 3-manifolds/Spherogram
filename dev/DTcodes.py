from spherogram import FatGraph, FatEdge

def sign(x):
    return 1 if x > 0 else -1 if x < 0 else 0

# To construct a knot projection from a DT code we first construct
# a fat graph, which may not be a planar surface.  There are only
# two possible orderings of edges at each vertex, since we know which
# pairs of edges are opposites.  The process of finding the
# projection consists of reversing the orderings of some vertices
# to get a planar surface.
#
# Each crossing in the diagram is traversed twice in building the DT
# code.  Our convention is to initially label the four edges at the
# crossing as South, East, North, West so that the first time a
# component passes through the vertex it enters at South and leaves at
# North; the second time it enters at East and leaves at West.
#
#
#        N
#        ^
#        |
#        |
# W <--------- E second         Initial vertex orientation 
#        |
#        |
#        S
#      first

South, East, North, West = 0, 1, 2, 3

class DTvertex:
    """
    A vertex of the 4-valent graph which is described by a DT code.
    Instantiate with an even-odd pair, in either order.
    """
    def __init__(self, pair, overcrossing=1):
        x, y = pair
        self.first = x if x < y else y
        self.second = y if x < y else x
        self.even = x if y%2 else y
        self.odd = y if y%2 else x
        self.even_over = True if overcrossing == -1 else False

    def __repr__(self):
        return str((self.first, self.second))

    def enter(self, N):
        if N == self.first: return South
        elif N == self.second: return East
        else: raise ValueError('%d is not a label of this vertex'%N)

    def exit(self, N):
        if N == self.first: return North
        elif N == self.second: return West
        else: raise ValueError('%d is not a label of this vertex'%N)

class DTPath:
    """
    An iterator which starts at a FatEdge and walks around the
    link component containing that edge.  A DTPath raises
    StopIteration when it returns to its starting edge.
    """
    def __init__(self, edge, graph, forward=True):
        self.first_edge = edge
        self.next_edge = None
        self.graph = graph
        self.end = 1 if forward else 0

    def __iter__(self):
        return self

    def next(self):
        try:
            slot = self.next_edge.slots[self.end]
            self.next_edge = self.graph(self.next_edge[self.end])[slot+2]
            if self.next_edge == self.first_edge:
                raise StopIteration
        except AttributeError:
            self.next_edge = self.first_edge
        return self.next_edge

class DTFatEdge(FatEdge):
    def __init__(self, x, y, twists=0):
        FatEdge.__init__(self, x, y, twists)
        self.marked = False

class DTFatGraph(FatGraph):
    edge_class = DTFatEdge

    def path(self, vertex, edge):
        """
        Iterates through the component containing the edge.
        """
        if vertex == edge[0]:
            forward = True
        elif vertex == edge[1]:
            forward = False
        else:
            raise ValueError('Vertex is not an endpoint.')
        return DTPath(edge, self, forward)

    def _outside_edges(self, edge, side):
        """
        Assume that the marked subFatGraph has been embedded in the
        plane.  This generator starts at a marked FatEdge and walks
        around one of its adjacent boundary curves (left=-1, right=1),
        yielding all of the (v,e) pairs where e is an unmarked edge
        adjacent to v that, if directed away from v, would point into the
        disk bounded by that curve if the embedding were extended to
        the entire FatGraph.
        """
        if not edge.marked:
            raise ValueError('Must begin at a marked edge.')
        first_vertex = vertex = edge[1] 
        while True:
            end = 0 if edge[0] == vertex else 1
            slot = edge.slots[end]
            for k in range(3):
                slot += side
                interior_edge = self(vertex)[slot]
                if not interior_edge.marked:
                    yield vertex, interior_edge
                else:
                    break
            if edge == interior_edge:
                raise ValueError('Dead end in marked subgraph.')
            else:
                edge = interior_edge
                vertex = edge(vertex)
            if vertex == first_vertex:
                break

    def left_edges(self, edge):
        return self._outside_edges(edge, side=-1)

    def right_edges(self, edge):
        return self._outside_edges(edge, side=1)

    def unmarked_edges(self, vertex):
        return [e for e in self(vertex) if not e.marked]

    def num_unmarked_edges(self, vertex):
        return len([e for e in self(vertex) if not e.marked])

    def clear(self):
        for e in self.edges:
            e.marked = False

    def flip(self, vertex, fixed=0):
        if fixed == 0 or fixed == 2:
            self.reorder(vertex, (0,3,2,1))
        elif fixed == 1 or fixed == 3:
            self.reorder(vertex, (2,1,0,3))
        else:
            raise ValueError('Invalid fixed direction.')

class DTcode:
    """
    Represents a DTcode of a link projection.  Instantiate either from
    a list of tuples or an alphabetical code.
    """
    def __init__(self, code):
        if isinstance(code,str):
            code = self.convert_alpha(code)
        overcrossings = [sign(x) for comp in code for x in comp]
        evens = [abs(x) for comp in code for x in comp]
        self.size = size = 2*len(evens)
        odds =  range(1, size, 2)
        self.pairs = pairs = zip(odds, evens)
        # Build a lookup table for vertices.
        # (DT codes are 1-based; we just waste the 0 entry.) 
        self.lookup = lookup = [None for n in range(1+size)]
        for pair, overcrossing in zip(pairs, overcrossings):
            m, n = pair
            V = DTvertex(pair, overcrossing)
            lookup[m] = lookup[n] = V
        # Now build the fatgraph determined by the DT code.
        self.fat_graph = DTFatGraph()
        N = start = 1
        last_odd = -1
        V = self[1]
        for component in code:
            last_odd += 2*len(component)
            # Walk around this component, adding edges.
            while N <= last_odd:
                W = self[N + 1]
                self.fat_graph.add_edge((V, V.exit(N)),
                                        (W, W.enter(N+1)))
                N += 1
                V = W
            # Close up this component and go to the next one.
            S = self[start]
            self.fat_graph.add_edge( (V, V.exit(N)),
                                     (S, S.enter(start)) )
            start = N
            V = W

    def __getitem__(self, n):
        """
        We can look up a vertex by either of its labels
        """
        return self.lookup[n]

    def find_circle(self, first_edge):
        """
        Follow a path, starting at the given (directed) edge, until the
        first time it crosses itself.  Mark each edge in the circle.
        Return the list of vertices and the list of edges traversed
        by the cycle.
        """
        edges = []
        vertices = [first_edge[0]]
        seen = set(vertices)
        for edge in self.fat_graph.path(first_edge[0], first_edge):
            vertex = edge[1]
            if vertex in seen:
                edges.append(edge)
                break
            else:
                seen.add(vertex)
                vertices.append(vertex)
                edges.append(edge)
        n = vertices.index(vertex)
        edges = edges[n:]
        vertices = vertices[n:]
        return vertices, edges

    def find_arc(self, vertex):
        """
        Starting at this vertex, find an unmarked edge and follow its
        component until it hits a vertex with at least one marked
        edge.  Remove loops to get an embedded arc. Return the list
        of edges traversed by the embedded arc.
        """
        G = self.fat_graph
        unmarked = G.unmarked_edges(vertex)
        if len(unmarked) == 0:
            raise ValueError('Vertex must have unmarked edges.')
        if len(unmarked) == 4:
            raise ValueError('Vertex must lie in the subgraph.')
        edges = []
        vertices = []
        seen = set()
        for edge in G.path(vertex, unmarked[0]):
            edges.append(edge)
            vertex = edge(vertex)
            if G.num_unmarked_edges(vertex) < 4:
                break
            if vertex in seen:
                n = vertices.index(vertex)
                edges = edges[:n+1]
                vertices = vertices[:n+1]
                seen = set(vertices)
            else:
                vertices.append(vertex)
                seen.add(vertex)
        return edges, vertex

    def mark(self, edgelist):
        for edge in edgelist:
            edge.marked = True

    def embed(self, edge=None):
        """
        Run the Dowker-Thistlethwaite algorithm to flip crossings
        in the FatGraph until it becomes planar.
        """
        G = self.fat_graph
        if edge is None: # just pick one
            for edge in G.edges: break
        vertices, edges = self.find_circle(edge)
        print 'circle'
        print vertices
        print edges
        self.mark(edges)
        # and so on ...

    def get_incomplete_vertex(self):
        """
        Return an incomplete vertex, preferably with only one unmarked edge. 
        """
        # This could be a table lookup.
        G = self.fat_graph
        # try to find a vertex with 1 unmarked edge
        for v in G.vertices:
            if G.num_unmarked_edges(v) == 1:
                return v
        for v in G.vertices:
            if G.num_unmarked_edges(v) == 2:
                return v
        raise RuntimeError('No more vertices.')

    def flip_test(self, v, w, edge):
        """
        Decide whether v needs to be flipped in order to add an arc from
        v to w starting with the edge.  If so, flip it.  If no such arc
        can be added, raise an exception.
        """
        G = self.fat_graph
        # from the given edge, go ccw to a marked edge
        slot = G(v).index(edge)
        for k in range(1,3):
            ccw_edge = G(v)[slot+k]
            if ccw_edge.marked:
                break
        if not ccw_edge.marked:
            raise ValueError('Invalid marking')
        print 'check', ccw_edge
        print 'right', list(G.right_edges(ccw_edge))
        print 'left', list(G.left_edges(ccw_edge))
        test_right = (v, edge) in G.right_edges(ccw_edge)
        test_left = (v, edge) in G.left_edges(ccw_edge)
        print 'left:', test_left, 'right:', test_right
        if not test_left and not test_right:
            raise RuntimeError('DT code is not realizable')
        print G.num_unmarked_edges(v), 'unmarked edges'
        if ( ( v == ccw_edge[0] and not test_right )
             or
             ( v == ccw_edge[1] and not test_left )
             ): 
            if G.num_unmarked_edges(v) == 2:
                print 'flipping %s'%v
                G.flip(v)
            elif G.num_unmarked_edges(w) == 2:
                print 'flipping %s'%w
                G.flip(w)
            else:
                raise RuntimeError('DT code is not realizable')

    def embed_arc(self):
        G = self.fat_graph
        # find a vertex with unmarked edges, preferably only 1
        v = self.get_incomplete_vertex()
        print 'working on %s'%v
        # find an arc from v to the marked subgraph
        arc_edges, last_vertex = self.find_arc(v)
        print arc_edges
        print last_vertex
        # decide if one of v or w needs to be flipped
        self.flip_test(v, last_vertex, arc_edges[0])
        self.mark(arc_edges)
        
