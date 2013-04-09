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
# crossing as South, East, North, West in such a way that the first
# time a component passes through the vertex it enters at South and
# leaves at North; the second time it enters at East and leaves at
# West.
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
#
# This determines an embedding of each crossing into the oriented
# plane, which may not extend to an embedding of the knot diagram.  In
# constructing the planar embedding of the diagram we preserve this
# embedding# of the vertices, but adjust the fat graph by
# interchanging the attaching points of a pair of edges entering at
# opposite sides of the crossing.

South, East, North, West = 0, 1, 2, 3

class FlippingError(Exception):
    pass

class DTvertex:
    """
    A vertex of the 4-valent graph which is described by a DT code.
    Instantiate with an even-odd pair, in either order.
    """
    # In keeping with the Spherogram graphs philosophy, vertices are
    # immutable.  All information is stored in the edges

    def __init__(self, pair, overcrossing=1):
        self._first = min(pair)
        self._second = max(pair)
        self._even_over = True if overcrossing == -1 else False

    def __repr__(self):
        return str((self._first, self._second))

    def enter(self, N):
        if N == self._first: return South
        elif N == self._second: return East
        else: raise ValueError('%d is not a label of this vertex'%N)

    def exit(self, N):
        if N == self._first: return North
        elif N == self._second: return West
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
    """
    A fat edge which can be marked.
    """
    def __init__(self, x, y, twists=0):
        FatEdge.__init__(self, x, y, twists)
        self.marked = False

class DTFatGraph(FatGraph):
    edge_class = DTFatEdge

    def path(self, vertex, edge):
        """
        Iterates through the edges of a component, starting at the
        given edge, in the direction determined by the vertex.
        """
        if not vertex in edge:
            raise ValueError('That vertex is not an endpoint of the edge.')
        forward = True if vertex == edge[0] else False
        return DTPath(edge, self, forward)

    def _outside_slots(self, edge, side):
        """
        Assume that the marked subFatGraph has been embedded in the
        plane.  This generator starts at a marked FatEdge and walks
        around one of its adjacent boundary curves (left=-1, right=1),
        yielding all of the pairs (v, s) where s is a slot of the
        vertex v which lies on the specified boundary curve.  To
        extend the embedding over an unmarked arc, the ending slots of
        both ends of the arc must lie on the same boundary curve.
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
                    # For lookups, slot values must be in 0,1,2,3
                    yield vertex, slot%4  
                else:
                    break
            if edge == interior_edge:
                raise ValueError('Found a dead end in the marked subgraph.')
            edge = interior_edge
            vertex = edge(vertex)
            if vertex == first_vertex:
                break

    def left_slots(self, edge):
        return self._outside_slots(edge, side=-1)

    def right_slots(self, edge):
        return self._outside_slots(edge, side=1)

    def unmarked_edges(self, vertex):
        return [e for e in self(vertex) if not e.marked]

    def num_unmarked_edges(self, vertex):
        return len([e for e in self(vertex) if not e.marked])

    def clear(self):
        for e in self.edges:
            e.marked = False
 
    def flip(self, vertex, slot):
        """
        Move the edge at this slot to the opposite slot, and
        move the edge in the opposite slot to this slot.
        """
        um = self.num_unmarked_edges(vertex)
        print 'Flipping %s[%s]'%(vertex, slot)
        if um == 1:
            msg = 'Trying to flip %s which has only 1 unmarked edge.'%vertex
            print msg
            raise FlippingError(msg)
        if slot == 0 or slot == 2:
            self.reorder(vertex, (2,1,0,3))
        elif slot == 1 or slot == 3:
            self.reorder(vertex, (0,3,2,1))
        else:
            raise ValueError('Invalid slot.')

class DTcode:
    """
    Represents a DTcode of a link projection.  Instantiate either from
    a list of tuples or an alphabetical code.
    """
    def __init__(self, code):
        if isinstance(code,str):
            code = self.convert_alpha(code)
        # Here we unpack the numerical code
        overcrossings = [sign(x) for comp in code for x in comp]
        evens = [abs(x) for comp in code for x in comp]
        self.size = size = 2*len(evens)
        pairs = zip(range(1, size, 2), evens)
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
            N += 1
            V = W
            start = N

    def __getitem__(self, n):
        """
        A DTcode can look up a vertex by either of its labels.
        """
        return self.lookup[n]

    def find_circle(self, first_edge):
        """
        Follow a component, starting at the given (directed) edge,
        until the first time it crosses itself.  Throw away the tail
        to get a cycle.  Return the list of vertices and the list of
        edges traversed by the cycle.
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

    def find_unmarked_arc(self, vertex):
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
        while True:
            if not self.embed_arc():
                break

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
        print 'No vertices with only 1 unmarked edge!'
        for v in G.vertices:
            if G.num_unmarked_edges(v) == 2:
                return v
        return None

    def do_flips(self, v, v_edge, w, w_edge):
        """
        Decide whether v and/or w needs to be flipped in order to add
        an arc from v to w starting with the v_edge and ending with
        the w_edge.  If flips are needed, make them.  If the embedding
        cannot be extended raise an exception.
        """
        G = self.fat_graph
        vslot = G(v).index(v_edge)
        wslot = G(w).index(w_edge)
        # starting from the v_edge, go ccw to a marked edge
        for k in range(1,3):
            ccw_edge = G(v)[vslot+k]
            if ccw_edge.marked:
                break
        if not ccw_edge.marked:
            raise ValueError('Invalid marking')
        # first look for w on the same side as the v slot
        if v == ccw_edge[0]: # the ccw_edge points out of v
            slot_bdry = G.right_slots(ccw_edge)
        else:                # the ccw_edge points into v
            slot_bdry = G.left_slots(ccw_edge)
        slots = [slot for vertex, slot in slot_bdry if vertex == w]
        if wslot in slots:
            return
        elif slots:
            G.flip(w, wslot)
            return
        # if that fails, look for w on the other side
        if v == ccw_edge[0]: # the ccw_edge points out of v
            other_bdry = G.left_slots(ccw_edge)
        else:                # the ccw_edge points into v
            other_bdry = G.right_slots(ccw_edge)
        slots = [slot for vertex, slot in other_bdry if vertex == w]
        if slots:
            G.flip(v, vslot)
            if wslot in slots:
                return
            else:
                G.flip(w, wslot)
                return
        raise FlippingError('Could not find any flips to do.')

    def embed_arc(self):
        G = self.fat_graph
        # get a vertex, preferably with 3 marked edges.
        v = self.get_incomplete_vertex()
        if v is None:
            return False
        arc_edges, last_vertex = self.find_unmarked_arc(v)
        print 'Arc from %s to %s.'%(v,last_vertex)
        print arc_edges
        # We may need to try from both ends of the arc.  We may not be
        # able to flip one end because it only has 1 unmarked
        # edge. But sometimes one end can not see the other end's
        # slots while the second one can see the slots of the first:
        #   ______ ______
        #  |      |      |
        #  |      |/v    |/w
        #  |      |\     |\
        #  |______|______|   
        #
        # Since v is more likely to have 3 marked edges, try it second.
        try:
            self.do_flips(last_vertex, arc_edges[-1], v, arc_edges[0])
        except FlippingError:
            print 'Trying the other end'
            self.do_flips(v, arc_edges[0], last_vertex, arc_edges[-1])
        self.mark(arc_edges)
        return True
