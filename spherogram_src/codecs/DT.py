from ..graphs import FatGraph, FatEdge
from .. import Link, Crossing
from ..links.links import CrossingEntryPoint
from ..links.ordered_set import OrderedSet
from .Base64LikeDT import (decode_base64_like_DT_code, encode_base64_like_DT_code)


def sign(x):
    return 1 if x > 0 else -1 if x < 0 else 0


def char_to_int(x):
    n = ord(x)
    if 96 < n < 123:
        return n - 96
    if 64 < n < 91:
        return 64 - n
    raise ValueError('Not an ascii letter.')


def string_to_ints(s):
    return [char_to_int(x) for x in s]


def partition_list(L, parts):
    assert sum(parts) == len(L)
    ans = []
    k = 0
    for p in parts:
        ans.append(tuple(L[k:k + p]))
        k += p
    return ans


DT_alphabet = '_abcdefghijklmnopqrstuvwxyzZYXWVUTSRQPONMLKJIHGFEDCBA'

# To reconstruct a knot projection from a DT code we first construct
# a fat graph, which may not be a planar surface.  There are only
# two possible orderings of edges at each vertex, since we know which
# pairs of edges are opposites.  The process of finding the
# projection consists of changing the orderings of some vertices
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
# This arbitrary choice determines some embedding of each crossing
# into the oriented plane.  During the process of constructing a
# planar embedding of the diagram we preserve this planar embedding of
# the vertices, but we may adjust the fat graph by interchanging the
# attaching points of a pair of edges attached at opposite sides of
# the crossing.
#
# In fact, interchanging either pair of opposite connections is
# equivalent to interchanging the other pair.  So in practice we
# only interchange the North-South pair.

South, East, North, West = 0, 1, 2, 3


class FlippingError(Exception):
    pass


class EmbeddingError(Exception):
    pass


class DTvertex(tuple):
    """
    A vertex of the 4-valent graph which is described by a DT code.
    Instantiate with an even-odd pair, in either order.
    """
    # In keeping with the philosophy of Spherogram.graphs, vertices
    # should never be changed by graph methods.  The DTcodec can
    # do whatever it wants with them, but in this implementation
    # it never changes vertices either.

    def __new__(self, pair, overcrossing=1):
        even_over = bool(overcrossing == -1)
        return tuple.__new__(self, (min(pair), max(pair), even_over))

    def __repr__(self):
        return str((self[0], self[1]))

    def entry_slot(self, N):
        if N == self[0]:
            return South
        elif N == self[1]:
            return East
        else:
            raise ValueError('%d is not a label of %s' % (N, self))

    def exit_slot(self, N):
        if N == self[0]:
            return North
        elif N == self[1]:
            return West
        else:
            raise ValueError('%d is not a label of %s' % (N, self))

    def upper_pair(self):
        first, second, even_over = self
        return (0, 2) if bool(first % 2) ^ even_over else (1, 3)


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
            self.next_edge = self.graph(self.next_edge[self.end])[slot + 2]
            if self.next_edge == self.first_edge:
                raise StopIteration
        except AttributeError:
            self.next_edge = self.first_edge
        return self.next_edge

    __next__ = next


class DTFatEdge(FatEdge):
    """
    A fat edge which can be marked and belongs to a link component.
    """
    def __init__(self, x, y, twists=0, component=0):
        FatEdge.__init__(self, x, y, twists)
        self.marked = False
        self.component = component

    def PD_index(self):
        """
        The labelling of vertices when building a DT code also
        determines a labelling of the edges, which is needed
        for generating a PD description of the diagram.
        This method returns the edge label.
        """
        v = self[0]
        if self.slot(v) % 2 == 0:
            return v[0]
        else:
            return v[1]


class DTFatGraph(FatGraph):
    edge_class = DTFatEdge

    def __init__(self):
        self.vertices = OrderedSet()
        self.edges = OrderedSet()
        self.incidence_dict = {}
        self.Edge = self.__class__.edge_class
        self.marked_valences = {v: 0 for v in self.vertices}
        self.stack = []
        self.pushed = False

    def add_edge(self, x, y):
        # Adds keys to the marked_valences dict as vertices are added.
        # This will cause trouble if edges are added while edges are marked!
        edge = FatGraph.add_edge(self, x, y)
        self.marked_valences[x[0]] = 0
        self.marked_valences[y[0]] = 0
        return edge

    def mark(self, edgelist):
        """
        Mark all edges in the list.
        """
        vertices = set()
        for edge in edgelist:
            edge.marked = True
            vertices.update(edge)
        for v in vertices:
            self.marked_valences[v] = self.marked_valence(v)

    def marked_valence(self, vertex):
        """
        Compute the marked valence of a vertex.
        """
        valence = 0
        for e in self.incidence_dict[vertex]:
            if e.marked:
                valence += 1
        return valence

    def clear(self):
        """
        Remove all edge markings.
        """
        for e in self.edges:
            e.marked = False
        self.marked_valences = {v: 0 for v in self.vertices}

    def push(self, flips):
        """
        Save the state of this DTFatGraph before doing the flips.
        The flips will be done by the pop.
        """
        # Ignore the first push -- the first arc always is ambiguous
        # but the choice is irrelevant, up to reflecting the plane.
        if not self.pushed:
            self.pushed = True
            return
        self.stack.append(
            [flips,
             [((e[0], e.slots[0]), (e[1], e.slots[1]), e.component)
              for e in self.edges if e.marked],
             [((e[0], e.slots[0]), (e[1], e.slots[1]), e.component)
              for e in self.edges if not e.marked]])

    def pop(self):
        """
        Restore the state of this DTFatGraph and perform the saved flips.
        """
        self.edges = set()
        flips, marked, unmarked = self.stack.pop()
        for x, y, component in marked:
            edge = self.add_edge(x, y)
            edge.component = component
            edge.marked = True
        for x, y, component in unmarked:
            edge = self.add_edge(x, y)
            edge.component = component
        return flips

    def clear_stack(self):
        """
        Reset the state stack.
        """
        self.stack = []
        self.pushed = False

    def path(self, vertex, edge):
        """
        Return an iterator which iterates through the edges of a
        component, starting at the given edge, in the direction
        determined by the vertex.
        """
        if vertex not in edge:
            raise ValueError('That vertex is not an endpoint of the edge.')
        forward = bool(vertex == edge[0])
        return DTPath(edge, self, forward)

    def marked_arc(self, vertex):
        """
        Given a vertex with marked valence 2, find the maximal marked
        arc containing the vertex for which all interior edges have
        marked valence 2.  If the marked subgraph is a circle, or a
        dead end is reached, raise :class:`ValueError`.  Return a list of
        edges in the arc.
        """
        left_path, right_path, vertices = [], [], set()
        vertices.add(vertex)
        try:
            left, right = (e for e in self(vertex) if e.marked)
        except ValueError:
            raise RuntimeError('Vertex must have two marked edges.')
        for edge, path in (left, left_path), (right, right_path):
            V = vertex
            while True:
                path.append(edge)
                V = edge(V)
                if V == vertex:
                    raise ValueError('Marked graph is a circle')
                edges = [e for e in self(V) if e.marked and e != edge]
                if len(edges) == 0:
                    raise ValueError('Marked graph has a dead end at %s.' % V)
                if len(edges) > 1:
                    break
                else:
                    vertices.add(V)
                    edge = edges.pop()
        left_path.reverse()
        return left_path + right_path

    def unmarked_arc(self, vertex):
        """
        Starting at this vertex, find an unmarked edge and follow its
        component until we run into a vertex with at least one marked
        edge.  Remove loops to get an embedded arc. Return the list
        of edges traversed by the embedded arc.
        """
        valence = self.marked_valence(vertex)
        if valence == 4:
            raise ValueError('Vertex must have unmarked edges.')
        if valence == 0:
            raise ValueError('Vertex must be in the marked subgraph.')
        edges, vertices, seen = [], [], set()
        for first_edge in self(vertex):
            if not first_edge.marked:
                break
        for edge in self.path(vertex, first_edge):
            edges.append(edge)
            vertex = edge(vertex)
            if self.marked_valence(vertex) > 0:
                break
            # Remove loops as they appear
            if vertex in seen:
                n = vertices.index(vertex)
                edges = edges[:n + 1]
                vertices = vertices[:n + 1]
                seen = set(vertices)
            else:
                vertices.append(vertex)
                seen.add(vertex)
        return edges, vertex

    def bridge(self, marked_arc):
        """
        Try to find an embedded path of unmarked edges joining a
        vertex in the given arc to a vertex of the marked subgraph
        which lies in the complement of the interior of the arc.  This
        uses a depth-first search, and raises :class:`ValueError` on failure.
        Returns a triple (first vertex, last vertex, edge path).

        Suppose the marked subgraph has no vertices with marked
        valence 3 and has a unique planar embedding.  Choose a maximal
        arc with all interior vertices having marked valence 2. Then
        adding a bridge from that arc produces a new graph with a
        unique planar embedding.

        In the case where the diagram is prime, a bridge will exist
        for *some* vertex on the arc.  But, e.g. in the case where the
        arc snakes through a twist region, it will only exist for the
        two extremal vertices in the arc.
        """
        # Find the interior vertices on the arc.
        e0, e1 = marked_arc[:2]
        v = e0[0] if e0[1] in e1 else e0[1]
        vertex_list = []
        for edge in marked_arc[:-1]:
            v = edge(v)
            vertex_list.append(v)
        for start_vertex in vertex_list:
            # the vertex_set gets expanded to include vertices visited
            # by the bridge path.
            vertex = start_vertex
            vertex_set = set(vertex_list)
            edge_path, vertex_path, seen_edges = [], [], set()
            while True:
                edges = [e for e in self(vertex) if not e.marked
                         and e not in seen_edges
                         and e(vertex) not in vertex_set]
                try:
                    new_edge = edges.pop()
                    vertex = new_edge(vertex)
                    edge_path.append(new_edge)
                    if self.marked_valence(vertex) > 0:
                        return start_vertex, vertex, edge_path
                    seen_edges.add(new_edge)
                    vertex_path.append(vertex)
                    vertex_set.add(vertex)
                except IndexError:  # Cannot continue: back up.
                    if len(edge_path) == 0:
                        break
                    edge_path.pop()
                    vertex_set.remove(vertex_path.pop())
                    try:
                        vertex = vertex_path[-1]
                    except IndexError:
                        vertex = start_vertex
        raise ValueError('Could not find a bridge.')

    def _boundary_slots(self, edge, side):
        """
        Assume that the marked subFatGraph has been embedded in the
        plane.  This generator starts at a marked FatEdge and walks
        around one of its adjacent boundary curves (left=-1, right=1),
        yielding all of the pairs (v, s) where s is a slot of the
        vertex v which lies on the specified boundary curve, or
        (v, None) if none of the slots at v lie on the curve.  (To
        extend the embedding over an unmarked arc, the ending slots of
        both ends of the arc must lie on the same boundary curve.
        Flipping may be needed to arrange this.)
        """
        if not edge.marked:
            raise ValueError('Must begin at a marked edge.')
        first_vertex = vertex = edge[1]
        while True:
            end = 0 if edge[0] is vertex else 1
            slot = edge.slots[end]
            for k in range(3):
                slot += side
                interior_edge = self(vertex)[slot]
                if not interior_edge.marked:
                    # For lookups, slot values must be in 0,1,2,3
                    yield vertex, slot % 4
                else:
                    break
            if k == 0:
                yield (vertex, None)
            if edge is interior_edge:
                raise ValueError('Marked subgraph has a dead end.')
            edge = interior_edge
            vertex = edge(vertex)
            if vertex is first_vertex:
                break

    def left_slots(self, edge):
        """
        Return the (vertex, slot) pairs on the left boundary curve.
        """
        return set(self._boundary_slots(edge, side=-1))

    def right_slots(self, edge):
        """
        Return the (vertex, slot) pairs on the right boundary curve.
        """
        return set(self._boundary_slots(edge, side=1))

    def flip(self, vertex, force=False):
        """
        Move the edge at the North slot to the South slot, and
        move the edge in the South  slot to the North slot.
        """
        if not force and self.marked_valences[vertex] > 2:
            msg = 'Cannot flip %s with marked valence %d.' % (
                vertex, self.marked_valences[vertex])
            raise FlippingError(msg)
        self.reorder(vertex, (North, East, South, West))

    def incoming_under(self, vertex):
        first, second, even_over = vertex
        incoming = [e.PD_index() for e in self(vertex) if e[1] is vertex]
        incoming.sort(key=lambda x: x % 2)
        return incoming[0] if even_over else incoming[1]

    def PD_tuple(self, vertex):
        """
        Return the PD labels of the incident edges in order, starting
        with the incoming undercrossing as required for PD codes.
        """
        edgelist = [e.PD_index() for e in self(vertex)]
        n = edgelist.index(self.incoming_under(vertex))
        return tuple(edgelist[n:] + edgelist[:n])

    def flipped(self, vertex):
        """
        Has this vertex been flipped?
        """
        return bool(len([e for e in self(vertex)
                         if e[1] is vertex and e.slots[1] in (2, 3)]) % 2)

    def sign(self, vertex):
        """
        The sign of the crossing corresponding to this vertex.
        See the documentation for Spherogram.link.
        """
        flipped = self.flipped(vertex)
        even_first = bool(vertex[0] % 2 == 0)
        return -1 if (flipped ^ vertex[2] ^ even_first) else 1

    def KLP_strand(self, vertex, edge):
        """
        Return the SnapPea KLP strand name for the given edge at the
        end opposite to the vertex.
        """
        W = edge(vertex)
        slot = edge.slot(W)
        return 'X' if (slot == 0 or slot == 2) ^ self.flipped(W) else 'Y'

    def KLP_dict(self, vertex, indices):
        """
        Return a dict describing this vertex and its neighbors
        in KLP terminology.

        The translation from our convention is as follows::

                    Y                    Y
                    3                    0
                    ^                    ^
                    |                    |
             0 -----+----> 2 X     1 ----+---> 3 X
                    |                    |
                    |                    |
                    1                    2
               not flipped           flipped

        The indices argument is a dict that assigns an integer
        index to each vertex of the graph.
        """
        KLP = {}
        flipped = self.flipped(vertex)
        edges = self(vertex)
        neighbors = self[vertex]
        strands = [self.KLP_strand(vertex, edge) for edge in edges]
        ids = [indices[v] for v in neighbors]

        KLP['sign'] = 'R' if self.sign(vertex) == 1 else 'L'
        slot = 1 if flipped else 0
        KLP['Xbackward_neighbor'] = ids[slot]
        KLP['Xbackward_strand'] = strands[slot]
        slot = 3 if flipped else 2
        KLP['Xforward_neighbor'] = ids[slot]
        KLP['Xforward_strand'] = strands[slot]
        KLP['Xcomponent'] = edges[slot].component
        slot = 2 if flipped else 1
        KLP['Ybackward_neighbor'] = ids[slot]
        KLP['Ybackward_strand'] = strands[slot]
        slot = 0 if flipped else 3
        KLP['Yforward_neighbor'] = ids[slot]
        KLP['Yforward_strand'] = strands[slot]
        KLP['Ycomponent'] = edges[slot].component
        return KLP

# This assumes that the diagram has no loops, and that each component
# meets the next one (so in particular the diagram is connected.


class DTcodec:
    """
    Codec for DT codes of a link projection.
    """

    def __init__(self, input, flips=None):
        if isinstance(input, (bytes, str, list)):
            self.decode(input, flips)

    def __getitem__(self, n):
        """
        A DTcodec can look up a vertex by either of its labels.
        """
        return self.lookup[n]

    def decode(self, dt, flips=None):
        """
        Accepts input of the following types:
        1) a dt code in either numeric or alphabetical form and a sequence
        of boolean values, one for each successive label in the DT code,
        indicating whether the crossing with that label needs to be
        flipped.  Alphabetical DT codes can also specify signs by appending
        a period followed by a sequence of 0's and 1's.
        2) a DT code with flips set to None.  In this case the flips are
        computed.
        3) a bytes object containing a compact signed DT code.  The
        signed DT code may be provided as a byte string, for which
        the last byte has bit 7 set, as produced by the signed_DT method.
        Alternatives, it can be hex encoded as a string
        beginning with '0x', as produced by the hex_signed_DT method.

        This method constructs a planar FatGraph from its input data.
        """
        if flips is not None:
            flips = [bool(flip) for flip in flips]
        self.flips = flips
        if isinstance(dt, (str,)):
            if dt[:2] == '0x':
                dt_bytes = [int(dt[n:n + 2], 16) for n in range(2, len(dt), 2)]
                self.code, self.flips = self.unpack_signed_DT(dt_bytes)
            elif ord(dt[-1]) & 1 << 7:
                dt_bytes = bytearray(dt)
                self.code, self.flips = self.unpack_signed_DT(dt)
            elif dt[0] in '123456789':
                self.code, self.flips = decode_base64_like_DT_code(dt)
            else:
                parts = dt.split('.')
                self.code = self.convert_alpha(parts[0])
                if len(parts) > 1:
                    self.flips = [d != '0' for d in parts[1]]
        elif isinstance(dt, bytes):
            self.code, self.flips = self.unpack_signed_DT(dt)
        else:
            self.code = dt
        code = self.code
        overcrossings = (sign(x) for comp in code for x in comp)
        evens = [abs(x) for comp in code for x in comp]
        self.size = size = 2 * len(evens)
        pairs = zip(range(1, size, 2), evens)
        # Build a lookup table for the vertices.
        # (DT codes are 1-based; we just waste the 0 entry.)
        self.lookup = lookup = [None for n in range(1 + size)]
        for pair, overcrossing in zip(pairs, overcrossings):
            V = DTvertex(pair, overcrossing)
            m, n = pair
            lookup[m] = lookup[n] = V
        # Now build the fatgraph determined by the DT code.
        self.fat_graph = G = DTFatGraph()
        N = start = 1
        last_odd = -1
        for c, component in enumerate(code):
            if len(component) == 0:
                continue
            last_odd += 2 * len(component)
            V = self[N]
            # Walk around this component, adding edges.
            while N <= last_odd:
                W = self[N + 1]
                edge = G.add_edge((V, V.exit_slot(N)),
                                  (W, W.entry_slot(N + 1)))
                edge.component = c
                N += 1
                V = W
            # Close up this component and go to the next one.
            S = self[start]
            edge = G.add_edge((V, V.exit_slot(N)),
                              (S, S.entry_slot(start)))
            edge.component = c
            start = N = N + 1
        # Now build the planar embedding
        labels = [abs(N) for component in code for N in component]
        if self.flips is None:
            self.embed()
            # Our convention is that the first crossing is positive.
            if G.sign(self[1]) != 1:
                for label in labels:
                    G.flip(self[label], force=True)
            self.flips = [G.flipped(self[label]) for label in labels]
        else:
            for label, flip in zip(labels, self.flips):
                if flip:
                    G.flip(self[label])

    def unpack_signed_DT(self, signed_dt):
        dt = []
        component = []
        flips = []
        for byte in bytearray(signed_dt):
            flips.append(bool(byte & 1 << 6))
            label = (1 + byte & 0x1f) << 1
            if byte & 1 << 5:
                label = -label
            component.append(label)
            if byte & 1 << 7:
                dt.append(tuple(component))
                component = []
        return dt, flips

    def convert_alpha(self, code):
        code = string_to_ints(code)
        num_crossings, components = code[:2]
        comp_lengths = code[2:2 + components]
        crossings = [x << 1 for x in code[2 + components:]]
        assert len(crossings) == num_crossings
        return partition_list(crossings, comp_lengths)

    def encode(self, header=True, alphabetical=True, flips=True):
        """
        Returns a string describing the DT code.  Options control
        whether to include the 'DT:' header, whether to use the
        numerical or alphabetical format, and whether to use the
        extended form, which adds flip information for each crossing.
        If flips is set to "auto", only include flips in large links
        (>26 crossings).

        >>> d = DTcodec([(-6,-8,-2,-4)])
        >>> A = d.encode()
        >>> A in ['DT:dadCDAB.0110', 'DT:dadCDAB.1001']
        True
        >>> N = d.encode(alphabetical=False)
        >>> N in ['DT:[(-6,-8,-2,-4)], [0,1,1,0]',
        ...  'DT:[(-6,-8,-2,-4)], [1,0,0,1]']
        True
        >>> d.encode(flips=False)
        'DT:dadCDAB'
        >>> d.encode(alphabetical=False, flips=False)
        'DT:[(-6,-8,-2,-4)]'
        """
        code = self.code
        result = 'DT:' if header else ''
        chunks = [len(component) for component in code]
        num_crossings = sum(chunks)
        is_large = (num_crossings > 26)
        if flips == 'auto':
            flips = is_large
        if alphabetical:
            if is_large:
                if flips:
                    result += encode_base64_like_DT_code(code, self.flips)
                else:
                    result += encode_base64_like_DT_code(code)
            else:
                prefix_ints = [num_crossings, len(code)]
                prefix_ints += chunks
                code_ints = [x for component in code for x in component]
                alphacode = ''.join(DT_alphabet[n >> 1]
                                    for n in code_ints)
                prefix = ''.join(DT_alphabet[n] for n in prefix_ints)
                if flips:
                    alphacode += '.' + ''.join(str(int(f)) for f in self.flips)
                result += (prefix + alphacode)
        else:
            result += str(code)
            if flips:
                result += ',  %s' % [int(f) for f in self.flips]
            result = result.replace(', ', ',')
        return result

    def embed(self, edge=None):
        """
        Try to flip crossings in the FatGraph until it becomes planar.
        """
        G = self.fat_graph
        # Add the marked_valence cache
        if edge is None:  # OK. No problem. Just pick one at random.
            for edge in G.edges:
                break
        # Find a circle and embed it.
        vertices, circle_edges = self.find_circle(edge)
        G.mark(circle_edges)
        # Add one arc, to get a theta graph.
        # The first arc needs to be bridge!
        first, last, arc_edges = G.bridge(circle_edges[:2])
        self.do_flips(first, arc_edges[0], last, arc_edges[-1])
        G.mark(arc_edges)
        # Keep adding arcs until the whole graph is embedded.
        while True:
            try:
                if not self.embed_arc():
                    break
            except EmbeddingError:
                flips = G.pop()
                for vertex in flips:
                    G.flip(vertex)
        # Clean up.
        self.fat_graph.clear_stack()

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

    def get_incomplete_vertex(self):
        """
        Return a vertex with some marked and some unmarked edges.
        If there are any, return a vertex with marked valence 3.
        """
        G = self.fat_graph
        vertices = [v for v in G.vertices if 0 < G.marked_valences[v] < 4]
        vertices.sort(key=lambda v: G.marked_valences[v])
        try:
            return vertices.pop()
        except IndexError:
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
        for k in range(1, 3):
            ccw_edge = G(v)[vslot + k]
            if ccw_edge.marked:
                break
        if not ccw_edge.marked:
            raise ValueError('Invalid marking')
        # Here are the slots and vertices in the two boundary curves:
        left_slots = G.left_slots(ccw_edge)
        right_slots = G.right_slots(ccw_edge)
        v_valence, w_valence = G.marked_valences[v], G.marked_valences[w]
        if (v, vslot) in left_slots:
            v_slot_side, v_other_side = left_slots, right_slots
        else:
            v_slot_side, v_other_side = right_slots, left_slots
        w_on_slot_side = w in [x[0] for x in v_slot_side]
        w_on_other_side = w in [x[0] for x in v_other_side]
        if not w_on_slot_side and not w_on_other_side:
            raise EmbeddingError('Embedding does not extend.')
        if (w, wslot) in v_slot_side:
            if v_valence == w_valence == 2:
                # This is an ambiguous situation.  We could either do
                # nothing or flip both vertices.  Push our state with
                # instructions to flip both vertices if we pop.
                G.push([v, w])
            return
        if w_valence != 2:
            G.flip(v)
            return
        elif v_valence != 2:
            G.flip(w)
            return
        # At this point we know that flips are needed *and* that both
        # vertices have marked valence 2.
        if w_on_slot_side and not w_on_other_side:
            G.flip(w)
            return
        if w_on_slot_side and w_on_other_side:
            # This is an ambiguous situation.  We know w does not have
            # a slot on the v slot side.  We are are about to flip v.
            # But we could leave v alone and flip w instead.  We
            # push our state with instructions to flip w if we
            # pop.
            G.push([w])
        G.flip(v)
        if (w, wslot) not in v_other_side:
            G.flip(w)

    def embed_arc(self):
        G = self.fat_graph
        v = self.get_incomplete_vertex()
        if v is None:
            return False
        if G.marked_valences[v] == 2:
            # This should work for any vertex if the diagram is prime.
            try:
                first, last, arc_edges = G.bridge(G.marked_arc(v))
            except ValueError:
                # If we can't find a bridge it means that the diagram
                # has a separating pair of edges.  We just take any
                # arc we can get, and see if we can find our way by
                # pushing and popping.
                arc_edges, last = G.unmarked_arc(v)
                first = v
            self.do_flips(first, arc_edges[0], last, arc_edges[-1])
        else:
            arc_edges, last_vertex = G.unmarked_arc(v)
            # Since v has 3 marked edges, we put it second here.
            self.do_flips(last_vertex, arc_edges[-1], v, arc_edges[0])
        G.mark(arc_edges)
        return True

    def signed_DT(self):
        """
        Return a byte sequence containing the signed DT code.

        >>> d = DTcodec([(-6,-8,-2,-4)])
        >>> d2 = DTcodec(d.signed_DT())
        >>> d2.code
        [(-6, -8, -2, -4)]
        """
        code_bytes = bytearray()
        it = iter(self.flips)
        for component in self.code:
            for label in component:
                byte = abs(label)
                byte = (byte >> 1) - 1
                if label < 0:
                    byte |= 1 << 5
                if next(it):
                    byte |= 1 << 6
                code_bytes.append(byte)
            code_bytes[-1] |= 1 << 7
        return bytes(code_bytes)

    def hex_signed_DT(self):
        """
        Return the hex encoding of the signed DT byte sequence.

        >>> d = DTcodec([(-6,-8,-2,-4)])
        >>> d2 = DTcodec(d.hex_signed_DT())
        >>> d2.code
        [(-6, -8, -2, -4)]
        """
        return '0x' + ''.join('%.2x' % b for b in bytearray(self.signed_DT()))

    def PD_code(self, KnotTheory=False):
        """
        Return a PD code for the projection described by this DT code,
        as a list of lists of 4 integers.  If KnotTheory is set to
        True, return a string that can be used as input to the Knot
        Theory package.

        >>> d = DTcodec([(-6,-8,-2,-4)], [0,1,1,0])
        >>> sorted(d.PD_code())
        [(2, 8, 3, 7), (4, 1, 5, 2), (6, 4, 7, 3), (8, 5, 1, 6)]
        """
        G = self.fat_graph
        PD = [G.PD_tuple(v) for v in G.vertices]
        if KnotTheory:
            return 'PD[%s]' % ', '.join('X%s' % repr(list(t)) for t in PD)
        return PD

    def link(self):
        G = self.fat_graph
        crossing_dict, slot_dict = {}, {}
        for v in G.vertices:
            c = Crossing(v[0])
            c.make_tail(0)
            if G.sign(v) == 1:
                c.make_tail(3)
            else:
                c.make_tail(1)
            c.orient()
            crossing_dict[v] = c
            if v.upper_pair() == (0, 2):
                slot_dict[v] = (3, 0, 1, 2)
            else:
                slot_dict[v] = (2, 3, 0, 1) if G.flipped(v) else (0, 1, 2, 3)
        for edge in G.edges:
            v0, v1 = edge
            s0, s1 = edge.slots
            a0, a1 = slot_dict[v0][s0], slot_dict[v1][s1]
            c0, c1 = crossing_dict[v0], crossing_dict[v1]
            c0[a0] = c1[a1]
        link = Link(list(crossing_dict.values()), check_planarity=False, build=False)
        assert link.all_crossings_oriented()
        component_starts = []
        i = 1
        for comp in self.code:
            c = self[i]
            if i == c[0]:
                e = slot_dict[c][2] if G.flipped(c) else slot_dict[c][0]
            else:
                e = slot_dict[c][1]
            ce = CrossingEntryPoint(crossing_dict[c], e)
            component_starts.append(ce)
            i += 2 * len(comp)
        link._build_components(component_starts)
        if not link.is_planar():
            raise ValueError('DT code does not seem to define a *planar* diagram')
        return link

    def KLPProjection(self):
        """
        Constructs a python simulation of a SnapPea KLPProjection
        (Kernel Link Projection) structure.  See DTFatGraph.KLP_dict
        and Jeff Weeks' SnapPea file link_projection.h for
        definitions.  Here the KLPCrossings are modeled by
        dictionaries.
        """
        G = self.fat_graph
        vertices = list(G.vertices)
        KLP_indices = {v: n for n, v in enumerate(vertices)}
        KLP_crossings = [G.KLP_dict(v, KLP_indices) for v in vertices]
        return len(G.vertices), 0, len(self.code), KLP_crossings

    def exterior(L):
        raise RuntimeError("SnapPy doesn't seem to be available.")
