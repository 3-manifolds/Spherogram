"""

Code for describing bands, adding them to links, and generating
reasonable collections of bands and the resulting banded links.

"""

import networkx as nx
from itertools import product, combinations
from ..links_base import CrossingStrand, CrossingEntryPoint, Strand, Crossing
from ..links import Link
from ..simplify import dual_graph_as_nx, Face
from .merge_links import are_isometric_as_links


def normalize_crossing_labels(link):
    for i, C in enumerate(link.crossings):
        C.label = i


def crossing_labels_are_normalized(link):
    return all(C.label == i for i, C in enumerate(link.crossings))


def dual_graph(link):
    """
    Returns the dual graph with the dual_edges relabeled by the
    stand_label of the corresponding original edge.  Eliminates the
    need to edit "simplify.py" as in earlier versions.
    """
    G = link.dual_graph()
    for edge in G.edges:
        edge.label = edge.interface[0].strand_label()
    return G


def print_crossings(link):
    """
    A crossing c has c.adjacent, which is a list of 4 tuples each
    (c,i) where c is a crossing label and i is a the index of the
    strand.  c.adjacent[j] is the crossing, strand tuple that crossing
    c strand j is glued to.
    """
    for c in link.crossings:
        print("Label = " +str(c.label))
        print("Sign = " +str(c.sign))
        print(c.adjacent)
        print(c.strand_labels)


def opposite(edge, face):
    if edge[0] == face:
        return edge[1]
    elif edge[1] == face:
        return edge[0]
    else:
        raise ValueError(f'Edge {edge} not on face {face}')


def edge_and_face_to_crossing_strands(edges, vertices):
    """
    Returns a dictionary D where::

      D[edge via strand_label, face as a vertex of dual graph]

    is the crossing strand corresponding to end of that edge that is
    in clockwise-first order with respect to the face; equivalently,
    this is the CrossingStrand whose orientation is compatible with
    the counterclockwise orientation of the boundary of the face.
    """

    edge_and_face_to_crossing = {}
    for face in vertices:
        for cs in face:
            edge_and_face_to_crossing[cs.strand_label(), face] = cs.opposite()
    return edge_and_face_to_crossing


def floyd_warshall(vertices, edges):
    # These are dual.vertices, as a list
    # and dual.edges, as a list and then sorted by label
    # dist[i, j] and nxt[i, j] where i and j are faces, as vertices
    # nxt[i, j] is an edge label
    dist = {}
    nxt = {}
    infty = 10*len(vertices)
    for face in vertices:
        for face2 in vertices:
            dist[face, face2] = infty
            nxt[face, face2] = None

    for face in vertices:
        dist[face, face] = 0

    for edge in edges:
        dist[edge[0], edge[1]] = 1
        dist[edge[1], edge[0]] = 1
        nxt[edge[0], edge[1]] = edge
        nxt[edge[1], edge[0]] = edge

    for k, i, j in product(vertices, repeat=3):
        sum_ik_kj = dist[i, k] + dist[k, j]
        if dist[i, j] > sum_ik_kj:
            dist[i, j] = sum_ik_kj
            nxt[i, j]  = nxt[i, k] # the edge you take next to get form i to j.

    for face1 in vertices:
        for face2 in vertices:
            if nxt[face1, face2] != None:
                nxt[face1, face2] = nxt[face1, face2].label
    return dist, nxt


def to_boolean_array(x, n):
    """
    Given an integer x return the corresponding boolean array of
    length n.  The convention is that i binary digit of x,
    corresponding to 2^i, gives the boolean value at index i.

    >>> [to_boolean_array(x, 2) for x in range(4)]
    [[False, False], [True, False], [False, True], [True, True]]
    >>> to_boolean_array(1 + 2**4, 6)
    [True, False, False, False, True, False]
    """
    return [ (x >> i) & 1 == 1 for i in range(n)]


def boolean_array_to_int(listlike):
    """
    >>> L = [1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1]
    >>> x = boolean_array_to_int(L); x
    1117
    >>> to_boolean_array(x, len(L)) == L
    True
    """
    x = 0
    for a in reversed(listlike):
        x = x << 1
        if a:
            x = x | 1
    return x


def crossing_strands_to_hex(crossing_strands):
    """
    For a knot with < 64 crossings, encode a sequence of n crossing
    strands as 2n hex digits, using a little endian strategy.  Thus,
    sequences with the same prefix map to strings with the same
    suffix.

    >>> cs = [(0,1), (4,3), (3,0), (7,1), (4, 0), (63, 3), (0, 0)]
    >>> a = crossing_strands_to_hex(cs); a
    '00ff101d0c1301'
    >>> b = crossing_strands_to_hex(cs[:4]); b
    '1d0c1301'
    >>> hex_to_crossing_strands(a) == cs
    True

    For knots with [64, 1,024) crossings, string is prefixed by 'Z' and
    uses 3n hex digits.

    >>> big_cs = [(111, 2), (231, 0), (1000, 3)]
    >>> encoded = crossing_strands_to_hex(big_cs); encoded
    'Zfa339c1be'
    >>> hex_to_crossing_strands(encoded) == big_cs
    True
    """
    n = len(crossing_strands)
    c = max(c for c, s in crossing_strands)
    if c >= 16384:
        raise ValueError('Too many crossings to encode')
    if c < 64:
        nibbles_per_cs = 2
        prefix = ''
    else:
        nibbles_per_cs = 3
        prefix = 'Z'
    x = 0
    for c, s in reversed(crossing_strands):
        x = (x << 4*nibbles_per_cs) + (c << 2) + s
    return prefix + '{:02x}'.format(x).rjust(nibbles_per_cs*n, '0')


def hex_to_crossing_strands(hex_string):
    if hex_string[0] == 'Z':
        nibbles_per_cs = 3
        hex_string = hex_string[1:]
        mask = 4095
    else:
        nibbles_per_cs = 2
        mask = 255

    n = len(hex_string) // nibbles_per_cs
    x = int(hex_string, 16)
    ans = []
    for i in range(n):
        y = x & mask
        s = y & 3
        c = y >> 2
        ans.append((c, s))
        x = x >> 4*nibbles_per_cs
    return ans


def compress_band_spec(along_top, arc_is_under, num_twist):
    cs = crossing_strands_to_hex(along_top)
    if not isinstance(arc_is_under, int):
        arc_is_under = boolean_array_to_int(arc_is_under)
    return cs + '_' + hex(arc_is_under)[2:] + '_' + repr(num_twist)


def uncompress_band_spec(spec):
    """
    >>> orig_spec = ([(0,1),(4,3),(4,0),(7,1)], [False,True], -2)
    >>> spec = compress_band_spec(*orig_spec)
    >>> uncompress_band_spec(spec) == orig_spec
    True
    """

    along_top, arc_is_under, num_twist = spec.split('_')
    along_top = hex_to_crossing_strands(along_top)
    arc_is_under = to_boolean_array(int(arc_is_under, base=16), len(along_top) - 2)
    num_twist = int(num_twist)
    return (along_top, arc_is_under, num_twist)


class Band:
    """
    A band in a link diagram is determined by the following data.

    * cs_along_top: A list specifying the CrossingStrands along the
      top of the band going from left to right.  The input should be a
      list of tuples (not actual CrossingStrands), where each crossing
      is indicated by its index (not an actual Crossing).

    * arc_is_under: a boolean array of length two less than
      cs_along_top. This determines how the band weaves over or under
      the successive stands of the original link. The "arc" in
      "arc_is_under" refers the strand of the link, not the band, so
      this parameter could also have been called "band_is_over".

    * num_twists: the number of twists to put in the band.

    Both cs_along_top and arc_is_under can be encoded as hex strings
    for brevity, and arc_is_under can also be an integer encoding the
    corresponding bitstring.

    >>> B0 = Band([(2,1), (7,3), (7,0), (4,1)], [False, True], -2); B0
    Band([(2, 1), (7, 3), (7, 0), (4, 1)], [False, True], -2)
    >>> B1 = Band([(2,1), (7,3), (7,0), (4,1)], 2, -2)
    >>> B2 = Band('111c1f09', '2', -2)
    >>> B3 = Band('111c1f09_2_-2')
    >>> B0 == B1 == B2 == B3
    True
    >>> B0.compressed_spec()
    '111c1f09_2_-2'
    >>> B4 = Band([(2,1), (7,0), (7,3)], [True], 1)
    >>> B4 < B0 and B4 <= B0 and B0 > B4 and B0 >= B4 and B0 != B4
    True
    """

    def __init__(self, cs_along_top, arc_is_under=None, num_twist=None):
        if arc_is_under is None and num_twist is None:
            cs_along_top, arc_is_under, num_twist = cs_along_top.split('_')
            num_twist = int(num_twist)

        if isinstance(cs_along_top, str):
            cs_along_top = hex_to_crossing_strands(cs_along_top)

        if isinstance(arc_is_under, str):
            arc_is_under = int(arc_is_under, 16)
        if isinstance(arc_is_under, int):
            arc_is_under = to_boolean_array(arc_is_under, len(cs_along_top) - 2)

        assert len(cs_along_top) >= 2
        assert len(arc_is_under) == len(cs_along_top) - 2

        self.cs_along_top =  cs_along_top
        self.arc_is_under = arc_is_under
        self.num_twist = num_twist
        self._spec = (cs_along_top, arc_is_under, num_twist)

    def compressed_spec(self):
        return compress_band_spec(*self._spec)

    def __len__(self):
        return len(self.cs_along_top)

    def __eq__(self, other):
        if isinstance(other, Band):
            return self._spec == other._spec

    def __neq__(self, other):
        if isinstance(other, Band):
            return self._spec != other._spec

    def __lt__(self, other):
        if isinstance(other, Band):
            return (len(self), self._spec) < (len(other), other._spec)

    def __le__(self, other):
        if isinstance(other, Band):
            return (len(self), self._spec) <= (len(other), other._spec)

    def __gt__(self, other):
        if isinstance(other, Band):
            return (len(self), self._spec) > (len(other), other._spec)

    def __ge__(self, other):
        if isinstance(other, Band):
            return (len(self), self._spec) >= (len(other), other._spec)

    def __repr__(self):
        return f'Band{self._spec}'

    def is_nonminimal(self, link=None):
        """
        Do various tests to see if self is obviously equivalent to a
        band other with other < self.
        """
        cs_along_top = self.cs_along_top
        band_is_over = self.arc_is_under
        if len(cs_along_top) > 2:
            # Check if can slide an attaching point of the band
            # through crossing to get a strictly shorter band.

            def can_slide(cs_attach, cs_next, next_over):
                (A, a), (B, b) = cs_attach, cs_next
                if A != B:
                    return False
                attach_over = (a % 2 == 1)
                return ((attach_over and next_over) or
                        (not attach_over and not next_over))

            if (can_slide(cs_along_top[0], cs_along_top[1], band_is_over[0]) or
                can_slide(cs_along_top[-1], cs_along_top[-2], band_is_over[-1])):
                return True

            if link is not None:
                crossings = link.crossings
                for i, j in [(0, 1), (-1, -2)]:
                    (A, a), (B, b) = cs_along_top[i], cs_along_top[j]
                    cs_a = CrossingStrand(crossings[A], a).opposite()
                    cs_b = CrossingStrand(crossings[B], b).opposite()
                    if can_slide(cs_a, cs_b, band_is_over[i]):
                        return True

        return False


def add_one_band(link, band):
    """
    Creates a new link by banding its input along the specified
    band. See the documentation Band for our conventions.

    >>> L = Link('K8n1')
    >>> L1 = add_one_band(L, Band([(0,1),(4,3),(4,0),(7,1)], [False,True], -2))
    >>> are_isometric_as_links(L1, Link('L9n9'))    #doctest: +SNAPPY
    True
    >>> L2 = add_one_band(L, Band([(1,1),(1,2),(4,0),(3,3)], [False,True], 1))
    >>> are_isometric_as_links(L2, Link('L9a16'))   #doctest: +SNAPPY
    True

    One can also give the band by its compressed spec.

    >>> L3 = add_one_band(L, '111c1f09_2_-2')
    >>> are_isometric_as_links(L3, Link('L10n36'))  #doctest: +SNAPPY
    True

    Double-checking that arc_is_under being True means the *strand* of
    the knot is under the *band*.  Equivalently, can think of
    arc_is_under as band_is_over.

    >>> L = Link('K8n1')
    >>> band_below_strand = Band([(0, 1), (0, 2), (7, 0)], [False], 0)
    >>> E = add_one_band(L, band_below_strand)
    >>> are_isometric_as_links(E, Link('L8n2'))    #doctest: +SNAPPY
    True
    >>> band_above_strand = Band([(0, 1), (0, 2), (7, 0)], [True], 0)
    >>> L_new = add_one_band(L, band_above_strand)
    >>> success = L_new.simplify('global'); len(L_new.crossings)
    2
    """
    if isinstance(band, str):
        band = Band(band)
    elif not isinstance(band, Band):
        band = Band(*band)

    cs_along_top = band.cs_along_top
    arc_is_under = band.arc_is_under
    num_twist = band.num_twist
    assert crossing_labels_are_normalized(link)
    L = link.copy()

    positive_twist = num_twist >= 0
    num_twist = abs(num_twist)
    along_top = [CrossingStrand(L.crossings[c], s) for c, s in cs_along_top]

    X = along_top[0]
    Z = along_top[-1]
    Y, W = X.opposite(), Z.opposite()

    A = along_top[1:-1]
    D = [a.opposite() for a in A]
    num_arcs = len(A)

    crossing_label = len(L.crossings)
    B = [Crossing() for i in range(num_arcs)]
    C = [Crossing() for i in range(num_arcs)]
    for i in range(num_arcs):
        B[i].label = crossing_label + i
        C[i].label = crossing_label + i + num_arcs

    # Here we deal with vertical strands between A, B, and C type
    # crossings recall that A[i] and D[i] are tuples (crossing,
    # index), recording the loose strands from above and below whereas
    # B[i] and C[i] are just crossings.

    for i in range(num_arcs):
        # A[i] is a crossing, index tuple, with A[i][0] the crossing,
        # and A[i][1] the index
        if arc_is_under[i]:
            A[i][0].adjacent[A[i][1]] = (B[i], 2)
            D[i][0].adjacent[D[i][1]] = (C[i], 0)
            B[i].adjacent[2] = (A[i][0], A[i][1])
            B[i].adjacent[0] = (C[i], 2)
            C[i].adjacent[0] = (D[i][0], D[i][1])
            C[i].adjacent[2] = (B[i], 0)
        else:
            A[i][0].adjacent[A[i][1]] = (B[i], 3)
            D[i][0].adjacent[D[i][1]] = (C[i], 1)
            B[i].adjacent[3] = (A[i][0], A[i][1])
            B[i].adjacent[1] = (C[i], 3)
            C[i].adjacent[1] = (D[i][0], D[i][1])
            C[i].adjacent[3] = (B[i], 1)

    # Now we add the E crossings which are for the twists, depending on
    # whether they are right handed or left handed twists.

    crossing_label += 2*num_arcs

    E = [Crossing() for i in range(num_twist)]
    for i in range(num_twist):
        E[i].label = crossing_label+i


    # now we add the horizonal arcs between the E[i] crossings

    if positive_twist:
        for i in range(num_twist - 1):
            E[i].adjacent[1] = (E[i+1], 0)
            E[i].adjacent[2] = (E[i+1], 3)
    else:
        for i in range(num_twist - 1):
            E[i].adjacent[1] = (E[i+1], 2)
            E[i].adjacent[0] = (E[i+1], 3)

    # Connecting the left strands X and Y to the E crossings

    upper = X
    lower = Y

    if positive_twist:
        for i in range(num_twist):
            upper[0].adjacent[upper[1]] = (E[i],3)
            lower[0].adjacent[lower[1]] = (E[i],0)
            E[i].adjacent[3] = (upper[0], upper[1])
            E[i].adjacent[0] = (lower[0], lower[1])
            upper = (E[i],2)
            lower = (E[i],1)
    else:
        for i in range(num_twist):
            upper[0].adjacent[upper[1]] = (E[i],2)
            lower[0].adjacent[lower[1]] = (E[i],3)
            E[i].adjacent[2] = (upper[0], upper[1])
            E[i].adjacent[3] = (lower[0], lower[1])
            upper = (E[i],1)
            lower = (E[i],0)


    # Now adding in the horizontal gluings in the B and C region

    for i in range(num_arcs):
        if arc_is_under[i]:
            upper[0].adjacent[upper[1]] = (B[i],3)
            lower[0].adjacent[lower[1]] = (C[i],3)
            B[i].adjacent[3] = (upper[0], upper[1])
            C[i].adjacent[3] = (lower[0], lower[1])

            upper = (B[i],1)
            lower = (C[i],1)
        else:
            upper[0].adjacent[upper[1]] = (B[i],0)
            lower[0].adjacent[lower[1]] = (C[i],0)
            B[i].adjacent[0] = (upper[0], upper[1])
            C[i].adjacent[0] = (lower[0], lower[1])

            upper = (B[i],2)
            lower = (C[i],2)

    upper[0].adjacent[upper[1]] = Z
    lower[0].adjacent[lower[1]] = W
    Z[0].adjacent[Z[1]] = (upper[0], upper[1])
    W[0].adjacent[W[1]] = (lower[0], lower[1])

    L.crossings = L.crossings + B + C + E
    L._rebuild()
    assert L.is_planar()
    if L.name is not None:
        L.name += '+band'
    return L


def min_len_bands(link, max_twists=2, max_band_len=None):
    """
    Produces bands on the input link where the band must start and end
    on the same link component with its ends oriented so that the
    result of banding has one more component than the original. As the
    name suggests, it only considers minimal length paths in the dual
    graph to the projection.

    Does not modify the input link.

    >>> L = Link('K6a3')
    >>> len(min_len_bands(L))
    138
    >>> L = Link('K8n1')
    >>> some_bands = min_len_bands(L)
    >>> len(some_bands)
    315
    >>> some_bands[0]
    Band([(0, 0), (1, 1)], [], -2)
    >>> some_bands[-1].compressed_spec()
    '15030212_2_2'
    >>> L = Link('L5a1')
    >>> len(min_len_bands(L))
    54

    Can also specify a maximum length for the band, as measured in the
    length of the cs_along_top (so every band has length at least two).

    >>> L = Link('L8n1')
    >>> len(min_len_bands(L, max_twists=1, max_band_len=2))
    28

    """
    # dual.vertices is a set of Face objects, which are list of
    # CrossingStrand objects dual.edges is a set of Edges, each edge e
    # has e[0] and e[1] which are faces
    L = link
    if not crossing_labels_are_normalized(L):
        raise ValueError('Link needs normalized crossing labels')
    normalize_crossing_labels(L)
    dual = dual_graph(L)
    vertices = sorted(dual.vertices, key=lambda x:x.label)
    edges = sorted(dual.edges, key=lambda x:x.label)

    # a dictionary that given an edge (via strand_label) and a face
    # (as vertex of dual graph, gives you the two crossing strands as
    # ((c1, i1), (c2,i2)) going counter clockwise.
    edge_and_face_to_cs = edge_and_face_to_crossing_strands(edges, vertices)

    # Compute the shortest path between each pair of vertices in the dual.
    dist, nxt = floyd_warshall(vertices, edges)

    # Now interate over each pair of arcs in the diagram that we can
    # band together.

    ans = []
    for component in L.link_components:
        strands = [cs.strand_label() for cs in component]
        for arc1, arc2 in combinations(strands, 2):
            # First, determine the shortest path from arc1 to arc2,
            # there are two faces adjacent to each, so four
            # possibilities.
            faces_poss = list(product(edges[arc1], edges[arc2]))
            num_arcs = min(dist[x] for x in faces_poss)
            if max_band_len != None and num_arcs > max_band_len - 2:
                continue
            # There can be more than one minimum length path. Based on
            # K13a2084, it pays to try them all.
            for mini, minj in [x for x in faces_poss if dist[x] == num_arcs]:
                X = edge_and_face_to_cs[arc1, mini].opposite()
                Z = edge_and_face_to_cs[arc2, minj]
                A = []
                i = mini
                while dist[i, minj] > 0:
                    next_edge = edges[nxt[i, minj]]
                    A.append(edge_and_face_to_cs[next_edge.label, i])
                    i = opposite(next_edge, i)

                along_top = [X] + A + [Z]
                along_top = [(cs.crossing.label, cs.strand_index) for cs in along_top]

                X_orient = (X == X.oriented())
                Z_orient = (Z == Z.oriented())
                correct_parity = 1 if X_orient == Z_orient else 0
                good_twists = [twist for twist in range(-max_twists, max_twists + 1)
                               if twist % 2 == correct_parity]

                for num_twist in good_twists:
                    # arc_is_under is used here to iterate over bitstrings
                    for arc_is_under in range(2**num_arcs):
                        band = Band(along_top, arc_is_under, num_twist)
                        if not band.is_nonminimal(L):
                            ans.append(band)
    ans.sort()
    return ans


def simple_bands(link, max_twists=2, max_band_len=None):
    """
    Produces bands on the input link where the band must start and end
    on the same link component with its ends oriented so that the
    result of banding has one more component than the original.

    The "simple" in "simple_bands" refers to the fact that they come
    from simple (i.e. embedded) paths in the dual graph, even though
    they can be quite long.

    Does not modify the input link.

    >>> L = Link('K6a3')
    >>> len(simple_bands(L, max_band_len=3))
    138
    >>> L = Link('K8n1')
    >>> some_bands = simple_bands(L, max_band_len=2)
    >>> len(some_bands)
    106
    >>> L = Link('L5a1')
    >>> some_bands = simple_bands(L)
    >>> len(some_bands)
    393
    >>> some_bands[-1].compressed_spec()
    '0a0504111013_e_1'

    Can also specify a maximum length for the band, as measured in the
    length of the cs_along_top (so every band has length at least two).

    >>> L = Link('L8n1')
    >>> len(simple_bands(L, max_twists=1, max_band_len=2))
    28

    """
    # dual.vertices is a set of Face objects, which are list of
    # CrossingStrand objects dual.edges is a set of Edges, each edge e
    # has e[0] and e[1] which are faces
    L = link
    if not crossing_labels_are_normalized(L):
        raise ValueError('Link needs normalized crossing labels')
    normalize_crossing_labels(L)
    # nodes are integers, recording the face number.
    dual = dual_graph_as_nx(L, graph_class=nx.MultiGraph)
    cutoff = max_band_len - 2 if max_band_len is not None else None

    # Mapping from CrossingStrands to faces, namely the face that the
    # CrossingStrand orients clockwise.

    cs_to_face = dict()
    for F0, F1, attrs in dual.edges(data=True):
        interface = attrs['interface']
        cs0, cs1 = interface[F0], interface[F1]
        cs_to_face[cs0] = F0
        cs_to_face[cs1] = F1

    # Now iterate over each pair of arcs in the diagram that we can
    # band together.
    ans = []
    for component in L.link_components:
        for cs0, cs1 in combinations(component, 2):
            # There are two faces adjacent to each, so four
            # possibilities.

            faces0 = [cs_to_face[cs0], cs_to_face[cs0.opposite()]]
            faces1 = [cs_to_face[cs1], cs_to_face[cs1.opposite()]]

            for F0, F1 in product(faces0, faces1):
                if F0 == F1:
                    paths = [[F0]]
                else:
                    paths = nx.all_simple_paths(dual, F0, F1, cutoff=cutoff)
                for path in paths:
                    # First augment the path to include the dual edges
                    # corresponding to cs0 and cs1.
                    P0 = [F for F in faces0 if F != F0][0]
                    P1 = [F for F in faces1 if F != F1][0]
                    path = [P0] + path + [P1]
                    if not nx.is_simple_path(dual, path):
                        continue
                    # Extract the crossing strands along the top.
                    # Code looks complicated because this is a
                    # MultiGraph.
                    along_top = []
                    for A, B in zip(path, path[1:]):
                        for edge, info in dual[A][B].items():
                            if B in info['interface']:
                                along_top.append(info['interface'][B])
                                break
                    assert len(along_top) >= 2
                    X, Z = along_top[0], along_top[-1]
                    along_top = [(cs.crossing.label, cs.strand_index) for cs in along_top]

                    X_orient = (X == X.oriented())
                    Z_orient = (Z == Z.oriented())
                    correct_parity = 1 if X_orient == Z_orient else 0
                    good_twists = [twist for twist in range(-max_twists, max_twists + 1)
                                   if twist % 2 == correct_parity]
                    for num_twist in good_twists:
                        # arc_is_under is used here to iterate over bitstrings
                        num_arcs = len(along_top) - 2
                        for arc_is_under in range(2**num_arcs):
                            band = Band(along_top, arc_is_under, num_twist)
                            if not band.is_nonminimal(L):
                                ans.append(band)
    ans.sort()
    return ans



def banded_links(link, max_twists=2, max_band_len=None, paths='shortest'):
    """
    Produces bands on the input link where the band must start and end
    on the same link component with its ends oriented so that the
    result of banding has one more component than the original.

    Does not modify the input link. Considers only the bands produced
    by ``min_len_bands`` or ``simple_bands`` according to the
    ``paths`` argument.


    >>> L = Link('K6a3')
    >>> len(list(banded_links(L)))
    138
    >>> L = Link('K8n1')
    >>> len(list(banded_links(L)))
    315
    >>> L = Link('L5a1')
    >>> len(list(banded_links(L, paths='simple')))
    393
    """
    if not crossing_labels_are_normalized(link):
        link = link.copy()
        normalize_crossing_labels(link)

    if paths == 'shortest':
        paths = min_len_bands(link, max_twists, max_band_len)
    elif paths == 'simple':
        paths = simple_bands(link, max_twists, max_band_len)
    else:
        raise ValueError("Path type can be either 'shortest' or 'simple'")

    for band in paths:
        B = add_one_band(link, band)
        assert len(B.link_components) == len(link.link_components) + 1
        normalize_crossing_labels(B)
        yield B, band.compressed_spec()
