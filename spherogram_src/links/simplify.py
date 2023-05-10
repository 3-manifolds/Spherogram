"""
Simplifying link diagrams.

Important notes:

* The link diagram is modified in place.  All the relevant parts of the
  data structure are updated at each step.

* Unknot components which are also unlinked may be silently discarded.
"""

from .links import Link, Strand, Crossing, CrossingStrand
from .ordered_set import OrderedSet
from .. import graphs
import random
import networkx as nx
import collections


def remove_crossings(link, eliminate):
    """
    Deletes the given crossings. Assumes that they have already been
    disconnected from the rest of the link, so this just updates
    link.crossings and link.link_components.
    """

    if len(eliminate):
        for C in eliminate:
            link.crossings.remove(C)
        new_components = []
        for component in link.link_components:
            for C in eliminate:
                for cep in C.entry_points():
                    try:
                        component.remove(cep)
                    except ValueError:
                        pass
            if len(component):
                new_components.append(component)
        components_removed = len(link.link_components) - len(new_components)
        link.unlinked_unknot_components += components_removed
        link.link_components = new_components


def reidemeister_I(link, C):
    """
    Does a type-1 simplification on the given crossing C if possible.

    Returns the pair: {crossings eliminated}, {crossings changed}
    """
    elim, changed = set(), set()
    for i in range(4):
        if C.adjacent[i] == (C, (i + 1) % 4):
            (A, a), (B, b) = C.adjacent[i + 2], C.adjacent[i + 3]
            elim = set([C])
            if C != A:
                A[a] = B[b]
                changed = set([A, B])

    remove_crossings(link, elim)
    return elim, changed


def reidemeister_I_and_II(link, A):
    """
    Does a type-1 or type-2 simplification at the given crossing A if
    possible.

    Returns the pair: {crossings eliminated}, {crossings changed}
    """
    eliminated, changed = reidemeister_I(link, A)
    if not eliminated:
        for a in range(4):
            (B, b), (C, c) = A.adjacent[a], A.adjacent[a + 1]
            if B == C and (b - 1) % 4 == c and (a + b) % 2 == 0:
                eliminated, changed = reidemeister_I(link, B)
                if eliminated:
                    break
                else:
                    W, w = A.adjacent[a + 2]
                    X, x = A.adjacent[a + 3]
                    Y, y = B.adjacent[b + 1]
                    Z, z = B.adjacent[b + 2]
                    eliminated = set([A, B])
                    if W != B:
                        W[w] = Z[z]
                        changed.update(set([W, Z]))
                    if X != B:
                        X[x] = Y[y]
                        changed.update(set([X, Y]))
                    remove_crossings(link, eliminated)
                    break

    return eliminated, changed


def basic_simplify(link, build_components=True, to_visit=None,
                   force_build_components=False):
    """
    Do Reidemeister I and II moves until none are possible.
    """
    if to_visit is None:
        to_visit = set(link.crossings)
    eliminated = set()
    while to_visit:
        crossing = to_visit.pop()
        elim, changed = reidemeister_I_and_II(link, crossing)
        assert not elim.intersection(changed)
        eliminated.update(elim)
        to_visit.difference_update(elim)
        to_visit.update(changed)

    success = len(eliminated) > 0

    # Redo the strand labels (used for DT codes)
    if (success and build_components) or force_build_components:
        component_starts = []
        for component in link.link_components:
            assert len(component) > 0
            if len(component) > 1:
                a, b = component[:2]
            else:
                a = component[0]
                b = a.next()
            if a.strand_label() % 2 == 0:
                component_starts.append(a)
            else:
                component_starts.append(b)
        link._build_components(component_starts)
    return success


def possible_type_III_moves(link):
    """
    Returns all triples of crossings where a type III move is possible.

    In this example, one type III move is forbidden since a crossing
    repeats twice.

    >>> L = Link([(2,1,3,2),(3,8,4,1),(4,6,5,5),(6,8,7,7)])
    >>> len(possible_type_III_moves(L))
    1
    """
    ans = []
    for face in link.faces():
        if len(face) == 3:
            if sum(ce.strand_index % 2 for ce in face) in [1, 2]:
                while(face[1][1] % 2 != 0 or face[2][1] % 2 != 1):    # renumber face_list
                    face = [face[1], face[2], face[0]]
                if len(set(e.crossing for e in face)) == 3:  # No repeated crossings
                    ans.append(face)
    return ans


def insert_strand(X, x):
    Y, y = X.adjacent[x]
    S = Strand()
    S[0], S[1] = X[x], Y[y]
    return S


def reidemeister_III(link, triple):
    """
    Performs the given type III move.  Modifies the given link but doesn't
    update its lists of link components.
    """
    A, B, C = [t.crossing for t in triple]
    a, b, c = [t.strand_index for t in triple]
    # We insert Strands around the border of the triple to make the code more
    # transparent and eliminate some special cases.
    old_border = [(C, c - 1), (C, c - 2), (A, a - 1),
                  (A, a - 2), (B, b - 1), (B, b - 2)]
    border_strands = [insert_strand(*P) for P in old_border]
    new_boarder = [(A, a), (B, b + 1), (B, b), (C, c + 1), (C, c), (A, a + 1)]
    for i, (X, x) in enumerate(new_boarder):
        X[x] = border_strands[i][0]
    A[a - 1], B[b - 1], C[c - 1] = B[b + 2], C[c + 2], A[a + 2]
    [S.fuse() for S in border_strands]


def simplify_via_level_type_III(link, max_consecutive_failures=100):
    """
    Applies a series of type III moves to the link, simplifying it via type
    I and II moves whenever possible.
    """
    failures, success = 0, False
    if basic_simplify(link):
        success = True
    while failures < max_consecutive_failures:
        poss_moves = possible_type_III_moves(link)
        if len(poss_moves) == 0:
            break
        reidemeister_III(link, random.choice(poss_moves))
        if basic_simplify(link):
            failures = 0
            success = True
        else:
            failures += 1

    link._build_components()
    return success


def common_element(X, Y):
    return next(iter(set(X) & set(Y)))


class Face(tuple):
    """
    A complementary region of the link diagram.
    """
    def __new__(cls, edges, label=None):
        ans = tuple.__new__(cls, edges)
        ans.label = label
        return ans

    def __repr__(self):
        return "<F%d>" % self.label


class DualGraphOfFaces(graphs.Graph):
    """
    The dual graph to a link diagram D, whose vertices correspond to
    complementary regions (faces) of D and whose edges are dual to the
    edges of D.
    """
    def __init__(self, link):
        graphs.Graph.__init__(self)
        faces = [Face(face, i) for i, face in enumerate(link.faces())]
        self.edge_to_face = to_face = {}
        for face in faces:
            for edge in face:
                to_face[edge] = face

        for edge, face in to_face.items():
            neighbor = to_face[edge.opposite()]
            if face.label < neighbor.label:
                dual_edge = self.add_edge(face, neighbor)
                dual_edge.interface = (edge, edge.opposite())
                dual_edge.label = len(self.edges) - 1

        # assert self.is_planar()

    def two_cycles(self):
        """
        Find all two cycles and yield them as a pair of CrossingStrands which
        are dual to the edges in the cycle.

        The crossing strands are
        oriented consistently with respect to one of the faces which a
        vertex for the cycle.
        """
        for face0 in self.vertices:
            for dual_edge0 in self.incident(face0):
                face1 = dual_edge0(face0)
                if face0.label < face1.label:
                    for dual_edge1 in self.incident(face1):
                        if dual_edge0.label < dual_edge1.label and dual_edge1(face1) == face0:
                            yield (common_element(face0, dual_edge0.interface),
                                   common_element(face0, dual_edge1.interface))


def dual_graph_as_nx(link):
    corners = OrderedSet([CrossingStrand(c, i)
                          for c in link.crossings for i in range(4)])
    faces = []
    to_face_index = dict()
    while len(corners):
        count = len(faces)
        first_cs = corners.pop()
        to_face_index[first_cs] = count
        face = [first_cs]
        while True:
            next = face[-1].next_corner()
            if next == face[0]:
                faces.append(Face(face, count))
                break
            else:
                to_face_index[next] = count
                corners.remove(next)
                face.append(next)

    G = nx.Graph()
    to_face = {edge: faces[f] for edge, f in to_face_index.items()}

    for edge, face in to_face.items():
        opp_edge = edge.opposite()
        neighbor = to_face[opp_edge]
        if face.label < neighbor.label:
            G.add_edge(face.label, neighbor.label,
                       interface={face.label: edge, neighbor.label: opp_edge})

    G.edge_to_face = to_face_index
    return G


def deconnect_sum(link):
    """
    Warning: Destroys the original link.
    """
    for cs0, cs1 in DualGraphOfFaces(link).two_cycles():
        A, a = cs0.opposite()
        B, b = cs0
        C, c = cs1.opposite()
        D, d = cs1
        A[a] = D[d]
        B[b] = C[c]
    link._build_components()
    return link.split_link_diagram(destroy_original=True)


def dual_edges(overstrand, graph):
    """
    Find the set of crossings and edges of the dual graph encountered
    by moving along the link starting at startcep for length crossings.
    Also returns the next crossing entry point immediately after.
    """

    edges_crossed = []
    for cep in overstrand:
        f1 = graph.edge_to_face[cep]
        f2 = graph.edge_to_face[cep.opposite()]
        edges_crossed.append((f1, f2))

    # want one more edge
    endpoint = overstrand[-1].next()
    final_f1 = graph.edge_to_face[endpoint]
    final_f2 = graph.edge_to_face[endpoint.opposite()]
    edges_crossed.append((final_f1, final_f2))

    return edges_crossed


def extend_strand_forward(kind, strand, end_cep):
    """
    Extend the strand by adding on end_cep and what comes after it
    until you hit a crossing which is not of the given kind
    (over/under).
    """
    cep = end_cep.next()
    strand.append(end_cep)
    start_cep = strand[0].previous()
    while getattr(cep, 'is_' + kind + '_crossing')():
        if cep.next() == start_cep:  # prevents extending too far
            break
        strand.append(cep)
        cep = cep.next()
        if cep == strand[0]:
            break


def extend_strand_backward(kind, strand, start_cep):
    """
    Extend the strand by adding on end_cep and what comes before it
    until you hit a crossing which is not of the given kind
    (over/under).
    """
    cep = start_cep.previous()
    strand.insert(0, start_cep)
    end_cep = strand[-1].next()
    while getattr(cep, 'is_' + kind + '_crossing')():
        if cep.previous() == end_cep:  # prevents extending too far
            break
        strand.insert(0, cep)
        cep = cep.previous()
        if cep == strand[-1]:
            break


def pickup_strand(link, dual_graph, kind, strand):
    """
    Simplify the given (over/under)crossing strand by erasing it from
    the diagram and then finding a path that minimizes the number of
    edges it has to cross over to connect the same endpoints. Returns
    number of crossings removed.
    """
    init_link_cross_count = len(link.crossings)
    G = dual_graph
    startcep = strand[0].previous()

    if startcep == strand[-1]:
        # Totally overcrossing loop, must be totally unlinked and
        # unknotted
        remove_strand(link, strand)
        return len(strand)
    if startcep == strand[-1].next() and startcep.other() in strand:
        # We have a figure-8 curve with a single crossing in front
        # of the rest of the components.
        remove_strand(link, [startcep] + strand)
        return len(strand)
    crossing_set = set(cep.crossing for cep in strand)
    endpoint = strand[-1].next()

    if endpoint.crossing in crossing_set:
        # Strand crosses itself underneath
        extend_strand_forward(kind, strand, endpoint)
        return pickup_strand(link, G, kind, strand)
    if startcep.crossing in crossing_set:
        # Strand crosses itself over
        extend_strand_backward(kind, strand, startcep)
        return pickup_strand(link, G, kind, strand)

    edges_crossed = dual_edges(strand, G)
    source = edges_crossed[0][0]
    dest = edges_crossed[-1][0]

    nx.set_edge_attributes(G, 1, 'weight')
    for (f0, f1) in edges_crossed:
        G[f0][f1]['weight'] = 0

    path = nx.shortest_path(G, source, dest, weight='weight')
    new_len = sum(G[f0][f1]['weight'] for f0, f1 in zip(path, path[1:]))
    crossingsremoved = len(crossing_set) - new_len
    if crossingsremoved == 0:
        return 0

    # creating a new list of crossings from which to rebuild the link,
    # remove old overcross
    newcrossings = []

    removed = remove_strand(link, strand)
    loose_end = startcep.rotate(2)

    # find new sequence of overcrossings to create
    for i in range(len(path) - 1):
        label = 'new%d' % i
        f1, f2 = path[i:i + 2]
        edge = G[f1][f2]
        if edge['weight'] > 0:
            cep_to_cross = G[f1][f2]['interface'][f1]
            new_crossing, loose_end = cross(link, cep_to_cross, kind, loose_end, label)
            newcrossings.append(new_crossing)

    ec, ecep = endpoint.crossing, endpoint.strand_index
    ec[ecep] = loose_end

    link.crossings.extend(newcrossings)
    active = set()
    for C in removed:
        for i in range(4):
            D = C.adjacent[i][0]
            if D not in removed:
                active.add(D)

    active.update(newcrossings)
    basic_simplify(link, force_build_components=True, to_visit=active)

    final_cross_removed = init_link_cross_count - len(link.crossings)
    assert final_cross_removed >= crossingsremoved
    return final_cross_removed


def strand_pickup(link, kind):
    """
    Simplifies link by optimizing the path of the longest sequence of overcrossings.
    Returns a new link and the number of crossings removed.
    """
    G = None
    strands = over_or_under_strands(link, kind)
    for strand in randomize_within_lengths(strands):
        if len(strand) == 1:
            continue
        if G is None:
            G = dual_graph_as_nx(link)
        crossings_removed = pickup_strand(link, G, kind, strand)
        if crossings_removed != 0:
            return crossings_removed
    return 0


def remove_strand(link, strand):
    """
    Delete an overstrand or understrand from a link.  If the strand is
    a loop, it doesn't leave any loose strands and removes the
    loop. Otherwise, there will be two strands left in the link not
    attached to anything.  This function assumes that the start and
    end of the strand are not places where strands crosses itself.
    """
    # only add bridge strands for the places where the strand doesn't cross itself
    crossings_seen = [s.crossing for s in strand]
    crossing_set = set()
    for c in crossings_seen:
        if c in crossing_set:
            crossing_set.remove(c)
        else:
            crossing_set.add(c)

    start_cep = strand[0].previous()
    end_cep = strand[-1].next()
    bridge_strands = {c: Strand('strand' + str(c.label)) for c in crossing_set}
    for cep in strand:
        c = cep.crossing
        if c not in crossing_set:
            continue
        right_cs = cep.rotate(1).opposite()
        left_cs = cep.rotate(3).opposite()
        if right_cs.crossing in crossing_set:
            signs_equal = (c.sign == right_cs.crossing.sign)
            bridge_strands[c][0] = bridge_strands[right_cs.crossing][signs_equal]
        else:
            bridge_strands[c][0] = right_cs.crossing[right_cs.strand_index]
        if left_cs.crossing in crossing_set:
            signs_equal = (c.sign == left_cs.crossing.sign)
            bridge_strands[c][1] = bridge_strands[left_cs.crossing][1 - signs_equal]
        else:
            bridge_strands[c][1] = left_cs.crossing[left_cs.strand_index]
    remove_crossings(link, set(crossings_seen))

    for s in bridge_strands.values():
        s.fuse()

    return set(crossings_seen)


def cross(link, cep_to_cross, kind, loose_end, label):
    """
    Create a new (over/under) crossing along the edge defined by
    cep_to_cross and its opposite, and attach one side to a given
    position loose_end.
    """

    ic = cep_to_cross
    ico = ic.opposite()
    while ic.crossing not in link.crossings:
        ic = ic.next()
    while ico.crossing not in link.crossings:
        ico = ico.next()
    left_to_right = not ic.crossing.is_incoming(ic.strand_index)

    if kind == 'over':
        if left_to_right:
            a, b, c, d = 2, 3, 0, 1
        else:
            a, b, c, d = 0, 1, 2, 3
    else:
        assert kind == 'under'
        a, b, c, d = 3, 0, 1, 2

    new_crossing = Crossing(label)
    new_crossing[b] = loose_end
    new_crossing[c] = ic
    new_crossing[a] = ico
    new_crossing.make_tail(b)
    if left_to_right:
        new_crossing.make_tail(c)
    else:
        new_crossing.make_tail(a)
    new_crossing.orient()
    return new_crossing, new_crossing.crossing_strands()[d]


def over_or_under_strands(link, kind):
    """
    Returns a list of the sequences of (over/under) crossings (which
    are lists of CrossingEntryPoints), sorted in descending order
    of length.
    """

    def criteria(cep):
        return getattr(cep, 'is_' + kind + '_crossing')()

    ceps = OrderedSet(
        [cep for cep in link.crossing_entries() if criteria(cep)])
    strands = []
    while ceps:
        cep = ceps.pop()
        start_crossing = cep.crossing
        is_loop = False
        forward_strand = [cep]
        forward_cep = cep.next()
        while criteria(forward_cep):
            if forward_cep.crossing == start_crossing:
                is_loop = True
                break
            forward_strand.append(forward_cep)
            forward_cep = forward_cep.next()
        backwards_strand = []
        backwards_cep = cep.previous()
        if is_loop:
            strand = forward_strand
        else:
            while criteria(backwards_cep):
                backwards_strand.append(backwards_cep)
                backwards_cep = backwards_cep.previous()
            strand = backwards_strand[::-1]
            strand.extend(forward_strand)

        strands.append(strand)
        ceps.difference_update(strand)

    return sorted(strands, key=len, reverse=True)


def randomize_within_lengths(items):
    by_lens = collections.defaultdict(list)
    for item in items:
        by_lens[len(item)].append(item)
    ans = []
    for length, some_items in sorted(by_lens.items(), reverse=True):
        random.shuffle(some_items)
        ans += some_items
    return ans


def reverse_type_I(link, crossing_strand, label, hand, rebuild=False):
    """
    Add a loop on the strand at crossing_strand with a given label and
    'handedness' hand (twisting left or right).
    """
    D = Crossing(label)
    link.crossings.append(D)
    cs1 = crossing_strand
    cs2 = cs1.opposite()

    if hand == 'left':
        D[1] = D[2]
        cs1ec, cs1cep = cs1.crossing, cs1.strand_index
        D[0] = cs1ec[cs1cep]
        cs2ec, cs2cep = cs2.crossing, cs2.strand_index
        D[3] = cs2ec[cs2cep]
    else:
        D[2] = D[3]
        cs1ec, cs1cep = cs1.crossing, cs1.strand_index
        D[0] = cs1ec[cs1cep]
        cs2ec, cs2cep = cs2.crossing, cs2.strand_index
        D[1] = cs2ec[cs2cep]
    if rebuild:
        link._rebuild(same_components_and_orientations=True)


def random_reverse_type_I(link, label, rebuild=False):
    """
    Randomly adds a loop in a strand, adding one crossing with given label
    """
    cs = random.choice(link.crossing_strands())
    lr = random.choice(['left', 'right'])
    reverse_type_I(link, cs, label, lr, rebuild=rebuild)


def reverse_type_II(link, c, d, label1, label2, rebuild=False):
    """
    Cross two strands defined by two crossing strands c and d in the same face.
    """

    new1, new2 = Crossing(label1), Crossing(label2)
    c_cross, c_ep = c.crossing, c.strand_index
    cop_cross, cop_ep = c.opposite().crossing, c.opposite().strand_index
    d_cross, d_ep = d.crossing, d.strand_index
    dop_cross, dop_ep = d.opposite().crossing, d.opposite().strand_index
    new1[2], new1[3] = new2[0], new2[3]
    new1[0], new1[1] = dop_cross[dop_ep], c_cross[c_ep]
    new2[1], new2[2] = cop_cross[cop_ep], d_cross[d_ep]

    link.crossings.append(new1)
    link.crossings.append(new2)

    if rebuild:
        link._rebuild(same_components_and_orientations=True)


def random_reverse_type_II(link, label1, label2, rebuild=False):
    """
    Randomly crosses two strands, adding two crossings, with labels
    label1 and label2
    """

    faces = link.faces()
    while True:
        face = random.choice(faces)
        if len(face) > 1:
            break
    c, d = random.sample(face, 2)
    reverse_type_II(link, c, d, label1, label2, rebuild=rebuild)


def random_reverse_move(link, t, n):
    """
    Performs a crossing increasing move of type t, where t is 1, 2, or 3
    n is for labeling the new crossings
    """
    if t == 1:
        random_reverse_type_I(link, 'new' + str(n))
    elif t == 2:
        random_reverse_type_II(link, 'new' + str(n), 'new' + str(n + 1))
    else:
        poss_moves = possible_type_III_moves(link)
        if len(poss_moves) != 0:
            reidemeister_III(link, random.choice(poss_moves))


def backtrack(link, num_steps=10, prob_type_1=.3, prob_type_2=.3):
    """
    Randomly perform a series of Reidemeister moves which increase or preserve
    the number of crossings of a link diagram, with the number of such moves
    num_steps.  Can set the probability of type I or type II moves, so the
    remainder is then the probability of type III. Use the method backtrack
    in the Link class.
    """
    if len(link) == 0:
        return

    n = 0
    if prob_type_1 + prob_type_2 > 1:
        raise ValueError("Probabilities add to more than 1")
    p1 = prob_type_1
    p2 = p1 + prob_type_2
    for i in range(num_steps):
        x = random.uniform(0, 1)
        if x < p1:
            t = 1
        elif p1 < x < p2:
            t = 2
        else:
            t = 3

        n += t % 3

        random_reverse_move(link, t, n)

    link._rebuild(same_components_and_orientations=True)


def clear_orientations(link):
    """
    Resets the orientations on the crossings of a link to default values
    """
    link.link_components = None
    for i in link.crossings:
        i.sign = 0
        i.directions.clear()


def relabel_crossings(link):
    """
    Relabel the crossings as integers
    """
    for i, cr in enumerate(link.crossings):
        cr.label = str(i)


def pickup_simplify(link, type_III=0):
    """
    Optimizes the overcrossings on a diagram, then the undercrossings,
    simplifying in between until the process stabilizes.
    """
    L = link
    init_num_crossings = len(L.crossings)
    if type_III:
        simplify_via_level_type_III(link, type_III)
    else:
        basic_simplify(L, build_components=False)
    stabilized = init_num_crossings == 0

    while not stabilized:
        old_cross = len(L.crossings)
        strand_pickup(L, 'over')
        if type_III:
            simplify_via_level_type_III(link, type_III)

        strand_pickup(L, 'under')
        if type_III:
            simplify_via_level_type_III(link, type_III)

        new_cross = len(L.crossings)
        stabilized = new_cross == 0 or new_cross == old_cross

    L._rebuild()
    return len(L.crossings) != init_num_crossings
