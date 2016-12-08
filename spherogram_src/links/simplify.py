"""
Simplifying link diagrams.

Important notes:

* The link diagram is modified in place.  All the relavent parts of the
data structure are updated at each step.  

* Unknot components which are also unlinked may be silently discarded.
"""

from future.utils import iteritems
from .links import Link, Strand, Crossing
from .. import graphs
import random
import networkx as nx

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
        components_removed = len(link.link_components)-len(new_components) 
        link.unlinked_unknot_components += components_removed
        link.link_components = new_components
        
def add_crossings(link, crossings_to_add, tail_dict):
    return

def reidemeister_I(link, C):
    """
    Does a type-1 simplification on the given crossing C if possible.

    Returns the pair: {crossings eliminated}, {crossings changed}
    """
    elim, changed = set(), set()
    for i in range(4):
        if C.adjacent[i] == (C, (i+1)%4):
            (A, a), (B, b) = C.adjacent[i+2], C.adjacent[i+3]
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
            (B, b), (C, c) = A.adjacent[a], A.adjacent[a+1]
            if B == C and (b-1) % 4 == c and (a+b) % 2 == 0:
                eliminated, changed = reidemeister_I(link, B)
                if eliminated:
                    break
                else:
                    W, w = A.adjacent[a+2]
                    X, x = A.adjacent[a+3]
                    Y, y = B.adjacent[b+1]
                    Z, z = B.adjacent[b+2]
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

def basic_simplify(link):
    """
    Do Reidemeister I and II moves until none are possible.
    """
    to_visit, eliminated = set(link.crossings), set()
    while to_visit:
        crossing = to_visit.pop()
        elim, changed = reidemeister_I_and_II(link, crossing)
        assert not elim.intersection(changed)
        eliminated.update(elim)
        to_visit.difference_update(elim)
        to_visit.update(changed)

    success = len(eliminated) > 0

    # Redo the strand labels (used for DT codes)
    if success:
        component_starts = []
        for component in link.link_components:
            a, b = component[:2]
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
                while(face[1][1]% 2 != 0 or face[2][1]% 2 != 1):    # renumber face_list
                    face = [face[1], face[2], face[0]]
                if len(set([e.crossing for e in face])) == 3:  # No repeated crossings
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
    a, b, c =  [t.strand_index for t in triple]
    # We insert Strands around the border of the triple to make the code more
    # transparent and eliminate some special cases.
    old_border =  [(C, c-1), (C, c-2), (A, a-1), (A, a-2), (B, b-1), (B, b-2)]
    border_strands = [insert_strand(*P) for P in old_border]
    new_boarder = [(A,a), (B, b+1), (B, b), (C, c+1), (C, c), (A, a+1)]
    for i, (X,x) in enumerate(new_boarder):
        X[x] = border_strands[i][0]
    A[a-1], B[b-1], C[c-1] = B[b+2], C[c+2], A[a+2]
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
    return list(set(X) & set(Y))[0]

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

        for edge, face in iteritems(to_face):
            neighbor = to_face[edge.opposite()]
            if face.label < neighbor.label:
                dual_edge = self.add_edge(face, neighbor)
                dual_edge.interface = (edge, edge.opposite())
                dual_edge.label= len(self.edges) - 1

        #assert self.is_planar()

    def two_cycles(self):
        """
        Finds all two cycles and returns them as a pair of CrossingStrands which
        are dual to the edges in the cycle.  The crossing strands are
        oriented consistently with respect to one of the faces which a
        vertex for the cycle.
        """
        cycles = []
        for face0 in self.vertices:
            for dual_edge0 in self.incident(face0):
                face1 = dual_edge0(face0) 
                if face0.label < face1.label:
                    for dual_edge1 in self.incident(face1):
                        if dual_edge0.label < dual_edge1.label and dual_edge1(face1) == face0:
                            cycles.append( (common_element(face0, dual_edge0.interface),
                                            common_element(face0, dual_edge1.interface)))
        return cycles

        
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
        edges_crossed.append( graph.edges_between(f1,f2).pop() )

    #want one more edge
    endpoint = overstrand[-1].next()
    final_f1 = graph.edge_to_face[endpoint]
    final_f2 = graph.edge_to_face[endpoint.opposite()]
    edges_crossed.append( graph.edges_between(final_f1,final_f2).pop() )

    return edges_crossed

def extend_overstrand_forward(overstrand,end_cep):
    """
    Starting at a crossing (under or over), extend overstrand until you hit
    the next undercrossing.
    """
    cep = end_cep.next()
    overstrand.append(end_cep)
    while cep.is_over_crossing():
        overstrand.append(cep)
        cep = cep.next()
        if cep == overstrand[0]:
            break

def extend_overstrand_backward(overstrand, start_cep):
    """
    Starting at a crossing (under or over), extend overstrand in the backwards
    direction until you hit another undercrossing.
    """
    cep = start_cep.previous()
    overstrand.insert(0,start_cep)
    while start_cep.is_over_crossing():
        overstrand.insert(0,cep)
        cep = cep.previous()
        if cep == overstrand[-1]:
            break

def pickup_overstrand(link,overstrand):
    """
    Simplify the given overcrossing strand by erasing from the diagram and
    then finding a path that minimizes the number of edges it has to cross
    over to connect the same endpoints. Returns number of crossings removed.
    """
    startcep = overstrand[0].previous()

    if startcep == overstrand[-1]:
        #Totally overcrossing loop, must be totally unlinked and unknotted
        remove_overstrand(overstrand)
        return len(overstrand)
    length = len(overstrand)
    G = link.dual_graph()
    crossing_set = set([cep.crossing for cep in overstrand])
    endpoint = overstrand[-1].next()

    if endpoint.crossing in crossing_set:
        #strand crosses itself underneath
        extend_overstrand_forward(overstrand,endpoint)
        return pickup_overstrand(link,overstrand)
    if startcep.crossing in crossing_set:
        #strand crosses itself over
        extend_overstrand_backward(overstrand,startcep)
        return pickup_overstrand(link,overstrand)
        
    edges_crossed = dual_edges(overstrand, G)

    #create a networkx graph with these edges, and find the connected components
    cross_graph = nx.Graph()
    cross_graph.add_edges_from(edges_crossed)
    components = list(nx.connected_components(cross_graph))

    #collapse the connected components in original dual graph
    Gx = G.to_networkx()
    component_dict = {}
    for i,comp in enumerate(components):
        merge_vertices(Gx,comp)
        for vert in comp:
            component_dict[vert]=tuple(comp)
    Gx_nodes = Gx.nodes()

    #get components which are the start and end nodes in Gx
    source = component_dict[edges_crossed[0][0]]
    dest = component_dict[edges_crossed[-1][0]]

    path = nx.shortest_path(Gx,source,dest)

    crossingsremoved = len(crossing_set) - (len(path) - 1)

    if crossingsremoved == 0:
        return 0

    #force all elements of path to be represented as tuples (to account for single elements)
    for i in range(len(path)):
        if(not type(path[i]) == tuple):
            path[i] = path[i],


    #creating a new list of crossings from which to rebuild the link, remove old overcross
    newcrossings = link.crossings

    cr = remove_overstrand(link, overstrand)
    loose_end = startcep.rotate(2)
    #find new sequence of overcrossings to create
    for i in range(len(path)-1):
        label = 'new%d' % i

        for f1 in path[i]:
            for f2 in path[i+1]:
                if G.edges_between(f1,f2):
                    current_face = f1
                    next_edge = G.edges_between(f1,f2).pop()

        nextedge = next_edge
        if next_edge.interface[0] in current_face:
            cep_to_cross = next_edge.interface[0]
        else:
            cep_to_cross = next_edge.interface[1]

        new_crossing, loose_end = cross_over(link, cep_to_cross, loose_end, label)
        newcrossings.append(new_crossing)

    lec, lecep = loose_end.crossing, loose_end.strand_index
    ec, ecep = endpoint.crossing, endpoint.strand_index
    ec[ecep] = lec[lecep]

    link._rebuild()

    return crossingsremoved


def strand_pickup(link,overstrands):
    """
    Simplifies link by optimizing the path of the longest sequence of overcrossings.
    Returns a new link and the number of crossings removed.
    """
    for overstrand in overstrands:
        if len(overstrand) == 1: continue
        crossings_removed = pickup_overstrand(link,overstrand)
        if crossings_removed != 0:
            return crossings_removed
    return 0


def remove_overstrand(link,overstrand):
    """
    Delete an overstrand from a link.  If the overstrand is a loop, it doesn't
    leave any loose strands and removes the loop. Otherwise, there will be 
    two strands left in the link not attached to anything.  This function
    assumes that the start and end of the overstrand are not places where 
    overstrands crosses itself.
    """
    #only add bridge strands for the places where the strand doesn't cross itself
    crossings_seen = [s.crossing for s in overstrand]
    crossing_set = set()
    for c in crossings_seen:
        if c in crossing_set:
            crossing_set.remove(c)
        else:
            crossing_set.add(c)

    start_cep = overstrand[0].previous()
    end_cep = overstrand[-1].next()
    bridge_strands = dict([(c, Strand('strand'+str(c.label))) for c in crossing_set])
    for cep in overstrand:
        c = cep.crossing
        if c not in crossing_set:
            continue
        strand_index = cep.strand_index
        right_cs = cep.rotate(1).opposite()
        left_cs = cep.rotate(3).opposite()
        if right_cs.crossing in crossing_set:
            signs_equal = (c.sign == right_cs.crossing.sign)
            bridge_strands[c][0] = bridge_strands[right_cs.crossing][signs_equal]
        else:
            bridge_strands[c][0] = right_cs.crossing[right_cs.strand_index]
        if left_cs.crossing in crossing_set:
            signs_equal = (c.sign == left_cs.crossing.sign)
            bridge_strands[c][1] = bridge_strands[left_cs.crossing][1-signs_equal]
        else:
            bridge_strands[c][1] = left_cs.crossing[left_cs.strand_index]
    remove_crossings(link,set(crossings_seen))

    for s in bridge_strands.values():
        s.fuse()

    return len(crossing_set)

def cross_over(link, cep_to_cross, loose_end, label):
    """
    Create a new crossing crossing over the edge defined by cep_to_cross
    and its opposite, and attach one side to a given position loose_end.
    """
    e = cep_to_cross
    new_crossing = Crossing(label)
    lec, lecep = loose_end.crossing, loose_end.strand_index
    new_crossing[1] = lec[lecep]
    ic,icep = e.crossing,e.strand_index
    ico,icoep = e.opposite().crossing, e.opposite().strand_index
    while ic not in link.crossings:
        temp = ic.crossing_strands()[icep]
        ic,icep = temp.next().crossing,temp.next().strand_index
    while ico not in link.crossings:
        temp = ico.crossing_strands()[icoep]
        ico,icoep = temp.next().crossing,temp.next().strand_index
    new_crossing[2] = ic[icep]
    new_crossing[0] = ico[icoep]
    
    return new_crossing, new_crossing.crossing_strands()[3]

def merge_vertices(graph,vertices):
    """
    Merges list of vertices of networkx graph and throws together all
    edges of all the merged vertices.
    """

    v = tuple(vertices)
    graph.add_node(v)
    for i in range(len(v)):
        edgelist = graph.edges(v[i])
        for j in range(len(edgelist)):
            graph.add_edge(v,edgelist[j][0])
            graph.add_edge(v,edgelist[j][1])
        graph.remove_node(v[i])
    return


def reverse_type_I(link,crossing_strand,label,hand,rebuild=False):
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
        comp_sts = [comp[0] for comp in link.link_components]
        link._rebuild(component_starts=comp_sts)

def random_reverse_type_I(link,label,rebuild=False):
    """
    Randomly adds a loop in a strand, adding one crossing with given label
    """
    
    cs = random.choice(link.crossing_strands())
    lr = random.choice(['left','right'])
    reverse_type_I(link,cs,label,lr,rebuild=rebuild)

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
        comp_sts = [comp[0] for comp in link.link_components]
        link._rebuild(component_starts = comp_sts)


def random_reverse_type_II(link, label1, label2, rebuild=False):
    """
    Randomly crosses two strands, adding two crossings, with labels 
    label1 and label2
    """

    faces = link.faces()
    while True:
        face = random.choice(faces)
        if len(face)>1:
            break
    c, d = random.sample(face,2)
    reverse_type_II(link,c,d,label1,label2,rebuild=rebuild)

def random_reverse_move(link,t,n):
    """
    Performs a crossing increasing move of type t, where t is 1, 2, or 3
    n is for labeling the new crossings
    """
    if t == 1:
        random_reverse_type_I(link,'new'+str(n))
    elif t == 2:
        random_reverse_type_II(link,'new'+str(n),'new'+str(n+1))
    else:
        poss_moves = possible_type_III_moves(link)
        if len(poss_moves) != 0:
            reidemeister_III(link, random.choice(poss_moves))


def backtrack(link, num_steps = 10, prob_type_1 = .3, prob_type_2 = .3):
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
        raise Exception("Probabilities add to more than 1")
    p1 = prob_type_1
    p2 = p1 + prob_type_2
    for i in range(num_steps):
        x = random.uniform(0,1)
        if x < p1:
            t = 1
        elif p1 < x < p2:
            t = 2
        else:
            t = 3

        n += t%3

        random_reverse_move(link,t,n)

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
    for i,cr in enumerate(link.crossings):
        cr.label = str(i)

    
def pickup_simplify(link, type_III=0):
    """
    Performs optimize_overcrossings on a diagram, flips, and performs the
    same process on the other side of the diagram, simplifying in between
    until the process stabilizes. The boolean full_simplify indicates 
    whether or not to perform a simplification that includes Reidemeister III
    type moves.
    """
    L = link
    init_num_crossings = len(L.crossings)
    stabilized = init_num_crossings == 0

    def intermediate_simplify(a_link):
        if type_III:
            simplify_via_level_type_III(a_link, type_III)
        else:
            basic_simplify(a_link)

    intermediate_simplify(link)

    while not stabilized:
        overcrossingsremoved = L.optimize_overcrossings()
        intermediate_simplify(L)

        if len(L.crossings) == 0:
            break
        mirror = L.mirror()
        undercrossingsremoved = mirror.optimize_overcrossings()
        L = mirror.mirror()
        intermediate_simplify(L)
        stabilized = ((overcrossingsremoved == 0) and (undercrossingsremoved == 0)) or (len(L.crossings) == 0)

    link.crossings = L.crossings
    link.labels = L.labels
    link.link_components = L.link_components
    link.name = L.name
    return len(L.crossings) != init_num_crossings

