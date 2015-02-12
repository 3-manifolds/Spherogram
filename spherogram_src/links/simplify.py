"""
Simplifying link diagrams.

Important notes:

* The link diagram is modified in place.  All the relavent parts of the
data structure are updated at each step.  

* Unknot components which are also unlinked may be silently discarded.
"""

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
        link.link_components = new_components
        
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
            if sum(ce.entry_point % 2 for ce in face) in [1, 2]:
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
    a, b, c =  [t.entry_point for t in triple]
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

        for edge, face in to_face.iteritems():
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
        
        
def strand_pickup(link,overcrossingstrand):
    """
    Simplifies link by optimizing the path of the longest sequence of overcrossings.
    Returns a new link and the number of crossings removed.
    """
    for overcross in overcrossingstrand:
        startcep = overcross[0]
        length = overcross[1]
        G = link.dual_graph()

        #finding all crosses traversed by the overcrossing, accounting for possible self-intersection
        endpoint = startcep.next()
        crossingset = set([endpoint.crossing])
        for i in range(1,length):
            endpoint = endpoint.next()
            crossingset.add(endpoint.crossing)
        endpoint = endpoint.next()

        #creating list of edges of the dual graph corresponding to segments of the strand overcross
        crossgraph = nx.Graph()
        edgescrossed = []
        s = startcep


        for i in range(length+1):
            edge = (s.rotate(2),s.next())
            listOfEdges = list(G.edges)
            for j in range(len(listOfEdges)):
                edgereverse = (edge[1],edge[0])
                if(listOfEdges[j].interface == edge or listOfEdges[j].interface == edgereverse ):
                    edge = listOfEdges[j]
                    break

            edgescrossed.append(edge)
            s = s.next()



        #create a networkx graph with these edges, and find the connected components
        for i in range(len(edgescrossed)):
            crossgraph.add_edge(edgescrossed[i][0],edgescrossed[i][1])


        components = list(nx.connected_components(crossgraph))

        #collapse the connected components in original dual graph
        Gx = G.to_networkx()
        for i in range(len(components)):
            merge_vertices(Gx,components[i])
        Gx_nodes = Gx.nodes()

        #find shortest path between start and end points
        source = None
        dest = None

        for i in range(len(Gx)):
            for j in range(len((Gx_nodes[i]))):

                if Gx_nodes[i][j] == edgescrossed[0][0]:

                    source = Gx_nodes[i]

                if Gx_nodes[i][j] == edgescrossed[-1][0]:

                    dest = Gx_nodes[i]
        path = nx.shortest_path(Gx,source,dest)

        crossingsremoved = length - (len(path) - 1)

        if crossingsremoved == 0:
            continue

        #force all elements of path to be represented as tuples (to account for single elements)
        for i in range(len(path)):
            if(not type(path[i]) == tuple):
                path[i] = path[i],

        #creating a new list of crossings from which to rebuild the link, remove old overcross
        newcrossings = list(link.crossings)
        for i in newcrossings:   #remove old orientations
            i.sign = 0
            i.directions.clear()
        toremove = startcep.next()
        for i in range(len(crossingset)):
            loose1 = toremove.rotate(1).opposite()
            loose2 = toremove.rotate(3).opposite()

            lc1, lc1ep = loose1.crossing, loose1.entry_point
            lc2, lc2ep = loose2.crossing, loose2.entry_point

            while lc1 not in newcrossings:
                lc1, lc1ep = lc1.rotate(2).opposite().crossing, lc1.rotate(2).opposite().entry_point
            while lc2 not in newcrossings:
                lc2, lc2ep = lc2.rotate(2).opposite().crossing, lc2.rotate(2).opposite().entry_point
            lc1[lc1ep] = lc2[lc2ep]
            newcrossings.remove(toremove.crossing)
            toremove = toremove.next()


        looseend = startcep.rotate(2)

        #find new sequence of overcrossings to create
        for i in range(len(path)-1):
            nextedge = None
            label = 'new%d' % i
            crossingtoadd = Crossing(label)
            first = None
            for j in range(len(path[i])):
                found = False
                idict = G.incidence_dict[path[i][j]]
                for l in range(len(idict)):
                    if(idict[l][0] != path[i][j]):
                        totest = idict[l][0]
                    else:
                        totest = idict[l][1]

                    if(totest in path[i+1]):
                        found = True
                        nextedge = idict[l]
                        first = path[i][j]


            for i in first:
                if i == nextedge.interface[0] or i == nextedge.interface[1]:
                    lec, lecep = looseend.crossing, looseend.entry_point
                    crossingtoadd[1] = lec[lecep]
                    ic,icep = i.crossing,i.entry_point
                    ico,icoep = i.opposite().crossing, i.opposite().entry_point
                    while ic not in newcrossings:
                        temp = ic.crossing_strands()[icep]
                        ic,icep = temp.rotate(2).opposite().crossing,temp.rotate(2).opposite().entry_point
                    while ico not in newcrossings:
                        temp = ico.crossing_strands()[icoep]
                        ico,icoep = temp.rotate(2).opposite().crossing,temp.rotate(2).opposite().entry_point

                    crossingtoadd[2] = ic[icep]
                    crossingtoadd[0] = ico[icoep]



            looseend = crossingtoadd.crossing_strands()[3]
            newcrossings.append(crossingtoadd)

        lec, lecep = looseend.crossing, looseend.entry_point
        ec, ecep = endpoint.crossing, endpoint.entry_point
        ec[ecep] = lec[lecep]
        return Link(newcrossings), crossingsremoved

    return link, 0

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


def random_reverse_type_I(link,label):
    """
    Randomly adds a loop in a strand, adding one crossing with given label
    """
    
    cs1 = random.choice(link.crossing_strands())
    D = Crossing(label)
    link.crossings.append(D)

    cs2 = cs1.opposite()
    D[2] = D[3]
    cs1ec, cs1cep = cs1.crossing, cs1.entry_point
    D[0] = cs1ec[cs1cep]
    cs2ec, cs2cep = cs2.crossing, cs2.entry_point
    D[1] = cs2ec[cs2cep]
    
    D.rotate(random.randint(0,1)) #choose whether over or under crossing

def random_reverse_type_II(link, label1, label2):
    """
    Randomly crosses two strands, adding two crossings, with labels label1 and label2
    """

    G = DualGraphOfFaces(link)
    while True:
        face = random.choice(list(G.vertices))
        if len(face)>1:
            break
    c, d = random.sample(face,2)
    new1, new2 = Crossing(label1), Crossing(label2)    
    c_cross, c_ep = c.crossing, c.entry_point
    cop_cross, cop_ep = c.opposite().crossing, c.opposite().entry_point
    d_cross, d_ep = d.crossing, d.entry_point
    dop_cross, dop_ep = d.opposite().crossing, d.opposite().entry_point
    new1[2], new1[3] = new2[0], new2[3]
    new1[0], new1[1] = dop_cross[dop_ep], c_cross[c_ep]
    new2[1], new2[2] = cop_cross[cop_ep], d_cross[d_ep]

    link.crossings.append(new1)
    link.crossings.append(new2)

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


def backtrack(link, num_steps = 10):
    """
    Randomly perform a series of Reidemeister moves which increase or preserve the
    number of crossings of a link diagram, with the number of such moves num_steps
    Use the method backtrack in the Link class.
    """
    if len(link) == 0:
        return link

    n = 0
    for i in range(num_steps):
        t = random.randint(1,3)
        n += t%3

        random_reverse_move(link,t,n)
        
        clear_orientations(link)

        L = Link(link.crossings)
        link = L
    clear_orientations(link)
    relabel_crossings(link)
    return Link(link.crossings)
    

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
        L, overcrossingsremoved = L.optimize_overcrossings()
        intermediate_simplify(L)

        if len(L.crossings) == 0:
            break
        mirror = L.mirror()
        mirror, undercrossingsremoved = mirror.optimize_overcrossings()
        L = mirror.mirror()
        intermediate_simplify(L)
        stabilized = ((overcrossingsremoved == 0) and (undercrossingsremoved == 0)) or (len(L.crossings) == 0)

    link.crossings = L.crossings
    link.labels = L.labels
    link.link_components = L.link_components
    link.name = L.name
    return len(L.crossings) != init_num_crossings


    
