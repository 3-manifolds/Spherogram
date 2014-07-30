"""
Simplifying link diagrams.

Important notes:

* The link diagram is modified in place.  All the relavent parts of the
data structure are updated at each step.  

* Unknot components which are also unlinked may be silently discarded.
"""

from .links import Link, Strand
from .. import graphs
import random

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
                try:
                    component.remove(C)
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

    return len(eliminated) > 0

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

def simplify(link, max_consecutive_failures=100):
    """
    Applies a series of type III moves to the link, simplifying it via type
    I and II moves whenever possible.
    """
    failures, success = 0, False
    if link.basic_simplify():
        success = True
    while failures < max_consecutive_failures:
        poss_moves = possible_type_III_moves(link)
        if len(poss_moves) == 0:
            break
        reidemeister_III(link, random.choice(poss_moves))
        if link.basic_simplify():
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
        
