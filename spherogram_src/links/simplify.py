"""
Simplifying link diagrams.

Important notes:

* The link diagram is modified in place.  All the relavent parts of the
data structure are updated at each step.  

* Unknot components which are also unlinked may be silently discarded.
"""

from . import links
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
    """
    ans = []
    for face in link.faces():
        if len(face) == 3:
            if sum(ce.entry_point % 2 for ce in face) in [1, 2]:
                while(face[1][1]% 2 != 0 or face[2][1]% 2 != 1):    # renumber face_list
                    face = [face[1], face[2], face[0]]
                ans.append(face)
    return ans

def insert_strand(X, x):
    Y, y = X.adjacent[x]
    S = links.Strand()
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



