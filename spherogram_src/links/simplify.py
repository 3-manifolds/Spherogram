"""
Simplifying link diagrams.

Important notes:

* The link diagram is modified in place.  All the relavent parts of the
data structure are updated at each step.  

* Unknot components which are also unlinked may be silently discarded.
"""

def remove_crossings(link, eliminate):
    """
    Deletes the given crossings. Assumes that they have already been
    disconnected from the rest of the link, so this just updates
    link.crossings and link.link_components.
    """
    if len(eliminate):
        link.crossings = [C for C in link.crossings if not C in eliminate]
        new_components = []
        for comp in link.link_components:
            new_comp = [ce for ce in comp if ce.crossing not in eliminate]
            if new_comp:
                new_components.append(new_comp)
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
            elim = {C}
            if C != A:
                A[a] = B[b]
                changed = {A, B}

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
                    eliminated = {A, B}
                    if W != B:
                        W[w] = Z[z]
                        changed.update({W, Z})
                    if X != B:
                        X[x] = Y[y]
                        changed.update({X, Y})
                    break

    remove_crossings(link, eliminated)
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

