"""
This file implement's Dror Bar-Natan's fast tangle-based algorithm
for computing the Alexander polynomial of a knot in S^3, as described
in

http://www.math.toronto.edu/drorbn/Talks/Aarhus-1507/

with some additional details drawn from conversations with him at said
conference.
"""

import random
from .links_base import CrossingEntryPoint
from ..sage_helper import _within_sage

if _within_sage:
    from sage.all import (PolynomialRing, matrix, block_diagonal_matrix, ZZ, QQ)

class ClosedComponentCreated(Exception):
    """
    We never want to create closed link components; the final answer
    is derived from when one has a created the "long knot" string
    link.
    """

def cep(crossing_strand):
    """
    Returns the CrossingEntryPoint corresponding to the given
    CrossingStrand in the same crossing; that is, it orients the
    CrossingStrand without changing the crossing.
    """
    if crossing_strand == crossing_strand.oriented():
        return crossing_strand
    else:
        return crossing_strand.rotate(2)

def entry_pts_ab(crossing):
    """
    The two entry points of a crossing with Dror's convention that the
    overcrossing ("a") is first and the undercrossing ("b") is second. 
    """
    verts = [1, 0] if crossing.sign == -1 else [3, 0]
    return [CrossingEntryPoint(crossing, v) for v in verts]

class StrandIndices(object):
    """
    A map from the crossings strands of a link L onto range(n).
    """
    def __init__(self, link, ordered_crossings=None):
        indices = dict()
        if ordered_crossings is None:
            ordered_crossings = link.crossings
        for i, C in enumerate(ordered_crossings):
            for j, ce in enumerate(entry_pts_ab(C)):
                indices[ce] = 2*i + j
        self.indices = indices

    def merge(self, a, b):
        if self[a] > self[b]:
            a, b = b, a
        i, j = self[a], self[b]
        if i == j:
            raise ClosedComponentCreated
        self[b] = i
        for c, k in self.indices.items():
            if k == j:
                self.indices[c] = i
            elif k > j:
                self.indices[c] = k - 1

    def __getitem__(self, cs):
        return self.indices[cep(cs)]

    def __setitem__(self, cs, val):
        self.indices[cep(cs)] = val

def strand_matrix_merge(A, a, b):
    """
    The main computations all happen in this method.  Here A is a
    square matrix and a and b are row (equivalently column) indices.
    """
    assert a != b
    alpha, beta = A[a, a], A[a, b]
    gamma, delta = A[b, a], A[b, b]
    mu = 1 - beta
    theta, epsilon = A.row(a), A.row(b)
    phi, psi = A.column(a), A.column(b)
    A = A + matrix(psi).transpose()*matrix(theta)/mu
    i, j = min(a, b), max(a, b)
    A[i] = epsilon + delta*theta/mu
    A[:, i] = phi + alpha*psi/mu
    A[i, i] = gamma + alpha*delta/mu
    A = A.delete_rows([j])
    A = A.delete_columns([j])
    return A

def test_meta_associativity():
    """
    Tests strand_matrix_merge for required invariance properties. 
    """
    def eval_merges(merges):
        R = PolynomialRing(QQ, 'x', 49).fraction_field()
        A = matrix(7, 7, R.gens())
        for a, b in merges:
            A = strand_matrix_merge(A, a, b)
        return A

    associative_merges =  [
        ([(0,1), (0,1)], [(1,2), (0,1)]),
        ([(0,1), (0,2)], [(1,3), (0,1)]),
        ([(2,5), (2,3)], [(5,3), (2,3)])
        ]
    for m1, m2 in associative_merges:
        assert eval_merges(m1) == eval_merges(m2)
        
class DrorDatum(object):
    """
    The (omega, A) pair which is the invariant defined in the first column of 
    http://www.math.toronto.edu/drorbn/Talks/Aarhus-1507/
    """
    def __init__(self, link, ordered_crossings):
        self.strand_indices = StrandIndices(link, ordered_crossings)
        self.ring = R = PolynomialRing(ZZ, 't').fraction_field()
        self.omega = R.one()        
        self.A = matrix(R, 0, 0)

    def add_crossing(self, crossing):
        indices = self.strand_indices
        t = self.ring.gen()
        a, b = entry_pts_ab(crossing)
        n = self.A.nrows()
        assert indices[a] == n and indices[b] == n + 1
        T = t if crossing.sign == 1 else t**-1
        B = matrix([[1, 1 - T], [0, T]])
        self.A = block_diagonal_matrix([self.A, B])

    def merge(self, cs_a, cs_b):
        indices, A = self.strand_indices, self.A
        a, b = indices[cs_a], indices[cs_b]
        if a == b:
            raise ClosedComponentCreated
        mu = 1 - A[a, b]
        self.omega *= mu
        self.A = strand_matrix_merge(A, a, b)
        indices.merge(cs_a, cs_b)

def num_overlap(crossing, frontier):
    neighbor_strands = set([cs.opposite() for cs in crossing.crossing_strands()])
    return len(neighbor_strands.intersection(frontier))
        
class Exhaustion(object):
    """
    An exhaustion of a link where crossings are added in one-by-one
    so that the resulting tangle is connected at every stage.  

    Starting at the given crossing, it uses a greedy algorithm to try
    to minimize the sizes of the frontiers of the intermediate tangles.

    If no initial crossing is specified, one is choosen at random.
    """
    def __init__(self, link, crossing=None):
        if crossing is None:
            crossing = random.choice(link.crossings)
        crossings = [crossing]
        gluings = [[]]
        frontier = set(crossing.crossing_strands())
        frontier_lengths = [4]
        while len(crossings) < len(link.crossings):
            choices = [(num_overlap(cs.opposite()[0], frontier), cs) for cs in frontier]
            overlap, cs = max(choices)
            C = cs.opposite().crossing
            assert C not in crossings
            crossings.append(C)
            C_gluings = []
            for cs in C.crossing_strands():
                opp = cs.opposite()
                if opp in frontier:
                    frontier.discard(opp)
                    b = cs.oriented()
                    a = b.opposite()
                    C_gluings.append((a, b))
                else:
                    frontier.add(cs)

            assert frontier_lengths[-1] + 4 - 2*overlap == len(frontier)
            frontier_lengths.append(len(frontier))
            gluings.append(C_gluings)

        self.link = link
        self.crossings = crossings
        self.frontier_lengths = frontier_lengths
        self.gluings = gluings
        self.width = max(frontier_lengths)//2

    def test_indices(self):
        indices = StrandIndices(self.link, self.crossings)
        all_gluings = sum(self.gluings, [])[:-1]
        for a, b in all_gluings[:-1]:
            indices.merge(a, b)
            
    def alexander_polynomial(self):
        D = DrorDatum(self.link, self.crossings)
        gluings = self.gluings[:]
        # The last crossing is special because we need to end at a
        # string link.
        for C, gluings in zip(self.crossings, gluings)[:-1]:
            D.add_crossing(C)
            for a, b in gluings:
                D.merge(a, b)
    
        C = self.crossings[-1]
        D.add_crossing(C)
        for a, b in self.gluings[-1][:-1]:
            D.merge(a, b)

        alex = D.omega
        p, q = alex.numerator(), alex.denominator()
        # make sure the denominator is +/- t^n
        assert [abs(c) for c in q.coefficients()] == [1]
        if p.leading_coefficient() < 0:
            p = -p
        t, e = p.parent().gen(), min(p.exponents())
        return p//t**e
        
def good_exhaustion(link, max_failed_tries=20):
    """
    Uses a random search to try to find an Exhaustion with small width. 
    """
    crossings = list(link.crossings)
    C = random.choice(crossings)
    crossings.remove(C)
    E_best = Exhaustion(link, C)
    tries = 0
    while len(crossings) and tries < max_failed_tries:
        C = random.choice(crossings)
        crossings.remove(C)
        E = Exhaustion(link, C)
        if E.width < E_best.width:
            E_best = E
        tries += 1
    return E_best

def alexander(K):
    c = len(K.crossings)
    if c < 100:
        E = Exhaustion(K)
    else:
        E = good_exhaustion(K, max(20, 0.15*c))
    return E.alexander_polynomial()
