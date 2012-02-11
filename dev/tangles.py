"""
Rational tangles, following

Classifying and Applying Rational Knots and Rational Tangles
http://homepages.math.uic.edu/~kauffman/VegasAMS.pdf

Note:  Tangles are all mutable, operations are typically in place,
and there's currently no way to copy one.  
"""

import link, copy

def join_strands( (a,i), (b,j) ):
    a.adjacent[i] = (b,j)
    b.adjacent[j] = (a,i)
    
class Tangle:
    def __init__(self, crossings=None, entry_points=None, label=None):
        if crossings == None:
            crossings = []

        self.crossings = crossings
        self.adjacent = [None, None, None, None]
        if entry_points:
            for i, e in enumerate(entry_points):
                join_strands( (self, i), e)
        self.label = label
        
    def __add__(self, other):
        "Join with self to right of other"
        A, B = self.copy(), other.copy()        
        a, b = A.adjacent, B.adjacent
        join_strands( b[3], a[0] )
        join_strands( b[2], a[1] )
        return Tangle(A.crossings + B.crossings, [b[0], b[1], a[2], a[3]])

    def __mul__(self, other):
        "Join with self *above* other"
        A, B = self.copy(), other.copy()
        a, b = A.adjacent, B.adjacent
        join_strands( b[1], a[0] )
        join_strands( b[2], a[3] )
        return Tangle(A.crossings + B.crossings, [b[0], a[1], a[2], b[3]])

    def __neg__(self):
        "Mirror image of self"
        T = self.copy()
        [c.rotate_by_90() for c in T.crossings if not isinstance(c, link.Strand)]
        return T

    def copy(self):
        return copy.deepcopy(self)
    
    def invert(self):
        "Rotate anticlockwise by 90 and take the mirror image"
        T = self.copy()
        T.adjacent = T.adjacent[1:] + T.adjacent[:1]
        for i, (o, j) in enumerate(T.adjacent):
            o.adjacent[j] = (T, i)
        return -T

    def numerator_closure(self):
        "The bridge picture closure"
        a, b, c, d = self.adjacent
        join_strands(a, d)
        join_strands(b, c)
        T = self.copy()
        return link.Link(T.crossings)

    def denominator_closure(self):
        "The bridge picture closure"
        a, b, c, d = self.adjacent
        join_strands(a, b)
        join_strands(c, d)
        T = self.copy()
        return link.Link(T.crossings)

    def __repr__(self):
        return "<Tangle: %s>" % self.label

class ZeroTangle(Tangle):
    def __init__(self):
        bot, top = link.Strand('B'), link.Strand('T')
        Tangle.__init__(self, [bot, top], [ (bot, 0), (top, 0), (top, 1), (bot, 1) ] )

class InfinityTangle(Tangle):
    def __init__(self):
        left, right = link.Strand('L'), link.Strand('R')
        Tangle.__init__(self, [left, right],  [ (right, 0), (right, 1), (left, 1), (left, 0) ] )

class MinusOneTangle(Tangle):
    def __init__(self):
        c = link.Crossing('-one')
        Tangle.__init__(self, [c], [(c,i) for i in range(4)])

class OneTangle(Tangle):
    def __init__(self):
        c = link.Crossing('one')
        Tangle.__init__(self, [c], [(c, (i+1) % 4) for i in range(4)])
    
class IntegerTangle(Tangle):
    def __init__(self, n):
        if n == 0:
            T = ZeroTangle()
        if n > 0:
            T = OneTangle()
            for i in range(n - 1):
                T += OneTangle()
        if n < 0:
            T = -IntegerTangle(-n)

        Tangle.__init__(self, T.crossings, T.adjacent)
    
def continued_fraction_expansion(a, b):
    """
    The continued fraction expansion of a/b. 
    """
    if b == 0: return []
    if b == 1: return [a]
    if b < 0:
        return continued_fraction_expansion(-a, -b)
    q, r = a//b, a % b
    if a < 0:
        return [q] + continued_fraction_expansion(r, b)[1:]
    return [q] + continued_fraction_expansion(b, r)

class RationalTangle(Tangle):
    def __init__(self, a, b):
        self.fraction = (a,b)
        self.partial_quotients = pqs = continued_fraction_expansion(a,b)
        T = InfinityTangle()
        for p in reversed(pqs):
            T = IntegerTangle(p) + T.invert()
        Tangle.__init__(self, T.crossings, T.adjacent)
