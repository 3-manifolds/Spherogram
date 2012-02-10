"""
Rational tangles, following

Classifying and Applying Rational Knots and Rational Tangles
http://homepages.math.uic.edu/~kauffman/VegasAMS.pdf

Note:  Tangles are all mutable, operations are typically in place,
and there's currently no way to copy one.  
"""

from link import Crossing, Link, Strand

def join_strands( (a,i), (b,j) ):
    a.adjacent[i] = (b,j)
    b.adjacent[j] = (a,i)
    
class Tangle:
    def __init__(self, crossings=None, entry_points=None):
        if crossings == None:
            crossings = []

        self.crossings = crossings
        self.adjacent = [None, None, None, None]
        if entry_points:
            for i, e in enumerate(entry_points):
                join_strands( (self, i), e)
        
    def __add__(self, other):
        "Join with self to right of other"
        S, O = self.adjacent, other.adjacent
        join_strands( O[3], S[0] )
        join_strands( O[2], S[1] )
        return Tangle(self.crossings + other.crossings, [O[0], O[1], S[2], S[3]])

    def __mul__(self, other):
        "Join with self *above* other"
        S, O = self.adjacent, other.adjacent
        join_strands( O[1], S[0] )
        join_strands( O[2], S[3] )
        return Tangle(self.crossings + other.crossings, [O[0], S[1], S[2], O[3]])

    def __neg__(self):
        "Mirror image of self"
        [c.rotate_by_90() for c in self.crossings if not isinstance(c, Strand)]
        return self

    def invert(self):
        "Rotate anticlockwise by 90 and take the mirror image"
        self.adjacent = self.adjacent[1:] + self.adjacent[:1]
        for i, (o, j) in enumerate(self.adjacent):
            o.adjacent[j] = (self, i)
        return -self

    def numerator_closure(self):
        "The bridge picture closure"
        a, b, c, d = self.adjacent
        join_strands(a, d)
        join_strands(b, c)
        return Link(self.crossings)

    def denominator_closure(self):
        "The bridge picture closure"
        a, b, c, d = self.adjacent
        join_strands(a, b)
        join_strands(c, d)
        return Link(self.crossings)

class ZeroTangle(Tangle):
    def __init__(self):
        bot, top = Strand('B'), Strand('T')
        Tangle.__init__(self, [bot, top], [ (bot, 0), (top, 0), (top, 1), (bot, 1) ] )

class InfinityTangle(Tangle):
    def __init__(self):
        left, right = Strand('L'), Strand('R')
        Tangle.__init__(self, [left, right],  [ (right, 0), (right, 1), (left, 1), (left, 0) ] )

class MinusOneTangle(Tangle):
    def __init__(self):
        c = Crossing('-one')
        Tangle.__init__(self, [c], [(c,i) for i in range(4)])

class OneTangle(Tangle):
    def __init__(self):
        c = Crossing('one')
        Tangle.__init__(self, [c], [(c, (i+1) % 4) for i in range(4)])
    
def integer_tangle(n):
    if n == 0:
        return ZeroTangle()
    if n > 0:
        T = OneTangle()
        for i in range(n - 1):
            T += OneTangle()
        return T
    if n < 0:
        return -integer_tangle(-n)

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

def rational_tangle(a, b):
    T = InfinityTangle()
    for p in reversed(continued_fraction_expansion(a,b)):
        T = integer_tangle(p) + T.invert()
    return T




