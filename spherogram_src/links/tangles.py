"""
Rational tangles, following

Classifying and Applying Rational Knots and Rational Tangles
http://homepages.math.uic.edu/~kauffman/VegasAMS.pdf

"""

from .links import Crossing, Strand, Link
#import planar_isotopy
from . import planar_isotopy

try:
    import cPickle as pickle
except ImportError: # Python 3
    import pickle

def join_strands(x, y):
    (a,i), (b,j) = x, y
    a.adjacent[i] = (b,j)
    b.adjacent[j] = (a,i)

def rotate_list(L, s):
    n = len(L)
    return [ L[(i + s) % n] for i in range(n) ]
    
class Tangle(object):
    def __init__(self, n=2, crossings=None, entry_points=None, label=None):
        """
        A tangle is a rectangular block of crossings, with n strands
        coming into the top and bottom.

        These are strands are numbered 0 to n-1 on the bottom,
        and n to 2*n - 1 on the top, both from left to right.  
        """
        if crossings == None:
            crossings = []

        self.crossings = crossings
        self.n = n
        self.adjacent = (2*n)*[None]
        if entry_points:
            for i, e in enumerate(entry_points):
                join_strands( (self, i), e)
        self.label = label
        
    def __add__(self, other):
        "Put self to right of other and fuse inside strands"
        A, B = self.copy(), other.copy()
        n, m = A.n, B.n
        a, b = A.adjacent, B.adjacent
        join_strands( a[n-1], b[0] )
        join_strands( a[2*n - 1], b[m] )
        entry_points = a[:n-1] + b[1:m] + a[n:2*n-1] + b[m+1:]
        return Tangle(n + m - 2, A.crossings + B.crossings, entry_points)

    def __mul__(self, other):
        "Join with self *above* other, as with braid multiplication"
        if self.n != other.n:
            raise ValueError("Tangles must have the same number of strands to multiply them")
        A, B = self.copy(), other.copy()
        a, b = A.adjacent, B.adjacent
        n = A.n
        for i in range(n):
            join_strands(a[n+i], b[i])
        return Tangle(n, A.crossings + B.crossings, a[:n] + b[n:])

    def __neg__(self):
        "Mirror image of self"
        T = self.copy()
        [c.rotate_by_90() for c in T.crossings if not isinstance(c, Strand)]
        return T

    def __or__(self, other):
        "Put self to right if other, no fusing of strands"
        A, B = self.copy(), other.copy()
        n, m = A.n, B.n
        a, b = A.adjacent, B.adjacent
        entry_points = a[:n] + b[:m] + a[n:] + b[m:]
        return Tangle(n + m, A.crossings + B.crossings, entry_points)

    def copy(self):
        return pickle.loads(pickle.dumps(self))

    def rotate(self, s):
        "Rotate anticlockwise by s*90 degrees"
        if self.n != 2:
            raise ValueError("Only n=2 tangles can be rotated")
        anticlockwise = [0, 1, 3, 2]
        rotate = dict(zip(anticlockwise, rotate_list(anticlockwise, s)))
        T = self.copy()
        T.adjacent = [T.adjacent[rotate[i]] for i in range(4)]
        for i, (o, j) in enumerate(T.adjacent):
            o.adjacent[j] = (T, i)
        return T
        
    def invert(self):
        "Rotate anticlockwise by 90 and take the mirror image"
        if self.n != 2:
            raise ValueError("Only n=2 tangles can be inverted")
        return -self.rotate(1)

    def numerator_closure(self):
        "The bridge picture closure"
        if self.n % 2 != 0:
            raise ValueError("To do bridge closure, n must be even")
        T = self.copy()
        for i in range(0, 2*self.n, 2):
            join_strands(T.adjacent[i], T.adjacent[i + 1])
        return Link(T.crossings, check_planarity=False)

    def denominator_closure(self):
        "The braid closure picture"
        T = self.copy()
        for i in range(0, self.n):
            join_strands(T.adjacent[i], T.adjacent[i + self.n])
        return Link(T.crossings, check_planarity=False)


    def circular_rotate(self,n):
        """
        Rotate a tangle in a circular fashion.
        """
        tangle_copy = self.copy()
        tangle_copy._fuse_strands()
        adj = tangle_copy.adjacent
        #reverse second half
        adj[len(adj)/2:] = reversed(adj[len(adj)/2:])
        rotated_adj = rotate_list(adj,n)
        #undo reversal of second half
        rotated_adj[len(rotated_adj)/2:] = reversed(rotated_adj[len(rotated_adj)/2:])
        tangle_copy.adjacent = rotated_adj
        return tangle_copy

    def circular_sum(self,other,n):
        """
        Glue two tangles together.  There are many ways to do this.  Choice
        is given by the choice of integer n
        """
        if len(self.adjacent) != len(other.adjacent):
            raise Exception("Tangles do not have the same number of strands")
        return (self*(other.circular_rotate(n))).denominator_closure()

    def isosig(self,root = None, over_or_under = False):
        """
        Return a bunch of data which encodes the planar isotopy class of the
        tangle.  Of course, this is just up to isotopy of the plane
        (no Reidemeister moves).  A root can be specified with a CrossingStrand
        and over_or_under toggles whether only the underlying shadow (4-valent
        planar map) is considered or the tangle with the over/under data at 
        each crossing.
        """
        copy = self.copy()
        copy._fuse_strands()
        return planar_isotopy.min_isosig(copy, root, over_or_under)

    def is_planar_isotopic(self, other, root = None, over_or_under = False):
        return self.isosig() == other.isosig()

    def _fuse_strands(self):
        for s in reversed(self.crossings):
            if isinstance(s,Strand):
                s.fuse()
                self.crossings.remove(s)

    def __repr__(self):
        return "<Tangle: %s>" % self.label

Tangle.bridge_closure = Tangle.numerator_closure
Tangle.braid_closure = Tangle.denominator_closure

class ZeroTangle(Tangle):
    def __init__(self):
        bot, top = Strand('B'), Strand('T')
        Tangle.__init__(self, 2, [bot, top], [ (bot, 0), (bot, 1), (top, 0), (top, 1) ] )

class InfinityTangle(Tangle):
    def __init__(self):
        left, right = Strand('L'), Strand('R')
        Tangle.__init__(self, 2, [left, right],  [ (left, 0), (right, 0), (left, 1), (right, 1) ] )

class MinusOneTangle(Tangle):
    def __init__(self):
        c = Crossing('-one')
        Tangle.__init__(self, 2, [c], [(c,3), (c, 0), (c, 2), (c, 1)])

class OneTangle(Tangle):
    def __init__(self):
        c = Crossing('one')
        Tangle.__init__(self, 2, [c], [(c,0), (c, 1), (c, 3), (c, 2)])
    
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

        Tangle.__init__(self, 2, T.crossings, T.adjacent)
    
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
    def __init__(self, a, b=1):
        if b == 1 and hasattr(a, 'numerator') and hasattr(a, 'denominator') and not isinstance(a, int):
            a, b = a.numerator(), a.denominator()
        if b < 0:
            a, b = -a, -b
        self.fraction = (a,b)
        self.partial_quotients = pqs = continued_fraction_expansion(abs(a),b)
        T = InfinityTangle()
        for p in reversed(pqs):
            T = IntegerTangle(p) + T.invert()
        if a < 0:
            T = -T
        Tangle.__init__(self, 2, T.crossings, T.adjacent)

#----------------------------------------------------
#
# Basic braids
#
#----------------------------------------------------

class IdentityBraid(Tangle):
    def __init__(self, n):
        strands = [Strand() for i in range(n)]
        entry_points = [ (s, 0) for s in strands] + [(s,1) for s in strands]
        Tangle.__init__(self, n, strands, entry_points)
        
            
            
