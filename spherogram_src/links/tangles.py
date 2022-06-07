"""
Rational tangles, following

Classifying and Applying Rational Knots and Rational Tangles
http://homepages.math.uic.edu/~kauffman/VegasAMS.pdf

"""
import pickle

from .links import Crossing, Strand, Link
from . import planar_isotopy


def join_strands(x, y):
    """
    Takes two (c, i) pairs where c is a Crossing, Strand, or Tangle and i is an index into
    c.adjacent and joins them by having them refer to each other at those positions.

    Having c be a Tangle is conceptually a special case since its c.adjacent is being
    used to record the boundary strands.

    This is the same as creating a Strand s with s.adjacent = [x, y] and then
    doing s.fuse()
    """
    (a, i), (b, j) = x, y
    a.adjacent[i] = (b, j)
    b.adjacent[j] = (a, i)


def rotate_list(L, s):
    """Rotate the list, putting L[s] into index 0."""
    n = len(L)
    return [L[(i + s) % n] for i in range(n)]

def decode_boundary(boundary):
    """The boundary is either an nonnegative integer or a pair of non-negative integers.

    If it is an integer n, then returns (n, n). If it is a pair (m, n) then returns it.
    """
    if isinstance(boundary, tuple):
        m, n = boundary
    else:
        m = n = boundary
    if m < 0: raise ValueError("Number of bottom boundary strands cannot be negative")
    if n < 0: raise ValueError("Number of top boundary strands cannot be negative")
    return (m, n)

class Tangle():
    def __init__(self, boundary=2, crossings=None, entry_points=None, label=None):
        """
        A tangle is a fragment of a Link with some number of boundary
        strands. Tangles can be composed in various ways along their boundary strands,
        including the horizontal and vertical compositions of the tangle category.

        Inputs:
        - When boundary is an integer, then the tangle has n strands coming into both
          the top and the bottom of the tangle. When boundary is a pair of integers
          (m, n), then the tangle has m strands coming into the bottom and n coming
          into the top.

          The strands are numbered 0 to m-1 on the bottom and m to m+n-1 on the
          top, both from left to right.

        - crossings is a list of Crossing or Strand that comprise the tangle.
        - entry_points is a list of pairs (c, i) where c is a Crossing or Strand
          and i indexes into c.adjacent. These pairs describe the boundary strands
          in order of the strand numbering.
        - label is an arbitrary label for the tangle.
        """

        m, n = decode_boundary(boundary)

        if crossings is None:
            crossings = []
        for c in crossings:
            if not isinstance(c, (Crossing, Strand)):
                raise ValueError("Every element of crossings must be a Crossing or a Strand")
        self.crossings = crossings

        # the pair for the number of lower strands and the number of upper strands
        self.boundary = (m, n)

        # a list of (c, i) pairs for the boundary strands. Each c will reciprocally
        # contain (self, j) where j is the strand number. The fact this is called
        # 'adjacent' means that the Tangle can take part in the joining protocol
        # implemented in join_strands.
        self.adjacent = (m + n) * [None]
        entry_points = entry_points or []
        if len(entry_points) != m + n:
            raise ValueError("The number of boundary strands is not equal to the length"
                             " of entry_points")
        for i, e in enumerate(entry_points):
            join_strands((self, i), e)

        self.label = label

    def __add__(self, other):
        "Put self to left of other and fuse inside strands"
        A, B = self.copy(), other.copy()
        (mA, nA), (mB, nB) = A.boundary, B.boundary
        if mA == 0 or mB == 0 or nA == 0 or nB == 0:
            raise ValueError("Tangles must have at least one top and bottom strand each.")
        a, b = A.adjacent, B.adjacent
        join_strands(a[mA - 1], b[0])
        join_strands(a[mA + nA - 1], b[mB])
        entry_points = a[:mA - 1] + b[1:mB] + a[mA:mA + nA - 1] + b[mB + 1:]
        return Tangle((mA + mB - 2, nA + nB - 2), A.crossings + B.crossings, entry_points)

    def __mul__(self, other):
        "Join with self *above* other, as with braid multiplication"
        A, B = self.copy(), other.copy()
        (mA, nA), (mB, nB) = A.boundary, B.boundary
        if mA != nB:
            raise ValueError("Tangles must have a compatible number of strands to multiply them")

        a, b = A.adjacent, B.adjacent
        for i in range(mA):
            join_strands(a[i], b[mB + i])
        return Tangle((mB, nA), A.crossings + B.crossings, b[:mB] + a[mA:])

    def __neg__(self):
        "Mirror image of self"
        T = self.copy()
        for c in T.crossings:
            if not isinstance(c, Strand):
                c.rotate_by_90()
        return T

    def __or__(self, other):
        "Put self to left of other, no fusing of strands"
        A, B = self.copy(), other.copy()
        (mA, nA), (mB, nB) = A.boundary, B.boundary
        a, b = A.adjacent, B.adjacent
        entry_points = a[:mA] + b[:mB] + a[mA:] + b[mB:]
        return Tangle((mA + mB, nA + nB), A.crossings + B.crossings, entry_points)

    def copy(self):
        return pickle.loads(pickle.dumps(self))

    def rotate(self, s):
        "Rotate anticlockwise by s*90 degrees"
        if self.boundary != (2, 2):
            raise ValueError("Only boundary=(2,2) tangles can be rotated")
        anticlockwise = [0, 1, 3, 2]
        rotate = dict(zip(anticlockwise, rotate_list(anticlockwise, s)))
        T = self.copy()
        T.adjacent = [T.adjacent[rotate[i]] for i in range(4)]
        for i, (o, j) in enumerate(T.adjacent):
            o.adjacent[j] = (T, i)
        return T

    def invert(self):
        "Rotate anticlockwise by 90 and take the mirror image"
        if self.boundary != (2,2):
            raise ValueError("Only boundary=(2,2) tangles can be inverted")
        return -self.rotate(1)

    def numerator_closure(self):
        "The bridge picture closure"
        m, n = self.boundary
        if m % 2 or n % 2:
            raise ValueError("To do bridge closure, both the top and bottom must have an even number of strands")
        T = self.copy()
        for i in range(0, m, 2):
            join_strands(T.adjacent[i], T.adjacent[i + 1])
        for j in range(0, n, 2):
            join_strands(T.adjacent[m + i], T.adjacent[m + i + 1])
        return Link(T.crossings, check_planarity=False)

    def denominator_closure(self):
        "The braid closure picture"
        m, n = self.boundary
        if m != n:
            raise ValueError("To do braid closure, both the top and bottom number of strands must be equal")
        T = self.copy()
        for i in range(0, n):
            join_strands(T.adjacent[i], T.adjacent[m + i])
        return Link(T.crossings, check_planarity=False)

    def link(self):
        "Get the tangle as a link if its boundary is (0, 0)."
        if self.boundary != (0, 0):
            raise ValueError("The boundary must be (0, 0)")
        return Link(self.copy().crossings, check_planarity=False)

    def reshape(self, boundary, displace=0):
        """Renumber the boundary strands so that the tangle has the new boundary
        shape. This is performed by either repeatedly moving the last strands from the
        bottom right to the top right or vice versa. Simultaneously, displace controls
        a rotation of the tangle where the tangle is rotated clockwise 'displace' units
        (so, for example, if 0 <= displace < m then that strand number becomes the new
        lower-left strand).
        """
        m, n = self.boundary
        Tm, Tn = decode_boundary(boundary)
        if (m, n) == (Tm, Tn):
            return self
        if m + n != Tm + Tn:
            raise ValueError("Reshaping requires the tangle have the same number of boundary"
                             " strands as in the new boundary.")
        T = self.copy()
        # The 'adjacent' array but in total counterclockwise order
        adj_ccw = T.adjacent[:m] + list(reversed(T.adjacent[m:]))
        adj_ccw = rotate_list(adj_ccw, displace)

        return Tangle((Tm, Tn), T.crossings,
                      adj_ccw[:Tm] + list(reversed(adj_ccw[Tm:])))

    def circular_rotate(self, n):
        """
        Rotate a tangle in a circular fashion clockwise, keeping the same boundary.
        """
        return self.reshape(self.boundary, n)

    def circular_sum(self, other, n=0):
        """
        Glue two tangles together to form a link by gluing them vertically and then taking
        the braid closure (the denominator_closure).
        The second tangle is rotated clockwise by n strands.
        """
        Am, An = self.boundary
        Bm, Bn = self.boundary
        if (Am, An) != (Bn, Bm):
            raise Exception("Tangles must have a compatible number of strands")
        return (self * (other.circular_rotate(n))).denominator_closure()

    def isosig(self, root=None, over_or_under=False):
        """
        Return a bunch of data which encodes the planar isotopy class of the
        tangle.  Of course, this is just up to isotopy of the plane
        (no Reidemeister moves).  A root can be specified with a CrossingStrand
        and ``over_or_under`` toggles whether only the underlying
        shadow (4-valent planar map) is considered or the tangle with the
        over/under data at each crossing.
        """
        copy = self.copy()
        copy._fuse_strands()
        return planar_isotopy.min_isosig(copy, root, over_or_under)

    def is_planar_isotopic(self, other, root=None, over_or_under=False) -> bool:
        return self.isosig() == other.isosig()

    def _fuse_strands(self, preserve_boundary=False):
        """Fuse all strands and delete them, even ones incident only to the boundary (unless
        preserve_boundary is True)."""
        for s in reversed(self.crossings):
            if isinstance(s, Strand):
                # check that the strand is not only incident to the boundary
                if preserve_boundary and all(a[0] == self for a in s.adjacent):
                    continue
                s.fuse()
                self.crossings.remove(s)

    def __repr__(self):
        return "<Tangle: %s>" % self.label

    def describe(self, fuse_strands=False):
        """Give a PD-like description of the tangle in the form
        Tangle[{lower strands}, {upper strands}, P and X codes].

        If fuse_strands is True, then fuse all internal Strand nodes first."""
        T = self.copy()
        if fuse_strands:
            T._fuse_strands(preserve_boundary = True)
        T.label = 0
        # give each crossing/strand a unique identifier
        for i, c in enumerate(T.crossings):
            c.label = i + 1
        arc_ids = {}
        parts = []
        lower = None
        upper = None
        for c in [T] + T.crossings:
            arcs = []
            for i, (d, j) in enumerate(c.adjacent):
                arc = tuple(sorted([(c.label, i), (d.label, j)]))
                # get arc id or assign a fresh one
                arc_id = arc_ids.setdefault(arc, len(arc_ids) + 1)
                arcs.append(arc_id)
            if isinstance(c, Crossing):
                parts.append("X[%s,%s,%s,%s]" % tuple(arcs))
            elif isinstance(c, Strand):
                parts.append("P[%s,%s]" % tuple(arcs))
            elif isinstance(c, Tangle):
                lower = "{" + ",".join(str(a) for a in arcs[:T.boundary[0]]) + "}"
                upper = "{" + ",".join(str(a) for a in arcs[T.boundary[0]:]) + "}"
            else:
                raise Exception("Unexpected entity")
        parts.sort()
        return f"Tangle[{lower}, {upper}{''.join(', ' + p for p in parts)}]"

Tangle.bridge_closure = Tangle.numerator_closure
Tangle.braid_closure = Tangle.denominator_closure

def CapTangle():
    """The unknotted (2,0) tangle."""
    cap = Strand("cap")
    return Tangle((2, 0), [cap], [(cap, 0), (cap, 1)])

def CupTangle():
    """The unknotted (0,2) tangle."""
    cup = Strand("cup")
    return Tangle((0, 2), [cup], [(cup, 0), (cup, 1)])

def ZeroTangle():
    bot, top = Strand('B'), Strand('T')
    return Tangle(2, [bot, top],
                  [(bot, 0), (bot, 1), (top, 0), (top, 1)],
                  "ZeroTangle")

def InfinityTangle():
    left, right = Strand('L'), Strand('R')
    return Tangle(2, [left, right],
                  [(left, 0), (right, 0), (left, 1), (right, 1)],
                  "InfinityTangle")

def MinusOneTangle():
    c = Crossing('-one')
    return Tangle(2, [c], [(c, 3), (c, 0), (c, 2), (c, 1)],
                  "MinusOneTangle")

def OneTangle():
    c = Crossing('one')
    return Tangle(2, [c], [(c, 0), (c, 1), (c, 3), (c, 2)],
                  "OneTangle")

def IntegerTangle(n):
    if n == 0:
        return ZeroTangle()
    elif n > 0:
        T = OneTangle()
        for i in range(n - 1):
            T += OneTangle()
        T.label = f"IntegerTangle({n})"
        return T
    elif n < 0:
        T = -IntegerTangle(-n)
        T.label = f"IntegerTangle({n})"
        return T
    else:
        raise ValueError("Expecting int")

def continued_fraction_expansion(a, b):
    """
    The continued fraction expansion of a/b.
    """
    if b == 0:
        return []
    if b == 1:
        return [a]
    if b < 0:
        return continued_fraction_expansion(-a, -b)
    q, r = a // b, a % b
    if a < 0:
        return [q] + continued_fraction_expansion(r, b)[1:]
    return [q] + continued_fraction_expansion(b, r)


class RationalTangle(Tangle):
    def __init__(self, a, b=1):
        if b == 1 and hasattr(a, 'numerator') and hasattr(a, 'denominator') and not isinstance(a, int):
            a, b = a.numerator(), a.denominator()
        if b < 0:
            a, b = -a, -b
        self.fraction = (a, b)
        self.partial_quotients = pqs = continued_fraction_expansion(abs(a), b)
        T = InfinityTangle()
        for p in reversed(pqs):
            T = IntegerTangle(p) + T.invert()
        if a < 0:
            T = -T
        Tangle.__init__(self, 2, T.crossings, T.adjacent,
                        f"RationalTangle({a}, {b})")

# ---------------------------------------------------
#
# Basic braids
#
# ---------------------------------------------------


def IdentityBraid(n):
    """
    The braid with n strands and no crossings
    """
    if n < 0:
        raise ValueError("Expecting non-negative int")
    strands = [Strand() for i in range(n)]
    entry_points = [(s, 0) for s in strands] + [(s, 1) for s in strands]
    return Tangle(n, strands, entry_points,
                  f"IdentityBraid({n})")
