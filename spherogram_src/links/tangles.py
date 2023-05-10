"""
A tangle is piece of a knot diagram in a disk where some of the
strands meet the boundary. Tangles can be composed by gluing them
along arcs in each boundary that have the same number of incident
strands.

This module gives a version of tangles where there are four distinguished
boundary arcs used for gluing: the bottom and top, which can have incident
strands, and the left and right, which cannot. Tangles can be glued
vertically using ``*`` and horizontally using ``|``. There is also a second
kind of horizontal composition using ``+`` where the rightmost strands of the top
and bottom of the first tangle are glued to the leftmost strands of the top and
bottom of the second tangle.

Rational tangles (created using ``RationalTangle``) are following the paper

Classifying and Applying Rational Knots and Rational Tangles
http://homepages.math.uic.edu/~kauffman/VegasAMS.pdf

See doc.pdf for conventions.
"""
import pickle

from .links import Crossing, Strand, Link
from . import planar_isotopy


def join_strands(x, y):
    """
    Input: two (c, i) pairs where c is a Crossing, Strand, or Tangle object and i is an index into
    c.adjacent. Joins the objects by having them refer to each other at those positions.

    When c is a Tangle it is conceptually a special case since its c.adjacent is being
    used to record the boundary strands.

    This function equivalent to creating a Strand s with s.adjacent = [x, y] and then
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
    """The boundary is either a nonnegative integer or a pair of non-negative integers.

    * When the input is an integer n, this returns (n, n).
    * When the input is a pair (m, n), then it returns (m, n).

    >>> decode_boundary(2)
    (2, 2)
    >>> decode_boundary((3,4))
    (3, 4)
    >>> decode_boundary(-2)
    Traceback (most recent call last):
        ...
    ValueError: Number of bottom boundary strands cannot be negative
    >>> decode_boundary((3,-2))
    Traceback (most recent call last):
        ...
    ValueError: Number of top boundary strands cannot be negative
    """
    if isinstance(boundary, tuple):
        m, n = boundary
    else:
        m = n = boundary
    if m < 0:
        raise ValueError("Number of bottom boundary strands cannot be negative")
    if n < 0:
        raise ValueError("Number of top boundary strands cannot be negative")
    return (m, n)


class Tangle():
    def __init__(self, boundary=2, crossings=None, entry_points=None, label=None):
        """
        A tangle is a fragment of a Link with some number of boundary
        strands. Tangles can be composed in various ways along their boundary strands,
        including the horizontal and vertical compositions of the tangle category.

        Inputs:

        * When boundary is an integer, then the tangle has n strands coming into both
          the top and the bottom of the tangle. When boundary is a pair of integers
          (m, n), then the tangle has m strands coming into the bottom and n coming
          into the top.

          The strands are numbered 0 to m-1 on the bottom and m to m+n-1 on the
          top, both from left to right.

        * crossings is a list of Crossing or Strand objects that comprise the tangle.
        * entry_points is a list of pairs (c, i) where c is a Crossing or Strand
          and i indexes into c.adjacent. These pairs describe the boundary strands
          in order of the strand numbering.
        * label is an arbitrary label for the tangle for informational purposes, which
          appears in the ``repr`` form of the tangle.

        Usually tangles should not be created directly using this constructor since the
        tangle operations and various primitive tangles are sufficient to create any tangle.
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
        """Put self to left of other and fuse the top-right strand of self to the top-left
        strand of other and the bottom-right strand of self to the bottom-left strand of other.

        >>> (IdentityBraid(2) + BraidTangle([1])).describe()
        'Tangle[{1,2}, {3,4}, P[1,3], X[2,4,5,5]]'
        """
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
        """Join with self *above* other, as with braid multiplication.
        (See doc.pdf)

        >>> BraidTangle([1,1]).describe()
        'Tangle[{1,2}, {3,4}, X[5,4,3,6], X[2,5,6,1]]'
        >>> (BraidTangle([1])*BraidTangle([1])).describe()
        'Tangle[{1,2}, {3,4}, X[5,4,3,6], X[2,5,6,1]]'
        """
        A, B = self.copy(), other.copy()
        (mA, nA), (mB, nB) = A.boundary, B.boundary
        if mA != nB:
            raise ValueError("Tangles must have a compatible number of strands to multiply them")

        a, b = A.adjacent, B.adjacent
        for i in range(mA):
            join_strands(a[i], b[mB + i])
        return Tangle((mB, nA), A.crossings + B.crossings, b[:mB] + a[mA:])

    def __neg__(self):
        """Mirror image of self.

        >>> (-BraidTangle([1])).describe()
        'Tangle[{1,2}, {3,4}, X[4,3,1,2]]'"""
        T = self.copy()
        for c in T.crossings:
            if not isinstance(c, Strand):
                c.rotate_by_90()
        return T

    def __or__(self, other):
        """Put self to left of other. This is like tangle addition but without the fusing of strands.

        >>> (IdentityBraid(1) | CupTangle()).describe()
        'Tangle[{1}, {2,3,4}, P[1,2], P[3,4]]'
        """
        A, B = self.copy(), other.copy()
        (mA, nA), (mB, nB) = A.boundary, B.boundary
        a, b = A.adjacent, B.adjacent
        entry_points = a[:mA] + b[:mB] + a[mA:] + b[mB:]
        return Tangle((mA + mB, nA + nB), A.crossings + B.crossings, entry_points)

    def copy(self):
        return pickle.loads(pickle.dumps(self))

    def rotate(self, s):
        """Rotate anticlockwise by s*90 degrees. This is only for (2,2) tangles.

        See ``Tangle.reshape()`` for a generalization to all tangle shapes."""
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
        """Rotate anticlockwise by 90 and take the mirror image. This is only for (2,2) tangles."""
        if self.boundary != (2, 2):
            raise ValueError("Only boundary=(2,2) tangles can be inverted")
        return -self.rotate(1)

    def numerator_closure(self):
        """The bridge closure, where consecutive pairs of strands at both the top and
        at the bottom are respectively joined by caps and cups. The numbers of
        strands at both the top and the bottom must be even. Returns a Link.

        A synonym for this is ``Tangle.bridge_closure()``.

        sage: BraidTangle([2,-1,2],4).numerator_closure().alexander_polynomial()
        t^2 - t + 1
        sage: BraidTangle([1,1,1]).rotate(1).numerator_closure().alexander_polynomial()
        t^2 - t + 1
        """
        m, n = self.boundary
        if m % 2 or n % 2:
            raise ValueError("To do bridge closure, both the top and bottom must have an even number of strands")
        T = self.copy()
        for i in range(0, m, 2):
            join_strands(T.adjacent[i], T.adjacent[i + 1])
        for i in range(0, n, 2):
            join_strands(T.adjacent[m + i], T.adjacent[m + i + 1])
        return Link(T.crossings, check_planarity=False)

    def denominator_closure(self):
        """The braid closure, where corresponding strands between the top and bottom
        are joined. The number of strands at the top must equal the number of strands at
        the bottom. Returns a Link.

        A synonym for this is ``Tangle.braid_closure()``.

        sage: BraidTangle([1,1,1]).braid_closure().alexander_polynomial()
        t^2 - t + 1
        sage: BraidTangle([1,-2,1,-2]).braid_closure().alexander_polynomial()
        t^2 - 3*t + 1
        >>> BraidTangle([1,-2,1,-2]).braid_closure().exterior().identify() # doctest: +SNAPPY
        [m004(0,0), 4_1(0,0), K2_1(0,0), K4a1(0,0), otet02_00001(0,0)]
        """
        m, n = self.boundary
        if m != n:
            raise ValueError("To do braid closure, both the top and bottom numbers of strands must be equal")
        T = self.copy()
        for i in range(n):
            join_strands(T.adjacent[i], T.adjacent[m + i])
        return Link(T.crossings, check_planarity=False)

    def link(self):
        """If its boundary is (0, 0), return this Tangle as a Link."""
        if self.boundary != (0, 0):
            raise ValueError("The boundary must be (0, 0)")
        return Link(self.copy().crossings, check_planarity=False)

    def reshape(self, boundary, displace=0):
        """Renumber the boundary strands so that the tangle has the new boundary
        shape. This is performed by either repeatedly moving the last strands from the
        bottom right to the top right or vice versa. Simultaneously, displace controls
        a rotation of the tangle where the tangle is rotated clockwise by ``displace`` steps
        (so, for example, if 0 <= displace < m then the strand numbered ``displace``
        becomes the new lower-left strand).

        This is a generalization of ``Tangle.rotate()``.
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

        This generalizes ``Tangle.rotate()``, and it is a mild specialization of ``Tangle.reshape()``.
        """
        return self.reshape(self.boundary, n)

    def circular_sum(self, other, n=0):
        """
        Glue two tangles together to form a link by gluing them vertically and then taking
        the braid closure (the ``Tangle.denominator_closure()``).
        The second tangle is rotated clockwise by n strands using ``Tangle.circular_rotate()``.
        """
        Am, An = self.boundary
        Bm, Bn = self.boundary
        if (Am, An) != (Bn, Bm):
            raise ValueError("Tangles must have compatible boundary shapes")
        return (self * (other.circular_rotate(n))).denominator_closure()

    def isosig(self, root=None, over_or_under=False):
        """
        Return a bunch of data which encodes the planar isotopy class of the
        tangle.  Of course, this is just up to isotopy of the plane
        (no Reidemeister moves).  A root can be specified with a CrossingStrand
        and ``over_or_under`` toggles whether only the underlying
        shadow (4-valent planar map) is considered or the tangle with the
        over/under data at each crossing.

        >>> BraidTangle([1]).isosig() == BraidTangle([1]).circular_rotate(1).isosig()
        True
        >>> BraidTangle([1]).isosig() == BraidTangle([-1]).isosig()
        True
        """
        copy = self.copy()
        copy._fuse_strands()
        return planar_isotopy.min_isosig(copy, root, over_or_under)

    def is_planar_isotopic(self, other, root=None, over_or_under=False) -> bool:
        return self.isosig() == other.isosig()

    def _fuse_strands(self, preserve_boundary=False, preserve_components=False):
        """Fuse all strands and delete them, even ones incident to only the boundary (unless
        ``preserve_boundary`` is True). This will eliminate Strands that are loops as well.

        If ``preserve_components`` is True, then do not fuse strands that have the
        ``component_idx`` attribute."""
        for s in reversed(self.crossings):
            if isinstance(s, Strand):
                # check that the strand is not only incident to the boundary
                if preserve_boundary and all(a[0] == self for a in s.adjacent):
                    continue
                if preserve_components and s.component_idx is not None:
                    continue
                s.fuse()
                self.crossings.remove(s)

    def __repr__(self):
        return "<Tangle: %s>" % self.label

    def describe(self, fuse_strands=True):
        """Give a PD-like description of the tangle in the form
        Tangle[{lower arcs}, {upper arcs}, P and X codes].

        If fuse_strands is True, then fuse all internal Strand nodes first.

        >>> BraidTangle([1]).describe()
        'Tangle[{1,2}, {3,4}, X[2,4,3,1]]'
        """
        T = self.copy()
        if fuse_strands:
            T._fuse_strands(preserve_boundary=True, preserve_components=True)
        T.label = 0
        # give each crossing/strand a unique identifier, which
        # is used for calculating ids for arcs
        for i, c in enumerate(T.crossings):
            c.label = i + 1
        arc_ids = {}

        def arc_key(c, i):
            """For the given entity c and index into c.adjacent,
            create a name for the incident arc. This gives something
            that's suitable for use as a dictionary key."""
            d, j = c.adjacent[i]
            return tuple(sorted([(c.label, i), (d.label, j)]))

        def arc_id(c, i):
            """Get the unique integer id associated to the arc, generating
            a fresh one if needed."""
            return arc_ids.setdefault(arc_key(c, i), len(arc_ids) + 1)
        m, n = T.boundary
        lower = "{" + ",".join(str(arc_id(T, i)) for i in range(m)) + "}"
        upper = "{" + ",".join(str(arc_id(T, i)) for i in range(m, m + n)) + "}"
        parts = []
        for c in T.crossings:
            arcs = [arc_id(c, i) for i in range(len(c.adjacent))]
            if isinstance(c, Crossing):
                parts.append("X[%s,%s,%s,%s]" % tuple(arcs))
            elif isinstance(c, Strand):
                if c.component_idx is not None:
                    parts.append(f"P[{arcs[0]},{arcs[1]}, component->{c.component_idx}]")
                else:
                    parts.append(f"P[{arcs[0]},{arcs[1]}]")
            else:
                raise TypeError("Unexpected entity")
        return f"Tangle[{lower}, {upper}{''.join(', ' + p for p in parts)}]"


Tangle.bridge_closure = Tangle.numerator_closure
Tangle.braid_closure = Tangle.denominator_closure


def ComponentTangle(component_idx):
    """The unknotted (1,1) tangle with a specified component index.
    The component index can be a negative number following the usual
    Python list indexing rules, so -1 means the component containing
    this tangle should be the last component when it is turned into
    a Link.

    >>> T=(RationalTangle(2,3)+IdentityBraid(1))|(RationalTangle(2,5)+ComponentTangle(-1))
    >>> T.describe()
    'Tangle[{1,2}, {3,4}, X[5,6,7,3], X[8,7,6,9], X[1,8,9,5], X[10,11,12,4], X[13,14,11,10], X[15,16,14,17], X[2,15,17,13], P[16,12, component->-1]]'

    >>> M=T.braid_closure().exterior() # doctest: +SNAPPY
    >>> M.dehn_fill([(1,0),(0,0)]) # doctest: +SNAPPY
    >>> M.filled_triangulation().identify() # doctest: +SNAPPY
    [m004(0,0), 4_1(0,0), K2_1(0,0), K4a1(0,0), otet02_00001(0,0)]

    >>> T=(RationalTangle(2,3)+IdentityBraid(1))|(RationalTangle(2,5)+ComponentTangle(0))
    >>> T.describe()
    'Tangle[{1,2}, {3,4}, X[5,6,7,3], X[8,7,6,9], X[1,8,9,5], X[10,11,12,4], X[13,14,11,10], X[15,16,14,17], X[2,15,17,13], P[16,12, component->0]]'

    >>> M=T.braid_closure().exterior() # doctest: +SNAPPY
    >>> M.dehn_fill([(0,0),(1,0)]) # doctest: +SNAPPY
    >>> M.filled_triangulation().identify() # doctest: +SNAPPY
    [m004(0,0), 4_1(0,0), K2_1(0,0), K4a1(0,0), otet02_00001(0,0)]

    >>> T=(RationalTangle(2,3)+ComponentTangle(0))|(RationalTangle(2,5)+ComponentTangle(0))
    >>> T.braid_closure()
    Traceback (most recent call last):
        ...
    ValueError: Two Strand objects in different components have the same component_idx values

    """
    s = Strand(component_idx=component_idx)
    return Tangle((1, 1), [s], [(s, 0), (s, 1)])


def CapTangle():
    """The unknotted (2,0) tangle."""
    cap = Strand("cap")
    return Tangle((2, 0), [cap], [(cap, 0), (cap, 1)])


def CupTangle():
    """The unknotted (0,2) tangle."""
    cup = Strand("cup")
    return Tangle((0, 2), [cup], [(cup, 0), (cup, 1)])


def ZeroTangle():
    """The zero tangle, equivalent to ``RationalTangle(0)`` or
    ``CupTangle() * CapTangle()``."""
    bot, top = Strand('B'), Strand('T')
    return Tangle(2, [bot, top],
                  [(bot, 0), (bot, 1), (top, 0), (top, 1)],
                  "ZeroTangle")


def InfinityTangle():
    """The infinity tangle, equivalent to ``RationalTangle(1, 0)`` or
    ``IdentityBraid(2)``."""
    left, right = Strand('L'), Strand('R')
    return Tangle(2, [left, right],
                  [(left, 0), (right, 0), (left, 1), (right, 1)],
                  "InfinityTangle")


def MinusOneTangle():
    """The minus one tangle, equivalent to ``RationalTangle(-1)``."""
    c = Crossing('-one')
    return Tangle(2, [c], [(c, 3), (c, 0), (c, 2), (c, 1)],
                  "MinusOneTangle")


def OneTangle():
    """The one tangle, equivalent to ``RationalTangle(1)``."""
    c = Crossing('one')
    return Tangle(2, [c], [(c, 0), (c, 1), (c, 3), (c, 2)],
                  "OneTangle")


def IntegerTangle(n):
    """The tangle equivalent to ``RationalTangle(n)``. It is
    ``n`` copies of the ``OneTangle`` joined by ``+`` when ``n`` is
    positive, and otherwise ``-n`` copies of ``MinusOneTangle``."""
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

    >>> continued_fraction_expansion(3141,1000)
    [3, 7, 10, 1, 5, 2]
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
    """
    A rational tangle. ``RationalTangle(a, b)`` gives the a/b rational tangle when ``a``
    and ``b`` are integers. If ``q`` is a rational, then ``RationalTangle(q)`` gives the
    corresponding rational tangle.

    This is a class that extends Tangle since it provides some additional information as
    attributes: ``fraction`` gives (a, b) and ``partial_quotients`` gives the continued
    fraction expansion of ``abs(a)/b``.

    >>> RationalTangle(2,5).braid_closure().exterior().identify() # doctest: +SNAPPY
    [m004(0,0), 4_1(0,0), K2_1(0,0), K4a1(0,0), otet02_00001(0,0)]
    """
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
    The braid with n strands and no crossings.

    >>> IdentityBraid(0).describe()
    'Tangle[{}, {}]'
    >>> IdentityBraid(1).describe()
    'Tangle[{1}, {2}, P[1,2]]'
    >>> IdentityBraid(2).describe()
    'Tangle[{1,2}, {3,4}, P[1,3], P[2,4]]'
    >>> IdentityBraid(-1)
    Traceback (most recent call last):
        ...
    ValueError: Expecting non-negative int
    """
    if n < 0:
        raise ValueError("Expecting non-negative int")
    strands = [Strand() for i in range(n)]
    entry_points = [(s, 0) for s in strands] + [(s, 1) for s in strands]
    return Tangle(n, strands, entry_points,
                  f"IdentityBraid({n})")


def BraidTangle(gens, n=None):
    """
    Create an (n,n) tangle from a braid word.

    Input:

    * gens is a list of nonzero integers, positive for the positive generator
      and negative for the negative generator
    * n is the number of strands. By default it is inferred to be the least
      number of strands that works for the given list of generators

    >>> BraidTangle([], 1)
    <Tangle: IdentityBraid(1)>
    >>> BraidTangle([1]).describe()
    'Tangle[{1,2}, {3,4}, X[2,4,3,1]]'
    >>> BraidTangle([-1]).describe()
    'Tangle[{1,2}, {3,4}, X[1,2,4,3]]'
    >>> BraidTangle([1],3).describe()
    'Tangle[{1,2,3}, {4,5,6}, P[3,6], X[2,5,4,1]]'
    >>> BraidTangle([2],3).describe()
    'Tangle[{1,2,3}, {4,5,6}, P[1,4], X[3,6,5,2]]'
    >>> BraidTangle([1,2]).describe()
    'Tangle[{1,2,3}, {4,5,6}, X[7,5,4,1], X[3,6,7,2]]'
    >>> BraidTangle([1,2,1]).describe()
    'Tangle[{1,2,3}, {4,5,6}, X[7,5,4,8], X[3,6,7,9], X[2,9,8,1]]'
    """
    if n is None:
        n = max(-min(gens), max(gens)) + 1
    # n is the braid index

    def gen(i):
        g = OneTangle() if i < 0 else MinusOneTangle()
        return IdentityBraid(abs(i) - 1) | g | IdentityBraid(n - abs(i) - 1)

    b = IdentityBraid(n)
    for i in gens:
        if i == 0:
            raise ValueError("Generators must be nonzero integers")
        if abs(i) >= n:
            raise ValueError("Generators must have magnitude less than n")
        b = b * gen(i)
    return b
