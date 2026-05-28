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

from collections import OrderedDict, Counter
from .ordered_set import OrderedSet
from .links import Crossing, Strand, Link, CrossingStrand, CrossingEntryPoint

class CyclicList(list):
    def __init__(self, iterable):
        super().__init__(iterable)

    def __getitem__(self, i):
        if isinstance(i, slice):
            return super().__getitem__(i)
        return super().__getitem__(i % len(self))

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

class ArcLabels(OrderedDict):
    def __init__(self, iterable = []):
        super().__init__(iterable)

        self.counter = 0

    def add(self, c, advance):
        if c not in self:
            self[c] = self.counter
            if advance:
                self.counter += 1
        else:
            raise ValueError(f"Each CEP should only be labeled once, but {c} is already labeled with {self[c]}")

class TangleComponents(list):
    def add(self, c):
        component = c.component()
        self.append(component)
        return component

class Tangle:
    def __init__(self, boundary=2, crossings=None, entry_points=None, build = True, label=None, start_orientations = None, component_starts = None):
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

        Tangles now support creation from PD_code, for example:
        
        >>> Tangle(3, [[0,4,1,5],[1,8,2,9],[2,7,3,6],[5,9,6,10]], [0,4,8,10,3,7], label = 'RIII')
        <Tangle: RIII: 3 comp; 4 cross; (3, 3) boundary>

        see doc of ``PD_code`` for more details.
        """
        if label is None:
            self.label = id(self)
        else:
            self.label = label

        m, n = decode_boundary(boundary)        
        component_starts = component_starts
        start_orientations = start_orientations
        self.strand_labels = CyclicList(m * [None] + n * [None])
        self.strand_components = CyclicList(m * [None] + n * [None])

        if crossings is None:
            crossings = []
        else:
            if isinstance(crossings, str) or isinstance(entry_points, str):
                raise NotImplementedError("Not Implemented. If you are trying to create a tangle from a PD code, input the PD code as a list instead.")
            
            if (len(crossings) > 0 and not isinstance(crossings[0], (Strand, Crossing)))\
                  or (entry_points is not None and len(entry_points) > 0 and not isinstance(entry_points[0], (CrossingStrand, list, tuple))):
                assert component_starts is None and start_orientations is None, "Specifying components_starts and start_orientations is not compatible with creating from PD codes"

                crossings, component_starts, entry_points = self._crossings_from_PD_code(crossings, entry_points)
                start_orientations = component_starts[:]

        # the pair for the number of lower strands and the number of upper strands
        self.boundary = (m, n)
        self.boundary_strands = CyclicList([])
        # -1 if entering Tangle, 1 if exiting Tangle, 0 if not yet oriented. 
        self.boundary_signs = CyclicList(m * [0] + n * [0])

        # a list of (c, i) pairs for the boundary strands. Each c will reciprocally
        # contain (self, j) where j is the strand number. The fact this is called
        # 'adjacent' means that the Tangle can take part in the joining protocol
        # implemented in join_strands.
        self.adjacent = CyclicList((m + n) * [None])
        entry_points = entry_points or []
        if len(entry_points) != m + n:
            raise ValueError("The number of boundary strands is not equal to the length"
                             " of entry_points")

        for i, e in enumerate(entry_points):
            this_strand = Strand(label = f'TSE({str(self)}, {i})')
            self.boundary_strands.append(this_strand)
            join_strands(e, (this_strand, 1))
            join_strands((self, i), (this_strand, 0))
            
        if not all(isinstance(c, (Crossing, Strand)) for c in crossings):
            raise ValueError("Every element of crossings must be a Crossing or a Strand")

        self.unlinked_unknot_components = 0
        # Note that crossings in Tangle can contain Strands for now
        # which will be removed after build
        self.crossings = crossings

        if build:
            self._build(start_orientations, component_starts)
            assert self.is_oriented(), 'Tangle is not oriented after build'

            # Remove all Strands from crossings and components. 
            # Note that this will not affect strands in boundary_strands
            for s in reversed(crossings):
                if isinstance(s, Strand):
                    comp = self.components[s.strand_component]

                    if isinstance(comp[0].crossing, Tangle):
                        for cep in reversed(comp):
                            if cep.crossing == s:
                                comp.remove(cep)
                                break
                        else:
                            raise RuntimeError(f"Component strand {s} not found in component {comp}")
                        
                        # Note that the components are always built following the orientation
                        # hence below always insists that the comp_id is labeled on the entrance strand
                        if s.component_idx is not None:
                            comp_id =  s.component_idx
                            if comp[1].crossing.component_idx is not None:
                                assert comp[1].crossing.component_idx == comp_id
                            else:
                                comp[1].crossing.component_idx = comp_id
                    if s.is_loop():
                        self.unlinked_unknot_components += 1
                    else:
                        s.fuse()
                    self.crossings.remove(s)

    def __getitem__(self, i):
        return (self, i % (self.boundary[0] + self.boundary[1]))

    def __setitem__(self, i, other):
        o, j = other
        self.adjacent[i % (self.boundary[0] + self.boundary[1])] = other
        o.adjacent[j] = (self, i)

    def is_upward(self):
        return self.boundary_signs == CyclicList([-1] * self.boundary[0] + [1] * self.boundary[1])
    
    def is_downward(self):
        return self.boundary_signs == CyclicList([1] * self.boundary[0] + [-1] * self.boundary[1])

    def is_oriented(self):
        return all(s != 0 for s in self.boundary_signs)
    
    def make_upward(self):
        """
        Change the orientation of the tangle, trying to make it upwardly oriented. 
        The order of components is preserved. 
        
        >>> T = BraidTangle([1,2,1])
        >>> T.reverse_orientation([1,2])
        >>> T.is_upward()
        False
        >>> T.make_upward()
        >>> T.is_upward()
        True

        Like alternating() for Links, this may fail silently if there is no orientation
        which makes the tangle upwardly oriented.

        >>> T = RationalTangle(-2,3)
        >>> T.make_upward()
        >>> T.is_upward()
        False
        """
        if self.is_upward():
            return
        
        assert self.is_oriented(), 'Tangle should be oriented to tell if it is upward'
        
        to_reverse = set()
        for i in range(self.boundary[0]):
            if self.boundary_signs[i] == 1:
                to_reverse.add(self.boundary_strands[i].strand_component)
        
        self.reverse_orientation(to_reverse)

    def entry_points(self):
        assert self.is_oriented(), 'Tangle should be oriented to tell the entry points'
        return [CrossingEntryPoint(self, i) for i in range(self.boundary[0] + self.boundary[1]) 
                if self.boundary_signs[i] == -1]
    
    def update_label(self, label):
        self.label = label
        for i, s in enumerate(self.boundary_strands):
            s.label = f'TSE({self}, {i})'

    def _build(self, start_orientations=None, component_starts=None):
        self._orient_crossings(start_orientations=start_orientations)
        self._build_components(component_starts=component_starts)

    def _clear(self):
        self.components = None
        for c in self.crossings:
            c._clear()
        for s in self.boundary_strands:
            s._clear()
        
        self.boundary_signs = CyclicList(self.boundary[0] * [0] + self.boundary[1] * [0])
        self.strand_labels = CyclicList(self.boundary[0] * [None] + self.boundary[1] * [None])
        self.strand_components = CyclicList(self.boundary[0] * [None] + self.boundary[1] * [None])

    def _rebuild(self, same_components_and_orientations = False):
        if same_components_and_orientations:
            # Hopefully we have enough of the original components left
            # to figure out what this is.  Otherwise, new choices will
            # be made as in the default algorithm.
            start_css = self._start_orientations()
            
        self._clear()
        if same_components_and_orientations:
            self._build(start_orientations=start_css,
                        component_starts=start_css)
        else:
            self._build()

    def all_crossings_oriented(self):
        return all(c.sign != 0 for c in self.crossings + self.boundary_strands)

    def _orient_crossings(self, start_orientations=None):
        if self.all_crossings_oriented():
            return
        if start_orientations is None:
            start_orientations = list()
        else: # copy as algorithm modifies this list
            start_orientations = list(start_orientations)

        remaining = OrderedSet(
            [(c, i) for c in self.crossings + list(reversed(self.boundary_strands)) 
             for i in range(c._adjacent_len) if c.sign == 0])
            
        while len(remaining):
            if len(start_orientations) > 0:
                c, i = start = start_orientations.pop()
            else:
                c, i = start = remaining.pop()

            is_reversed = False
            finished = False
            while not finished:
                if is_reversed:
                    c.make_tail(i)
                else:
                    c.make_head(i)
                remaining.discard((c, i))
                
                if not c.adjacent[i][0] == self:
                    d, j = c.adjacent[i]
                    remaining.discard((d, j))
                    s = d._adjacent_len // 2
                    c, i = d, (j + s) % (2 * s)

                    finished = (c, i) == start
                else:
                    boundary_index = c.adjacent[i][1]
                    assert self.boundary_signs[boundary_index] == 0, f'Boundary {boundary_index} is unexpectedly signed, something is wrong with the gluings'
                    
                    if is_reversed:
                        # Hit the boundary of the tangle from both sides, 
                        # done with this component
                        self.boundary_signs[boundary_index] = -1
                        finished = True
                    else:
                        # Hit the boundary of the tangle, 
                        # now go back and orient reversely
                        self.boundary_signs[boundary_index] = 1
                        is_reversed = True

                        c, i = start
                        s = c._adjacent_len // 2
                        c, i = c, (i + s) % (2 * s)

        for c in self.crossings:
            c.orient()

        for s in self.boundary_strands:
            s.orient()

    def _build_components(self, component_starts=None):
        """
        Each component is stored as a list of *entry points* to crossings. 
        If the component starts and ends at the boundary of tangles, 
        the corresponding CrossingEntryPoint(self, boundary_index) 
        will be put at the tail and the head of the list. 

        If provided, the component_starts must consist of one
        CrossingEntryPoint per component.

        >>> len(RationalTangle(-2, 3).components)
        2
        >>> len(Tangle(3, [[0,4,1,5],[1,8,2,9],[2,7,3,6],[5,9,6,10]], 
        ... [0,4,8,10,3,7], label = 'RIII').components)
        3
        >>> len(((RationalTangle(2,3)+IdentityBraid(1))|(RationalTangle(2,5)+ComponentTangle(-1))).components)
        2
        """
        if component_starts is not None:
            # Take all CrossingStrand and CrossingEntryPoint objects
            # and turn them into CrossingEntryPoints
            component_starts = [cs.crossing.entry_points()[cs.strand_index % 2 if isinstance(cs.crossing, Crossing) else 0]
                                for cs in component_starts]
        remaining, components = OrderedSet(self.crossing_entries()), TangleComponents()
        other_crossing_entries = []
        self.labels = labels = ArcLabels()
        for c in self.crossings:
            c._clear_strand_info()
        for s in self.boundary_strands:
            s._clear_strand_info()

        while len(remaining):
            if component_starts:
                d = component_starts[len(components)]
            elif len(components) == 0:
                d = remaining.pop()
            else: # prioritize labeling crossing strands that are adjacent to already labeled ones
                found, comp_index = False, 0
                while not found and comp_index < len(components):
                    others = other_crossing_entries[comp_index]
                    if others:
                        for j, d in enumerate(others):
                            if d.component_label() is None:
                                if labels[d.other()] % 2 == 0:
                                    d = d.next()
                                found = True
                                break
                        other_crossing_entries[comp_index] = others[j:]
                    comp_index += 1

                if not found:
                    d = remaining.pop()

            component = components.add(d)

            # Label arcs along the component
            for i, c in enumerate(component):
                if isinstance(c.crossing, Tangle) and self.boundary_signs[c.strand_index] == 0:
                    assert len(component) > 2

                    if i == 0:
                        self.boundary_signs[c.strand_index] = -1
                    else:
                        assert i == len(component) - 1
                        self.boundary_signs[c.strand_index] = 1
                # if is a Crossing or at the end of the component, advance the label
                # otherwise, don't advance the label since we will still be on the same arc
                if isinstance(c.crossing, Crossing) or i == len(component) - 1:
                    advance = True
                else:
                    advance = False

                labels.add(c, advance)

            others = []
            for c in component:
                c.label_crossing(len(components) - 1, labels)
                o = c.other()
                if o is not None and o.component_label() is None:
                    others.append(o)
            other_crossing_entries.append(others)
            remaining.difference_update(component)

        self.components = components

    def crossing_entries(self):
        ans = []
        for C in self.crossings:
            ans += C.entry_points()

        for s in reversed(self.boundary_strands):
            ans += s.entry_points()

        return ans

    def _crossings_from_PD_code(self, code, entry_points):
        """
        entry_points as INPUT: a list of labels of arcs left open in the tangle
                     as OUTPUT: a list of CrossingStrands
        """
        assert Counter(entry_points).most_common(1)[0][1] <= 2, "Each entry point label should appear at most twice"

        labels = set()
        for X in code:
            for i in X:
                labels.add(i)
        for x in entry_points:
            labels.add(x)

        gluings = OrderedDict()

        for c, X in enumerate(code):
            for i, x in enumerate(X):
                if x in gluings:
                    gluings[x].append((c, i))
                else:
                    gluings[x] = [(c, i)]

        if any(len(v) > 2 for v in gluings.values()):
            raise ValueError("PD code isn't consistent")

        crossings = [Crossing(i) for i, d in enumerate(code)]
        
        for item in gluings.values():
            if len(item) > 1:
                (c, i), (d, j) = item
                crossings[c][i] = crossings[d][j]

        entry_strands = []
        entry_dict = dict()

        for i, x in enumerate(entry_points):
            if x in gluings:
                entry_strands.append(crossings[gluings[x][0][0]].crossing_strands()[gluings[x][0][1]])
            else:
                this_strand = Strand(label = f'PDSE({str(self)}, {i})')
                crossings.append(this_strand)
                if x not in entry_dict:
                    entry_strands.append((this_strand, 0))
                    entry_dict[x] = (this_strand, 1)
                else:
                    entry_strands.append((this_strand, 1))
                    join_strands(entry_dict[x], (this_strand, 0))

        component_starts = self._component_starts_from_PD(
            code, labels, gluings, entry_dict)

        component_starts = [crossings[c].crossing_strands()[i] 
                            if not isinstance(c, Strand) else CrossingStrand(c, i)
                            for (c, i) in component_starts]
        
        return crossings, component_starts, entry_strands
    
    def PD_code(self, KnotTheory=False, min_strand_index = 0):
        """
        The planar diagram code for the tangle. Unlike for links, it returns two extra fields,
        boundary and entry_info in addition to the PD code of crossings, in order to specify
        how the boundary and entries of the tangle is arranged. The fields are ordered as follows:

        boundary, PD, entry_info

        so that they can be unpacked immediately for creating Tangles.

        >>> RationalTangle(-1,2).PD_code()
        ((2, 2), [(1, 5, 2, 4), (3, 1, 4, 0)], [0, 3, 2, 5])
        >>> BraidTangle([1,2,1]).PD_code()
        ((3, 3), [(7, 5, 8, 4), (6, 2, 7, 1), (3, 1, 4, 0)], [0, 3, 6, 8, 5, 2])
        >>> Tangle(*RationalTangle(-1,2).PD_code()).PD_code()
        ((2, 2), [(1, 5, 2, 4), (3, 1, 4, 0)], [0, 3, 2, 5])
        """
        PD = []
        entry_info = [s + min_strand_index for s in self.strand_labels]

        for c in self.crossings:
            if isinstance(c, Crossing):
                PD.append([s + min_strand_index for s in c.strand_labels])

        if KnotTheory:
            PD = "PD" + repr(PD).replace('[', 'X[')[1:]
            entry_info = "EP" + repr(entry_info)
        else:
            PD = [tuple(x) for x in PD]

        return self.boundary, PD, entry_info

    def _component_starts_from_PD(self, code, labels, gluings, entry_dict):
        """
        A PD code determines an order and orientation on the tangle
        components as follows, where we view the code as labels on the
        strands at the point where two crossings are stuck together.

        1.  The minimum label on each component is used to order the
            components.

        2.  Each component is oriented by finding its minimal label,
            looking at the labels of its two neighbors, and then
            orienting the component towards the smaller of those two.

        This is designed so that a PLink-generated PD_code results in a
        link with the same component order and orientation.
        """
        starts = []
        while labels:
            m = min(labels)
            labels.remove(m)

            if m not in gluings:
                next_label = m
                starts.append(entry_dict[m])
            elif len(gluings[m]) == 1:
                # entrance strand of the tangle
                [(c, index)] = gluings[m]

                j = (index + 2) % 4
                next_label = code[c][j]
                direction = (c, j)

                starts.append(direction)
            else:
                (c1, index1), (c2, index2) = gluings[m]
                if c1 == c2:
                    # loop at strand, take next strand to be next smallest label
                    # on crossing
                    next_label = min(set(code[c1]) - {m})
                    direction = (c1, code[c1].index(next_label))
                    starts.append(direction)
                else:
                    # strand connects two different crossings, take next strand to
                    # be next smallest label on two 'opposite' strands
                    j1, j2 = (index1 + 2) % 4, (index2 + 2) % 4
                    l1, l2 = code[c1][j1], code[c2][j2]
                    if l1 < l2:
                        next_label = l1
                        direction = (c1, j1)
                    elif l2 < l1:
                        next_label = l2
                        direction = (c2, j2)
                    else:
                        # We have a component of length 2, so now rely on
                        # the convention that the first position at a PD
                        # crossing is a directed entry point. (If both
                        # crossings are over or both under, the
                        # orientation is arbitrary anyway.)
                        next_label = l1

                        # The strand labeled m is oriented c2 --> c1 if
                        # and only if either l1 = l2 is the incoming
                        # understrand of c2 or m is incoming understrand
                        # at c1.
                        if code[c2][0] == l1 or code[c1][0] == m:
                            direction = (c1, j1)
                        else:
                            direction = (c2, j2)

                    starts.append(direction)

            # Component start recorded. Erase the rest of the component 
            # by traversing along it and remove the labels
            while next_label != m:
                labels.remove(next_label)
                g = gluings[next_label]
                if len(g) == 1:
                    break
                other_direction = g[1 - g.index(direction)]
                direction = (other_direction[0], (other_direction[1] + 2) % 4)
                next_label = code[direction[0]][direction[1]]

        return starts


    # The following operators always clear the current orientations on both tangles 
    # and recreate an orientation with default behaviour.
    def __add__(self, other):
        """Put self to left of other and fuse the top-right strand of self to the top-left
        strand of other and the bottom-right strand of self to the bottom-left strand of other.

        >>> (IdentityBraid(2) + BraidTangle([1])).describe()
        'Tangle[{1,2}, {3,4}, X[2,4,5,5], P[1,3]]'
        """
        A, B = self.copy(), other.copy()
        (mA, nA), (mB, nB) = A.boundary, B.boundary
        if mA == 0 or mB == 0 or nA == 0 or nB == 0:
            raise ValueError("Tangles must have at least one top and bottom strand each.")
        a, b = A.adjacent, B.adjacent
        join_strands(a[mA - 1], b[0])
        join_strands(a[mA + nA - 1], b[mB])
        entry_points = a[:mA - 1] + b[1:mB] + a[mA:mA + nA - 1] + b[mB + 1:]

        crossings = A.crossings + A.boundary_strands + B.crossings + B.boundary_strands

        for c in crossings:
            c._clear()

        return Tangle((mA + mB - 2, nA + nB - 2), 
                      crossings, 
                      entry_points)

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

        crossings = A.crossings + A.boundary_strands + B.crossings + B.boundary_strands

        for c in crossings:
            c._clear()

        return Tangle((mB, nA), 
                      crossings,
                      b[:mB] + a[mA:])

    def __neg__(self):
        """Mirror image of self.

        >>> (-BraidTangle([1])).describe()
        'Tangle[{1,2}, {3,4}, X[4,3,1,2]]'"""
        T = self.copy()
        for c in T.crossings:
            if not isinstance(c, Strand):
                c.rotate_by_90()
                c.orient()
        return T

    def _start_orientations(self):
        """
        Obtain the start orientations according to the current orientation 
        and components (the latter may be outdated)
        """
        start_css = []
        for comp in self.components:
            for cs in comp:
                if cs.crossing in self.crossings + self.boundary_strands:
                    s = cs.crossing._adjacent_len // 2
                    start_css.append(cs.rotate(s))
                    break

        return start_css

    def __or__(self, other):
        """
        Put self to left of other. This is like tangle addition but without the fusing of strands.
        Preserves the orientations of both tangles, since no gluing happens.

        >>> (IdentityBraid(1) | CupTangle()).describe()
        'Tangle[{1}, {2,3,4}, P[1,2], P[3,4]]'

        >>> T = BraidTangle([1,2,1])
        >>> T.reverse_orientation([1,2])
        >>> T.is_upward()
        False
        >>> T.boundary_signs
        [-1, 1, 1, -1, -1, 1]
        >>> TT = T | snappy.RationalTangle(1,2)
        >>> TT.is_upward()
        False
        >>> TT.boundary_signs
        [-1, 1, 1, -1, -1, -1, -1, 1, 1, 1]
        """
        A, B = self.copy(), other.copy()
        (mA, nA), (mB, nB) = A.boundary, B.boundary
        a, b = A.adjacent, B.adjacent
        entry_points = a[:mA] + b[:mB] + a[mA:] + b[mB:]
        crossings = A.crossings + A.boundary_strands + B.crossings + B.boundary_strands

        start_css = A._start_orientations() + B._start_orientations()

        return Tangle((mA + mB, nA + nB), 
                      crossings, 
                      entry_points,
                      start_orientations=start_css,
                      component_starts=start_css)

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

        T.boundary_strands = [T.boundary_strands[rotate[i]] for i in range(4)]
        T._rebuild(True)

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

        crossings = T.crossings + T.boundary_strands
        for c in crossings:
            c._clear()

        return Link(crossings, check_planarity=False)

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

        crossings = T.crossings + T.boundary_strands
        for c in crossings:
            c._clear()
            
        return Link(crossings, check_planarity=False)

    def link(self):
        """If its boundary is (0, 0), return this Tangle as a Link."""
        if self.boundary != (0, 0):
            raise ValueError("The boundary must be (0, 0)")
        
        crossings = self.copy().crossings
        for c in crossings:
            c._clear()

        return Link(crossings, check_planarity=False)

    def reshape(self, boundary, displace=0):
        """Renumber the boundary strands so that the tangle has the new boundary
        shape. This is performed by either repeatedly moving the last strands from the
        bottom right to the top right or vice versa. Simultaneously, displace controls
        a rotation of the tangle where the tangle is rotated clockwise by ``displace`` steps
        (so, for example, if 0 <= displace < m then the strand numbered ``displace``
        becomes the new lower-left strand).

        This is a generalization of ``Tangle.rotate()``.

        >>> T = BraidTangle([1,2,1])
        >>> T.PD_code() 
        ((3, 3), [(7, 5, 8, 4), (6, 2, 7, 1), (3, 1, 4, 0)], [0, 3, 6, 8, 5, 2])
        >>> T.reshape((4,2)).PD_code()
        ((4, 2), [(7, 5, 8, 4), (6, 2, 7, 1), (3, 1, 4, 0)], [0, 3, 6, 2, 8, 5])
        >>> T.reshape((4,2), displace = 1).PD_code()
        ((4, 2), [(7, 5, 8, 4), (6, 2, 7, 1), (3, 1, 4, 0)], [3, 6, 2, 5, 0, 8])
        """
        m, n = self.boundary
        Tm, Tn = decode_boundary(boundary)
        if (m, n) == (Tm, Tn) and displace == 0:
            return self
        if m + n != Tm + Tn:
            raise ValueError("Reshaping requires the tangle have the same number of boundary"
                             " strands as in the new boundary.")
        anticlockwise = [i for i in range(m)] + list(reversed([m + i for i in range(n)]))
        rotate = dict(zip(anticlockwise, rotate_list(anticlockwise, displace)))
        
        T = self.copy()

        displaced_adj = [T.adjacent[rotate[i]] for i in range(m + n)]
        # The 'adjacent' array but in total counterclockwise order
        adj_ccw = displaced_adj[:m] + list(reversed(displaced_adj[m:]))
        T.adjacent = adj_ccw[:Tm] + list(reversed(adj_ccw[Tm:]))
        T.boundary = (Tm, Tn)
        for i, (o, j) in enumerate(T.adjacent):
            o.adjacent[j] = (T, i)

        displaced_bd_strands = [T.boundary_strands[rotate[i]] for i in range(m + n)]
        bd_strands_ccw = displaced_bd_strands[:m] + list(reversed(displaced_bd_strands[m:]))
        T.boundary_strands = bd_strands_ccw[:Tm] + list(reversed(bd_strands_ccw[Tm:]))

        T._rebuild(True)

        return T

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

    def _to_old_tangle(self):
        from . import old_tangles
        copy = self.copy()
        
        return old_tangles.Tangle(copy.boundary, 
                         copy.crossings + copy.boundary_strands, 
                         copy.adjacent, 
                         copy.label)

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
        >>> BraidTangle([1,1]).isosig() == BraidTangle([-1,-1]).isosig()
        True
        >>> BraidTangle([1,1]).isosig(over_or_under=True) == BraidTangle([-1,-1]).isosig(over_or_under=True)
        False
        """

        return self._to_old_tangle().isosig(root = root,
                                           over_or_under=over_or_under)

    def reverse_orientation(self, component_index):
        """
        Reverse the orientation of components specified by component_index,
        changing the current tangle and the signs of crossings. 

        component_index: either a single index of component or a list of indices of components

        >>> T = BraidTangle([1,2,1])
        >>> T
        <Tangle: BraidTangle([1, 2, 1], 3): 3 comp; 3 cross; (3, 3) boundary>
        >>> T.PD_code()
        ((3, 3), [(7, 5, 8, 4), (6, 2, 7, 1), (3, 1, 4, 0)], [0, 3, 6, 8, 5, 2])
        >>> T.reverse_orientation(1)
        >>> T.PD_code()
        ((3, 3), [(7, 3, 8, 4), (6, 2, 7, 1), (4, 0, 5, 1)], [0, 5, 6, 8, 3, 2])
        >>> T.reverse_orientation([1,2])
        >>> T.PD_code()
        ((3, 3), [(6, 4, 7, 5), (7, 1, 8, 2), (3, 1, 4, 0)], [0, 3, 8, 6, 5, 2])
        >>> T.reverse_orientation([0,2])
        >>> T.PD_code()
        ((3, 3), [(7, 5, 8, 4), (6, 0, 7, 1), (3, 1, 4, 2)], [2, 3, 6, 8, 5, 0])
        >>> T.reverse_orientation([0])
        >>> T.PD_code()
        ((3, 3), [(7, 5, 8, 4), (6, 2, 7, 1), (3, 1, 4, 0)], [0, 3, 6, 8, 5, 2])
        """
        if not isinstance(component_index, (set, list, tuple)):
            component_index = [component_index]

        org_entries = []
        for comp in self.components:
            for cs in comp:
                if cs.crossing in self.crossings + self.boundary_strands:
                    org_entries.append(cs)
                    break

        new_starts = []
        for i, cs in enumerate(org_entries):
            if i not in component_index:
                c, e = cs.crossing, cs.strand_index
                s = c._adjacent_len // 2
                reversed_cs = CrossingStrand(c, (e + s) % (2 * s))
                new_starts.append(reversed_cs)
            else:
                new_starts.append(cs)

        self._clear()
        self._build(start_orientations = new_starts,
                    component_starts = new_starts)

    def faces(self):
        """
        The faces are the complementary regions of the tangle diagram in the disk,
        where the boundary of the disk is thought of as the cusp.

        Each face is given as a list of corners of crossings as one
        goes around *clockwise*.  These corners are recorded as
        CrossingStrands, where CrossingStrand(c, j) denotes the corner
        of the face abutting crossing c between strand j and j + 1;
        similarly, if c is the tangle itself, it denots the corner 
        as one stands at the j-th boundary entry and look *counterclockwisely*.

        Alternatively, the sequence of CrossingStrands can be regarded
        as the *heads* of the oriented edges of the face.    

        >>> len(IdentityBraid(2).faces())
        3
        >>> len(BraidTangle([1,2,1]).faces())
        7
        """
        corners = OrderedSet([CrossingStrand(c, i)
                              for c in self.crossings + self.boundary_strands for i in range(c._adjacent_len)])
        faces = []
        while len(corners):
            cs0 = corners.pop()
            face = [cs0]
            next = cs0
            while True:
                # Next two lines equiv to: next = next.next_corner()
                c, e = next.crossing, next.strand_index
                if isinstance(c, Tangle):
                    if e == 0:
                        next = CrossingStrand(*c.adjacent[c.boundary[0]])
                    elif e < c.boundary[0]:
                        next = CrossingStrand(*c.adjacent[e-1])
                    elif e < c.boundary[0] + c.boundary[1] - 1:
                        next = CrossingStrand(*c.adjacent[e+1])
                    else:
                        assert e == c.boundary[0] + c.boundary[1] - 1 
                        next = CrossingStrand(*c.adjacent[c.boundary[0]-1])
                else:
                    next = next.next_corner()

                if next == cs0:
                    faces.append(face)
                    break
                else:
                    corners.discard(next)
                    face.append(next)

        return faces

    def simplify(self, mode = 'basic', type_III_limit = 100):
        """
        Tries to simplify the tangle diagram. Returns whether it succeeded 
        in reducing the number of crossings. Modifies the tangle in place,
        and unknot components which are also unlinked may be silently discarded. 
        The ordering of ``components`` is not always preserved.

        The following strategies can be employed.

        1. In the default ``basic`` mode, it does Reidemeister I and II moves
           until none are possible.

        2. In ``level`` mode, it does random Reidemeister III moves, reducing
           the number of crossings via type I and II moves whenever possible.
           The process stops when it has done ``type_III_limit`` *consecutive*
           type III moves without any simplification.

        The ``pickup`` and ``global`` modes are currently not available for tangles.

        Some examples:

        >>> T = Tangle(2, [[0,3,1,4],[1,5,2,4]], [0,3,2,5], label = 'RII')
        >>> T
        <Tangle: RII: 2 comp; 2 cross; (2, 2) boundary>
        >>> T.simplify('basic')
        True
        >>> T
        <Tangle: RII: 2 comp; 0 cross; (2, 2) boundary>
        >>> T.simplify('basic') # Already done all it can
        False
        
        >>> T = Tangle(3, [[0,4,1,5],[1,8,2,9],[2,7,3,6],[5,9,6,10]], 
        ... [0,4,8,10,3,7], label = 'RIII')
        >>> T
        <Tangle: RIII: 3 comp; 4 cross; (3, 3) boundary>
        >>> T.simplify('basic')
        False
        >>> T # No change happens
        <Tangle: RIII: 3 comp; 4 cross; (3, 3) boundary>
        >>> T.simplify('level')
        True
        >>> T
        <Tangle: RIII: 3 comp; 2 cross; (3, 3) boundary>
        """
        from . import simplify
        if mode == 'basic':
            return simplify.basic_simplify(self)
        elif mode == 'level':
            return simplify.simplify_via_level_type_III(self, type_III_limit)
        else:
            raise NotImplementedError()

    def is_planar_isotopic(self, other, root=None, over_or_under=False) -> bool:
        return self.isosig(root = root, over_or_under=over_or_under) == other.isosig(root = root, over_or_under = over_or_under)

    def __repr__(self):
        return "<Tangle: %s: %d comp; %d cross; (%d, %d) boundary>" % (self.label, len(self.components), len(self.crossings), self.boundary[0], self.boundary[1])
    
    def __str__(self):
        return "<Tangle: %s>" % self.label

    def describe(self, fuse_strands=True):
        """Give a PD-like description of the tangle in the form
        Tangle[{lower arcs}, {upper arcs}, P and X codes].

        If fuse_strands is True, then fuse all internal Strand nodes first.

        >>> BraidTangle([1]).describe()
        'Tangle[{1,2}, {3,4}, X[2,4,3,1]]'
        """

        return self._to_old_tangle().describe(fuse_strands=fuse_strands)


Tangle.bridge_closure = Tangle.numerator_closure

Tangle.braid_closure = Tangle.denominator_closure


def ComponentTangle(component_idx):
    """The unknotted (1,1) tangle with a specified component index.
    The component index can be a negative number following the usual
    Python list indexing rules, so -1 means the component containing
    this tangle should be the last component when it is turned into
    a Link.

    >>> ComponentTangle(2)
    <Tangle: ComponentTangle(2): 1 comp; 0 cross; (1, 1) boundary>

    >>> T=(RationalTangle(2,3)+IdentityBraid(1))|(RationalTangle(2,5)+ComponentTangle(-1))
    >>> T.describe()
    'Tangle[{1,2}, {3,4}, X[5,3,6,7], X[8,5,7,9], X[1,8,9,6], X[10,4,11,12], X[13,14,12,11], X[14,15,16,10], X[17,16,15,13], P[2,17, component->-1]]'

    >>> M=T.braid_closure().exterior() # doctest: +SNAPPY
    >>> M.dehn_fill([(1,0),(0,0)]) # doctest: +SNAPPY
    >>> M.filled_triangulation().identify() # doctest: +SNAPPY
    [m004(0,0), 4_1(0,0), K2_1(0,0), K4a1(0,0), otet02_00001(0,0)]

    >>> T=(RationalTangle(2,3)+IdentityBraid(1))|(RationalTangle(2,5)+ComponentTangle(0))
    >>> T.describe()
    'Tangle[{1,2}, {3,4}, X[5,3,6,7], X[8,5,7,9], X[1,8,9,6], X[10,4,11,12], X[13,14,12,11], X[14,15,16,10], X[17,16,15,13], P[2,17, component->0]]'

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
    return Tangle((1, 1), [s], [(s, 0), (s, 1)], label = f'ComponentTangle({component_idx})')


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
        T.update_label(f"IntegerTangle({n})")
        return T
    elif n < 0:
        T = -IntegerTangle(-n)
        T.update_label(f"IntegerTangle({n})")
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

    >>> RationalTangle(-2,3)
    <Tangle: RationalTangle(-2, 3): 2 comp; 3 cross; (2, 2) boundary>

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

        crossings = T.crossings + T.boundary_strands

        for c in crossings:
            c._clear()

        Tangle.__init__(self, 2, 
                        crossings,
                        T.adjacent,
                        label = f"RationalTangle({a}, {b})")

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
    >>> IdentityBraid(5)
    <Tangle: IdentityBraid(5): 5 comp; 0 cross; (5, 5) boundary>
    >>> IdentityBraid(-1)
    Traceback (most recent call last):
        ...
    ValueError: Expecting non-negative int
    """
    if n < 0:
        raise ValueError("Expecting non-negative int")
    entry_points = 2* [i for i in range(n)] 
    return Tangle(n, [], entry_points,
                  label = f"IdentityBraid({n})")


def BraidTangle(gens, n=None):
    """
    Create an (n,n) tangle from a braid word.

    Input:

    * gens is a list of nonzero integers, positive for the positive generator
      and negative for the negative generator
    * n is the number of strands. By default it is inferred to be the least
      number of strands that works for the given list of generators

    >>> BraidTangle([], 1)
    <Tangle: BraidTangle([], 1): 1 comp; 0 cross; (1, 1) boundary>
    >>> BraidTangle([1]).describe()
    'Tangle[{1,2}, {3,4}, X[2,4,3,1]]'
    >>> BraidTangle([-1]).describe()
    'Tangle[{1,2}, {3,4}, X[1,2,4,3]]'
    >>> BraidTangle([1],3).describe()
    'Tangle[{1,2,3}, {4,5,6}, X[2,5,4,1], P[3,6]]'
    >>> BraidTangle([2],3).describe()
    'Tangle[{1,2,3}, {4,5,6}, X[3,6,5,2], P[1,4]]'
    >>> BraidTangle([1,2]).describe()
    'Tangle[{1,2,3}, {4,5,6}, X[7,5,4,1], X[3,6,7,2]]'
    >>> BraidTangle([1,2,1]).describe()
    'Tangle[{1,2,3}, {4,5,6}, X[7,5,4,8], X[3,6,7,9], X[2,9,8,1]]'
    >>> BraidTangle([1,2,1])
    <Tangle: BraidTangle([1, 2, 1], 3): 3 comp; 3 cross; (3, 3) boundary>
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

    b.make_upward()
    b.update_label(f'BraidTangle({gens}, {n})')

    return b
