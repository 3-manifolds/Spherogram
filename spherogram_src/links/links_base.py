from __future__ import print_function
"""
Links are made from Crossings.  The general model is that of
a PD diagram as described in 

  http://katlas.org/wiki/Planar_Diagrams

See the file "doc.pdf" for the conventions, and the file
"test.py" for some examples of creating links.
"""

no_snappy_msg = 'requires that SnapPy be installed.'

from .. import graphs
CyclicList = graphs.CyclicList
import  string, os, sys, re, collections
try:
    import cPickle as pickle
except ImportError: # Python 3
    import pickle

import copy

def is_iterable(obj):
    try:
        iter(obj)
        return True
    except TypeError:
        return False

class Crossing(object):
    """
    See "doc.pdf" for the conventions.  The sign of a crossing can be in {0,
    +1, -1}.  In the first case, the strands at the crossings can have
    any orientations, but crossings with sign +/-1 must be oriented as
    shown in "doc.pdf".
    
    Roles of some of the other attributes:
    
    * label: Arbitrary name used for printing the crossing.
    
    * directions: store the orientations of the link components passing
    through the crossing.  For a +1 crossing this is { (0, 2), (3, 1) }.
    Set with calls to make_tail.
    
    * adjacent: The other Crossings that this Crossing is attached to.
    
    * strand_labels: Numbering of the strands, used for DT codes and
    such.
    
    * strand_components: Which element of the parent
    Link.link_components each input's strand belongs to.
    """
    def __init__(self, label=None):
        self.label = label
        self.adjacent = CyclicList([None, None, None, None])
        self._clear()

    def _clear(self):
        self.sign, self.directions = 0, set()
        self._clear_strand_info()

    def _clear_strand_info(self):
        self.strand_labels = CyclicList([None, None, None, None])
        self.strand_components = CyclicList([None, None, None, None])

    def make_tail(self, a):
        """
        Orients the strand joining input "a" to input" a+2" to start at "a" and end at
        "a+2".
        """
        b = (a, (a + 2) % 4)
        if (b[1], b[0]) in self.directions:
            raise ValueError("Can only orient a strand once.")
        self.directions.add(b)

    def rotate(self, s):
        """
        Rotate the incoming connections by 90*s degrees anticlockwise.  
        """
        rotate = lambda v : (v + s) % 4
        new_adjacent = CyclicList(self.adjacent[s:] + self.adjacent[:s])
        for i, (o, j) in enumerate(new_adjacent):
            if o != self:
                o.adjacent[j] = (self, i)
                self.adjacent[i] = (o, j)
            else:
                self.adjacent[i] = (self, (j - s) % 4)

        self.directions = set( [ (rotate(a), rotate(b)) for a, b in self.directions] )

    def rotate_by_90(self):
        "Effectively switches the crossing"
        self.rotate(1)
        
    def rotate_by_180(self):
        "Effective reverses directions of the strands"
        self.rotate(2)
    
    def orient(self):
        if (2, 0) in self.directions:
            self.rotate_by_180()
        self.sign = 1 if (3, 1) in self.directions else -1

    def __getitem__(self, i):
        return (self, i%4)

    def entry_points(self):
        verts = [0,1] if self.sign == -1 else [0, 3]
        return [CrossingEntryPoint(self, v) for v in verts]

    def crossing_strands(self):
        return [CrossingEntryPoint(self, v) for v in range(4)]
    
    def __setitem__(self, i, other):
        o, j = other
        self.adjacent[i % 4] = other
        other[0].adjacent[other[1]] = (self, i % 4)
        
    def __repr__(self):
        return "%s" % self.label

    def info(self):
        def format_adjacent(a):
            return (a[0].label, a[1]) if a else None
        print( "<%s : %s : %s : %s>" % (self.label, self.sign, [format_adjacent(a) for a in self.adjacent], self.directions) )

    def DT_info(self):
        labels = self.strand_labels
        over = labels[3]+1 if self.sign == 1 else labels[1] + 1
        under = labels[0]+1
        return (under, -over) if over % 2 == 0 else (over, under)

    def peer_info(self):
        labels = self.strand_labels
        SW = labels[0] if self.sign == 1 else labels[1] 
        NW = labels[3] if self.sign == 1 else labels[0]
        if SW % 2 == 0:
            ans = SW, (-NW, -self.sign)
        else:
            ans = NW, (SW, self.sign)

        return ans

BasicCrossingStrand = collections.namedtuple('BasicCrossingStrand', ['crossing', 'entry_point'])

class CrossingStrand(BasicCrossingStrand):
    """
    One of the four incoming strands at a crossing.
    """
    def rotate(self, s=1):
        """
        The CrossingStrand *counter-clockwise* from self.
        """
        return CrossingStrand(self.crossing, (self.entry_point + s) % len(self.crossing.adjacent))

    def opposite(self):
        """
        The CrossingStrand at the other end of the edge from self
        """
        return CrossingStrand(* self.crossing.adjacent[self.entry_point] )

    def next_corner(self):
        return self.rotate().opposite()

    def previous_corner(self):
        return self.opposite().rotate(-1)

    def strand_label(self):
        return self.crossing.strand_labels[self.entry_point]

    def oriented(self):
        """
        Returns the one of {self, opposite} which is the *head* of the
        corresponding oriented edge of the link.
        """
        c, e = self.crossing, self.entry_point
        if (c.sign == 1 and e in [0, 3]) or (c.sign == -1 and e in [0, 1]):
            return self
        return self.opposite()

    def __repr__(self):
        return "<CS %s, %s>" % (self.crossing, self.entry_point)
    
class CrossingEntryPoint(CrossingStrand):
    """
    One of the two entry points of an oriented crossing
    """    
    def next(self):
        c, e = self.crossing, self.entry_point
        s = 1 if isinstance(c, Strand) else 2
        return CrossingEntryPoint(*c.adjacent[ (e + s) % (2*s)])

    def other(self):
        nonzero_entry_point = 1 if self.crossing.sign == -1 else 3
        other = nonzero_entry_point if self.entry_point == 0 else 0
        return CrossingEntryPoint(self.crossing, other)

    def is_under_crossing(self):
        return self.entry_point == 0

    def is_over_crossing(self):
        return self.entry_point != 0
        
    def component(self):
        ans = [self]
        while True:
            next = ans[-1].next()
            if next == self:
                break
            else:
                ans.append(next)

        return ans

    def component_label(self):
        return self.crossing.strand_components[self.entry_point]      

    def label_crossing(self, comp, labels):
        c, e = self.crossing, self.entry_point
        f = (e + 2) % 4
        c.strand_labels[e], c.strand_components[e] = labels[self], comp
        c.strand_labels[f], c.strand_components[f] = labels[self.next()], comp

    def __repr__(self):
        return "<CEP %s, %s>" % (self.crossing, self.entry_point)
        

class LinkComponents(list):
    def add(self, c):
        component = c.component()
        self.append(component)
        return component

class Labels(dict):
    def add(self, x):
        self[x] = len(self)

class Strand(object):
    """
    When constructing links, it's convenient to have strands as well
    as crossings.  These are stripped by the Link class when it
    pieces things together.  
    """
    def __init__(self, label=None):
        self.label = label
        self.adjacent = [None, None]

    def fuse(self):
        """
        Joins the incoming and outgoing strands and removes
        self from the picture.
        """
        (a, i), (b,j) = self.adjacent
        a.adjacent[i] = (b, j)
        b.adjacent[j] = (a, i)

    def __getitem__(self, i):
        return (self, i%2)

    def __setitem__(self, i, other):
        o, j = other
        self.adjacent[i % 2] = other
        other[0].adjacent[other[1]] = (self, i)

    def __repr__(self):
        return "%s" % self.label

    def info(self):
        def format_adjacent(a):
            return (a[0].label, a[1]) if a else None
        print( "<%s : %s>" % (self.label, [format_adjacent(a) for a in self.adjacent]) )
    
def enumerate_lists(lists, n=0, filter=lambda x:True):
    ans = []
    for L in lists:
        ans.append([n + i for i, x in enumerate(L) if filter(n+i)])
        n += len(L)
    return ans



class Link(object):
    """
    Links are made from Crossings.  The general model is that of the PD
    diagrams used in `KnotTheory <http://katlas.org/wiki/Planar_Diagrams>`_.
    
    See the file "doc.pdf" for the conventions, which can be accessed
    via "spherogram.pdf_docs()", and the `Spherogram tutorial
    <http://www.math.uic.edu/t3m/SnapPy/spherogram.html>`_
    for some examples of creating links.

    Here are two ways of creating the figure-8 knot, first via a PD code

    >>> K1 = Link([[8,3,1,4],[2,6,3,5],[6,2,7,1],[4,7,5,8]])

    and by directly gluing up Crossings:
    
    >>> a, b, c, d = [Crossing(x) for x in 'abcd']
    >>> a[0], a[1], a[2], a[3] = c[1], d[0], b[1], b[0]
    >>> b[2], b[3] = d[3], c[2]
    >>> c[3], c[0] = d[2], d[1]
    >>> K2 = Link([a,b,c,d])

    Some families of named links are available, such a torus knots

    >>> Link('T(4, 2)')
    <Link: 2 comp; 6 cross>

    and if you have SnapPy installed also the Rolfsen and Hoste-Thistlethwaite
    tables:

    >>> Link('8_20')
    <Link 8_20: 1 comp; 8 cross>
    >>> Link('K12a123')
    <Link K12a123: 1 comp; 12 cross>
    >>> Link('L12n123')
    <Link L12n123: 2 comp; 12 cross>
    """
    
    def __init__(self, crossings, check_planarity=True, build=True):
        # We check if crossings is a string
        self.name = None
        if isinstance(crossings, str):
            if crossings.startswith('T(' ):
                from . import torus
                crossings = torus.torus_knot(crossings).crossings
            else:
                try:
                    import snappy
                    self.name = crossings
                    ans = snappy.Manifold(self.name).link()
                    build = False
                    crossings = ans.crossings
                    self.link_components = ans.link_components
                except ImportError:
                    raise RuntimeError('creating a Link object with argument of type str '+no_snappy_msg)
        
        #If crossings is just a PD code rather than a list of Crossings,
        # we create the corresponding Crossings.
        if len(crossings) > 0 and not isinstance(crossings[0], (Strand, Crossing)):
            gluings = collections.defaultdict(list)
            for c, X in enumerate(crossings):
                for i, x in enumerate(X):
                    gluings[x].append( (c,i) )

            if set([len(v) for v in gluings.values()]) != set([2]):
                raise ValueError("PD code isn't consistent")
             
            crossings = [Crossing() for d in crossings]
            for (c,i), (d,j) in gluings.itervalues():
                crossings[c][i] = crossings[d][j]

        # Make sure everything is tied up. 
        if True in [ None in c.adjacent for c in crossings]:
            raise ValueError("No loose strands allowed")

        # Fuse the strands.  If there any components made up
        # only of strands, these are thrown out here.
        [s.fuse() for s in crossings if isinstance(s, Strand)]
        self.crossings = [c for c in crossings if not isinstance(c, Strand)]

        if build:
            self._build()

        if check_planarity and not self.is_planar():
            raise ValueError("Link isn't planar")

        # If the crossings aren't labeled the label them for
        # debugging purposes.
        if False not in [X.label is None for X in self.crossings]:
            for c, X in enumerate(self.crossings):
                X.label = c

    def all_crossings_oriented(self):
        return len([c for c in self.crossings if c.sign == 0]) == 0

    def _build(self):
        self._orient_crossings()
        self._build_components()

    def _rebuild(self):
        self.link_components = None
        for c in self.crossings:
            c._clear()
        self._build()
        
    def _orient_crossings(self):
        if self.all_crossings_oriented():
            return
        
        remaining = set( [ (c, i) for c in self.crossings for i in range(4)] )
        while len(remaining):
            c, i = start = remaining.pop()
            finished = False
            while not finished:
                d, j = c.adjacent[i]
                d.make_tail(j)
                remaining.discard( (c,i) ), remaining.discard( (d,j) )
                c, i = d, (j+2) % 4
                finished = (c, i) == start

        for c in self.crossings:
            c.orient()

    def crossing_entries(self):
        ans = []
        for C in self.crossings:
            ans += C.entry_points()
        return ans

    def crossing_strands(self):
        ans = []
        for C in self.crossings:
            ans += C.crossing_strands()
        return ans
    
    def _build_components(self, component_starts=None):
        """
        Each component is stored as a list of *entry points* to
        crossings.  The labeling of the entry points (equivalently
        oriented edges) is compatible with the DT convention that each
        crossing has both an odd and and even incoming strand.
        """
        remaining, components = set( self.crossing_entries() ), LinkComponents()
        other_crossing_entries = []
        self.labels = labels = Labels()
        for c in self.crossings:
            c._clear_strand_info()
        
        while len(remaining):
            if component_starts:
                d = component_starts[len(components)]
            elif len(components) == 0:
                d = remaining.pop()
            else:
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
            for c in component:
                labels.add( c )
            others = []
            for c in component:
                c.label_crossing(len(components) - 1, labels)
                o = c.other()
                if o.component_label() is None:
                    others.append(o)
            other_crossing_entries.append(others)
            remaining = remaining - set(component)

        self.link_components = components

    def digraph(self):
        """
        The underlying directed graph for the link diagram. 
        """
        G = graphs.Digraph()
        for component in self.link_components:
            for c in component:
                cs0 = CrossingStrand(c.crossing, c.entry_point)
                cs1 = cs0.opposite()
                e = G.add_edge(cs0.crossing, cs1.crossing)
        return G

    def copy(self):
        return pickle.loads(pickle.dumps(self))
    
    def is_planar(self):
        """
        Whether the 4-valent graph underlying the link projection is planar.
        Should always be True for any actual Link.

        >>> L = Link('5^2_1')
        >>> L.is_planar()
        True
        >>> c = Crossing()
        >>> c[0], c[1] = c[2], c[3]   # Punctured torus gluing
        >>> bad = Link([c], check_planarity=False)
        >>> bad.is_planar()
        False
        """
        G = self.digraph()
        if not G.is_weakly_connected():
            return False
        v = len(self.crossings)
        assert 2*v == len(G.edges)
        euler = -v + len(self.faces())
        return euler == 2 or v == 0

    def faces(self):
        """
        The faces are the complementary regions of the link diagram. Each face
        is given as a list of corners of crossings as one goes around
        *clockwise*.  These corners are recorded as CrossingStrands,
        where CrossingStrand(c, j) denotes the corner of the face
        abutting crossing c between strand j and j + 1.
        
        Alternatively, the sequence of CrossingStrands can be regarded
        as the *heads* of the oriented edges of the face.
        """
        corners = set([ CrossingStrand(c,i) for c in self.crossings for i in range(4) ])
        faces = []
        while len(corners):
            face = [corners.pop()]
            while True:
                next = face[-1].next_corner()
                if next == face[0]:
                    faces.append(face)
                    break
                else:
                    corners.remove(next)
                    face.append(next)

        return faces

    def dual_graph(self):
        """
        The dual graph to a link diagram D, whose vertices correspond to
        complementary regions (faces) of D and whose edges are dual to the
        edges of D.
        """
        from . import simplify
        return simplify.DualGraphOfFaces(self)

    def simplify(self, mode='basic', type_III_limit=100):
        """
        Tries to simplify the link projection.  Returns whether it succeeded
        in reducing the number of crossings.  Modifies the link in
        place, and unknot components which are also unlinked may be
        silently discarded.  The ordering of ``link_components`` is not
        always preserved.

        The following strategies can be employed.

        1. In the default ``basic`` mode, it does Reidemeister I and II moves
           until none are possible.

        2. In ``level`` mode, it does random Reidemeister III moves, reducing
           the number of crossings via type I and II moves whenever possible.
           The process stops when it has done ``type_III_limit`` *consecutive*
           type III moves without any simplification.

        3. In ``pickup`` mode, it also minimizes the number of crossings of
           strands which cross entirely above (or below) the diagram by
           finding the path crossing over the diagram with the least number of
           overcrossings (or undercrossings); this has the effect of doing
           "picking up" strands and putting them down elsewhere.

        4. Finally, the ``global`` mode is the combination of 3 and 4. 


        Some examples:

        >>> K = Link([(13,10,14,11),(11,5,12,4),(3,13,4,12),  \
                 (9,14,10,1),(1,7,2,6),(2,7,3,8),(5,9,6,8)])
        >>> K
        <Link: 1 comp; 7 cross>
        >>> K.simplify('basic')
        True
        >>> K
        <Link: 1 comp; 4 cross>
        >>> K.simplify('basic')  # Already done all it can
        False

        >>> L = Link([(5,0,6,1), (14,5,15,4), (10,2,11,3), (7,12,8,11), \
                (17,0,14,9), (12,9,13,8), (3,13,4,10), (1,16,2,15), (16,6,17,7)])
        >>> L
        <Link: 3 comp; 9 cross>
        >>> L.simplify('basic')
        False
        >>> L.simplify('level')
        True
        >>> L    # Trivial unlinked component has been discarded!
        <Link: 2 comp; 2 cross>

        >>> K = Link('K14n2345')
        >>> K.backtrack(30) 
        >>> K.simplify('global')
        True
        >>> K
        <Link: 1 comp; 14 cross>
        """
        from . import simplify
        if mode == 'basic':
            return simplify.basic_simplify(self)
        elif mode == 'level':
            return simplify.simplify_via_level_type_III(self, type_III_limit)
        elif mode == 'pickup':
            return simplify.pickup_simplify(self)
        elif mode == 'global':
            return simplify.pickup_simplify(self, type_III_limit)

    def backtrack(self, steps=10):
        """
        Performs a sequence of Reidemeister moves which increase or maintain
        the number of crossings in a diagram.  The number of such
        moves is the parameter steps.  The diagram is modified in place. 

        >>> K = Link('L14a7689')
        >>> K
        <Link L14a7689: 2 comp; 14 cross>
        >>> K.backtrack(steps = 50)
        >>> len(K.crossings) > 20
        True
        """
        from . import simplify
        L = simplify.backtrack(self,num_steps = steps)
        self.crossings = L.crossings
        self.labels = L.labels
        self.link_components = L.link_components
        self.name = L.name


    def sublink(self, components):
        """
        Returns the sublink consisting of the specified components; see the
        example below for the various accepted forms.

        Warnings: Components in the sublink that are both unknotted
        and unlinked may be silently thrown away.  The order of the
        components in the sublink need not correspond to their order
        in the original link.

        >>> L = Link('L14n64110')
        >>> L
        <Link L14n64110: 5 comp; 14 cross>
        >>> L.sublink([1,2,3,4])
        <Link: 4 comp; 10 cross>
        >>> comps = L.link_components
        >>> L.sublink([comps[0], comps[1]])
        <Link: 2 comp; 2 cross>

        If you just want one component you can do this:

        >>> L = Link('L11a127')
        >>> L.sublink(1)
        <Link: 1 comp; 7 cross>
        >>> L.sublink(L.link_components[0])
        <Link: 0 comp; 0 cross>

        The last answer is because the second component is unknotted
        and so thown away.
        """

        if components in self.link_components or not is_iterable(components):
            components = [components]
        indices = []
        for c in components:
            if is_iterable(c):
                c = self.link_components.index(c)
            else:
                try:
                    self.link_components[c]
                except IndexError:
                    raise ValueError('No component of that index')
            indices.append(c)

        def keep(C):
            return set([i in indices for i in C.strand_components]) == set([True])
        
        L = self.copy()
        final_crossings = []
        for C in L.crossings:
            if keep(C):
                final_crossings.append(C)
            else:
                for j in [0, 1]:
                    if C.strand_components[j] in indices:
                        A, a = C.adjacent[j]
                        B, b = C.adjacent[j + 2]
                        A[a] = B[b]

        return type(self)(final_crossings, check_planarity=False)


    def __len__(self):
        return len(self.crossings)

    def PD_code(self, KnotTheory=False, min_strand_index=0):
        """
        The planar diagram code for the link.  
        """
        PD = []
        for c in self.crossings:
            PD.append([s + min_strand_index for s in c.strand_labels[:]])
        if KnotTheory:
            PD = "PD" + repr(PD).replace('[', 'X[')[1:]
        else:
            PD = [tuple(x) for x in PD ]
        return PD

    def DT_code(self, DT_alpha=False):
        """
        The Dowker-Thistlethwaite code for the link in either numerical or
        alphabetical form.
        """
        DT_dict = dict( [ c.DT_info() for c in self.crossings] )
        odd_labels = enumerate_lists(self.link_components, n=1, filter=lambda x:x%2==1)
        DT = [ tuple([DT_dict[x] for x in component]) for component in odd_labels]

        if DT_alpha:
            if len(self) > 52:
                raise ValueError("Too many crossing for alphabetic DT code")
            DT_alphabet = '_abcdefghijklmnopqrstuvwxyzZYXWVUTSRQPONMLKJIHGFEDCBA'
            init_data = [len(self), len(DT)] +  [len(c) for c in DT]
            DT = ''.join([DT_alphabet[x] for x in init_data] + [DT_alphabet[x>>1] for x in sum(DT, tuple())])
            DT = "DT[" + DT + "]"

        return DT

    def peer_code(self):
        peer = dict( [ c.peer_info() for c in self.crossings] )
        even_labels = enumerate_lists(self.link_components, n=0, filter=lambda x:x%2==0)
        ans = '[' + ','.join([repr([peer[c][0] for c in comp])[1:-1].replace(',', '') for comp in even_labels])
        ans += '] / ' + ' '.join([ ['_', '+', '-'][peer[c][1]] for c in sum(even_labels, [])])
        return ans
        
    def KLPProjection(self):
        return python_KLP(self)

    def split_link_diagram(self, destroy_original=False):
        """
        Breaks the given link diagram into pieces, one for each connected
        component of the underlying 4-valent graph.

        >>> L = Link([(2,1,1,2),(4,3,3,4)], check_planarity=False)
        >>> L.split_link_diagram()
        [<Link: 1 comp; 1 cross>, <Link: 1 comp; 1 cross>]
        """
        link = self.copy() if not destroy_original else self
        return [type(self)(list(component), check_planarity=False)
                for component in link.digraph().weak_components()]

    def deconnect_sum(self, destroy_original=False):
        """
        Undoes all connect sums that are diagramatically obvious,
        i.e. those where there is a circle which meets the projection
        in two points.

        >>> K = Link('5_2')
        >>> L = K.connected_sum(K); L
        <Link: 1 comp; 10 cross>
        >>> L.deconnect_sum()
        [<Link: 1 comp; 5 cross>, <Link: 1 comp; 5 cross>]
        """
        from . import simplify
        link = self.copy() if not destroy_original else self
        return simplify.deconnect_sum(link)
    
    def exterior(self):
        raise RuntimeError("SnapPy doesn't seem to be available.  Try: from snappy import *")

    def __repr__(self):
        name = ' ' + self.name if self.name else ''
        return "<Link%s: %d comp; %d cross>" % (name, len(self.link_components), len(self.crossings))

    def writhe(self):
        """
        Finds the writhe of a knot.

        >>> K = Link( [(4,1,5,2),(6,4,7,3),(8,5,1,6),(2,8,3,7)] )  # Figure 8 knot
        >>> K.writhe()
        0
        """
        writhe_value=0
        for i in range(len(self.crossings)):
                writhe_value+=self.crossings[i].sign
        return writhe_value

    def linking_number(self):
        """
        Returns the linking number of self if self has two components;
        or the sum of the linking numbers of all pairs of components 
        in general.
        """
        n = 0
        for s in self.link_components:
            tally = [0]*len(self.crossings)
            for c in s:
                for i, x in enumerate(self.crossings):
                    if c[0] == x:
                        tally[i] += 1
            for i, m in enumerate(tally):
                if m == 1:
                    n += (self.crossings)[i].sign
        n = n/4
        return n

    def pieces(self):
        """
        Auxiliary function used by knot_group. Constructs the strands 
        of the knot from under-crossing to under-crossing. Needed for the
        Wirtinger Presentation.
        """
        pieces = []
        for s, x in enumerate(self.crossings):
            y = x
            l = 2
            go = True
            pieces.append([])
            pieces[s].append(x[l])
            while go:
                if y.adjacent[l][1] == 0:
                    pieces[s].append(y.adjacent[l][0][0])
                    break
                pieces[s].append(y.adjacent[l])
                if y.adjacent[l][1] == 1:
                    pieces[s].append(y.adjacent[l][0][3])
                if y.adjacent[l][1] == 3:
                    pieces[s].append(y.adjacent[l][0][1])
                lnew = (y.adjacent[l][1] + 2)%4
                y = y.adjacent[l][0]
                l = lnew
        return pieces

    def connected_sum(self, other_knot):
        """
        Returns the connected sum of two knots.                                                       
       
        >>> K = Link('4_1')
        >>> K.connected_sum(K)                                                                         
        <Link: 1 comp; 8 cross>
        """
        first=self.copy()
        second=other_knot.copy()
        (f1,i1) = (first.crossings[0],0)
        (f2,i2) = f1.adjacent[i1]
        (g1,j1) = (second.crossings[0],0)
        (g2,j2) = g1.adjacent[j1]
        f1[i1]=g2[j2]
        f2[i2]=g1[j1]
        for c in first.crossings: #Indicate which crossings came from which knot.                        
            c.label=(c.label, 1)
        for c in second.crossings:
            c.label=(c.label, 2)
            first.crossings.append(c)
        return type(self)(first.crossings)

    def copy(self):
        """
        Returns a copy of the link.  

        >>> K=Link('4_1')
        >>> copy=K.copy()
        >>> K.PD_code() == copy.PD_code()
        True

        If the link is very large (100s of crossings), you may get a
        recursion limit exception; to get around this, call e.g.
        ``sys.setrecursionlimit(10000)``.
        """
        return copy.deepcopy(self)

    def mirror(self):
        """
        Returns the mirror image of the link, preserving link orientations and
        component order.
        """
        # Basically, we just mirror every crossing, but the particular data 
        # structures used make this a little involved.

        new_crossings = dict()
        for C in self.crossings:
            C_new = Crossing(label=C.label)
            C_new.sign = -C.sign
            new_crossings[C] = C_new

        def convert(C, c):
            """
            Go from a crossing to the mirror, which requires the a rotation of the
            entry points; the direction of rotation depends on the sign of
            the crossing.
            """
            return new_crossings[C], (c + C.sign) % 4

        for A in self.crossings:
            entry_points = [CEP.entry_point for CEP in A.entry_points()]
            for a in entry_points:
                B, b = A.adjacent[a]
                B_new, b_new = convert(B, b)
                B_new[b_new] = convert(A, a)

        new_link = type(self)(new_crossings.values(),
                        check_planarity=False, build=False)

        # Build the link components, starting in the same place as
        # the original link.

        component_starts = []
        for component in self.link_components:
            C_new, c_new = convert(*component[0])
            component_starts.append( CrossingEntryPoint(C_new, c_new) )

        new_link._build_components(component_starts)
        return new_link

    def alternating(self):
        """
        Returns the alternating link with the same planar graph.  No attempt
        is made to preserve the order of the link components or ensure
        that the DT code of the result has all positive entries (as
        opposed to all negative).

        >>> L = Link('L14n12345')
        >>> A = L.alternating()
        >>> A.exterior().identify()
        [L14a5150(0,0)(0,0)]
        """
        L = self.copy()
        for C in L.crossings:
            a, b = C.DT_info()
            if a*b < 0:
                C.rotate_by_90()
        L._rebuild()
        return L        
            
    def optimize_overcrossings(self):
        """
        Minimizes the number of crossings of a strand which crosses entirely
        above the diagram by finding the path crossing over the diagram with
        the least number of overcrossings.  It begins with the longest
        overcrossing, and continues with smaller ones until it successfully
        reduces the number of crossings.
        """
        from . import simplify
        return simplify.strand_pickup(self, self.overcrossings())

    def overcrossings(self):
        """
        Returns a list of the sequences of overcrossings of the link,
        sorted in descending order of length
        """
        comp = self.link_components
        overcrosses = []
        for c in comp:
            origin = None
            for i in c:
                if i.is_under_crossing():
                    origin = i
                    break
            start = origin
            end = start.next()
            circled_around = False
            while True:
                length = 0
                alreadyseen = []
                while (end.is_over_crossing() or (end.crossing in alreadyseen)) and (end.crossing != start.crossing):
                    alreadyseen.append(end.crossing)
                    length += 1
                    end = end.next()
                    if end == origin:
                        circled_around = True
                if length >= 1:
                    overcrosses.append([start,length])
                start = end
                end = start.next()
                circled_around = circled_around or (start == origin)
                if circled_around:
                    break
        return sorted(overcrosses, key=lambda overcross: overcross[1], reverse=True)





# ---- building the link exterior if SnapPy is present --------

def vertex_to_KLP(c, v):
    i = v if c.sign == 1 else (v - 1) % 4
    return ['Ybackward', 'Xforward', 'Yforward', 'Xbackward'][i]

class KLPCrossing(object):
    """
    SnapPea uses a convention where the orientation
    of the strands is fixed in the master picture but
    which strand is on top varies.
    """
    def __init__(self, c):
        self.adjacent = adjacent = 4*[None]
        self.index = c._KLP_index
        if c.sign == 1:
            strands, self.sign = [3, 0], 'R'
        else:
            strands, self.sign = [0,1], 'L'
            
        components = [ c.strand_components[s] for s in strands]
        self.Xcomponent, self.Ycomponent = components
        self.strand, self.neighbor = {}, {}
        for v in range(4):
            d, w = c.adjacent[v]
            self.neighbor[vertex_to_KLP(c, v)] = d._KLP_index
            self.strand[vertex_to_KLP(c, v)] = vertex_to_KLP(d, w)[:1]

    def __getitem__(self, index):
        if index.find('_') == -1:
            return getattr(self, index)
        vertex, info_type = index.split('_')
        return getattr(self, info_type)[vertex]

def python_KLP(L):
    vertices = list(L.crossings)
    if len(vertices) == 0:  # Unknot by convention
        return [0, 1, 1, []]
    for i, v in enumerate(vertices):
        v._KLP_index = i
    return [len(vertices), 0, len(L.link_components), [KLPCrossing(c) for c in vertices]]
