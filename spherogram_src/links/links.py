from __future__ import print_function
"""
Links are made from Crossings.  The general model is that of
a PD diagram as described in 

  http://katlas.org/wiki/Planar_Diagrams

See the file "doc.pdf" for the conventions, and the file
"test.py" for some examples of creating links.

"""

#Check if being run inside sage; if not, some functionality disabled.
_within_sage = False
try:
    import sage.all
    import sage.graphs.graph
    _within_sage = True
    from sage.matrix.constructor import matrix
    from sage.symbolic.ring import SR
    from sage.groups.free_group import FreeGroup
    import sage.graphs.graph as graph
    from sage.symbolic.ring import var
except ImportError:
    pass
not_in_sage_msg = 'is only available when running Spherogram inside Sage.'
no_snappy_msg = 'requires that SnapPy be installed.'

from .. import graphs
CyclicList = graphs.CyclicList
import  string, os, sys, re, collections
try:
    import cPickle as pickle
except ImportError: # Python 3
    import pickle

import copy

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
        self.label, self.sign, self.directions = label, 0, set()
        self.adjacent = CyclicList([None, None, None, None])
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
    Links are made from Crossings.  The general model is that of a PD
    diagram as described in
    
    http://katlas.org/wiki/Planar_Diagrams
    
    See the file "doc.pdf" for the conventions, which can be accessed
    via "spherogram.pdf_docs()", and the file "test.py" for some
    examples of creating links.

    Here are two ways of creating the figure-8 knot, first via a PD code

    >>> K1 = Link([[8,3,1,4],[2,6,3,5],[6,2,7,1],[4,7,5,8]])

    and by directly gluing up Crossings:
    
    >>> a, b, c, d = [Crossing(x) for x in 'abcd']
    >>> a[0], a[1], a[2], a[3] = c[1], d[0], b[1], b[0]
    >>> b[2], b[3] = d[3], c[2]
    >>> c[3], c[0] = d[2], d[1]
    >>> K2 = Link([a,b,c,d])
    """
    
    def __init__(self, crossings, check_planarity=True, build=True):
        # We check if crossings is a string
        self.name = None
        if isinstance(crossings, str):
            try:
                import snappy
                self.name = crossings
                crossings = (snappy.Manifold(crossings)).link().crossings
            except ImportError:
                raise RunTimeError('creating a Link object with argument of type str '+no_snappy_msg)
        
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
            self._orient_crossings()
            self._build_components()

        if check_planarity and not self.is_planar():
            raise ValueError("Link isn't planar")

        # If the crossings aren't labeled the label them for
        # debugging purposes.
        if False not in [X.label is None for X in self.crossings]:
            for c, X in enumerate(self.crossings):
                X.label = c

    def all_crossings_oriented(self):
        return len([c for c in self.crossings if c.sign == 0]) == 0
    
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

    def basic_simplify(self):
        """
        Do Reidemeister I and II moves until none are possible.  Modifies the
        link in place, and unknot components which are also unlinked may
        be silently discarded.
        
        >>> K = Link([(13,10,14,11),(11,5,12,4),(3,13,4,12),  \
                 (9,14,10,1),(1,7,2,6),(2,7,3,8),(5,9,6,8)])
        >>> K
        <Link: 1 comp; 7 cross>
        >>> K.basic_simplify()
        True
        >>> K
        <Link: 1 comp; 4 cross>
        >>> K.basic_simplify()  # Already done all it can
        False
        """
        from . import simplify
        return simplify.basic_simplify(self)

    def simplify(self, max_consecutive_failures=100):
        """
        Applies a series of Reidemeister type III moves to the link, simplifying
        it via type I and II moves whenever possible.  Continues until there
        are no type III moves or it has done the specified number of
        consecutive type III moves without reducing the crossing number.

        >>> K = Link([(5,0,6,1), (14,5,15,4), (10,2,11,3), (7,12,8,11), \
                (17,0,14,9), (12,9,13,8), (3,13,4,10), (1,16,2,15), (16,6,17,7)])
        >>> K
        <Link: 3 comp; 9 cross>
        >>> K.basic_simplify()
        False
        >>> K.simplify()
        True
        
        Modifies the link in place, and unknot components which are also
        unlinked may be silently discarded.

        >>> K
        <Link: 2 comp; 2 cross>
        """
        from . import simplify
        return simplify.simplify(self)

    def __len__(self):
        return len(self.crossings)

    def PD_code(self, KnotTheory=False, min_strand_index=0):
        PD = []
        for c in self.crossings:
            PD.append([s + min_strand_index for s in c.strand_labels[:]])
        if KnotTheory:
            PD = "PD" + repr(PD).replace('[', 'X[')[1:]
        else:
            PD = [tuple(x) for x in PD ]
        return PD

    def DT_code(self, DT_alpha=False):
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
        return [Link(list(component), check_planarity=False)
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
        link = self.copy() if not destroy_original else self
        from . import simplify
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

    def linking_matrix(self):
        """
        Calcluates the linking number for each pair of link components.
        
        Returns a linking matrix, in which the (i,j)th component is the linking                                    
        number of the ith and jth link components.
        """
        if not _within_sage:
            raise RuntimeError('linking_matrix '+not_in_sage_msg)
        matrix = [ [0 for i in range(len(self.link_components)) ] for j in range(len(self.link_components)) ]
        for n1, comp1 in enumerate(self.link_components):
            for n2, comp2 in enumerate(self.link_components):
                tally = [ [0 for m in range(len(self.crossings)) ] for n in range(2) ]
                if not (comp1 == comp2):
                    for i, c in enumerate(self.crossings):
                        for x1 in comp1:
                            if x1[0] == c:
                                tally[0][i] += 1
                        for x2 in comp2:
                            if x2[0] == c:
                                tally[1][i] += 1
                for k, c in enumerate(self.crossings):
                    if (tally[0][k] == 1 and tally[1][k] == 1):
                        matrix[n1][n2] += 0.5 * (c.sign)
        for i1, m1 in enumerate(matrix):
            for i2, m2 in enumerate(m1):
                matrix[i1][i2] = int(matrix[i1][i2])
        return matrix

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
                #print 'building strands'                                                               
                #print y.adjacent[l][1]                                                                 
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

    def knot_group(self):
        """
        Computes the knot group using the Wirtinger presentation. 
        Returns a finitely presented group.
        """
        if not _within_sage:
            raise RuntimeError('knot_group '+not_in_sage_msg)

        n = len(self.crossings)
        F = FreeGroup(n)
        g = list(F.gens())
        rels = []
        pieces = self.pieces()

        for z in self.crossings:
            for m, p in enumerate(pieces):
                for t, q in enumerate(p):
                    if q[0] == z:
                        if t == 0:
                            j = m
                        elif t == len(p)-1:
                            i = m
                        else:
                            k = m
            i+=1; j+=1; k+=1
            if z.sign > 0:
                r = F([-k,i,k,-j])
                #r = (g[k]**-1)*g[i]*g[k]*(g[j]**-1)                                                    
            if z.sign < 0:
                r = F([k,i,-k,-j])
                #r = g[k]*g[i]*(g[k]**-1)*(g[j]**-1)                                                    
            rels.append(r)

        G = F/rels
        return G

    def alexander_matrix(self):
        """
        Returns the alexander matrix of self.
        
        >>> L = Link([(3,0,4,1),(5,2,0,3),(1,4,2,5)])  # Trefoil
        >>> L.alexander_matrix()   #doctest: +SKIP
        [    -1 -t + 1      t]
        [     t     -1 -t + 1]
        [-t + 1      t     -1]
        """
        if not _within_sage:
            raise RuntimeError('alexander_matrix '+not_in_sage_msg)

        G = self.knot_group()
        B = G.alexander_matrix()
        m = B.nrows()
        n = B.ncols()
        C = matrix(SR,n,m)
        t = var('t')
        for i in range(n):
            for j in range(m):
                for k in B[j,i].terms():
                    p = 0
                    x = k.leading_item()
                    for w in x[0].syllables():
                        p = p+w[1]
                    y = (t**p)*x[1]
                    C[j,i] = C[j,i] + y
        return C #should this be a k by k minor? or full matrix?                                       
    
    def alexander_poly(self, v='no', method='wirt'):
        """
        Calculates the alexander polynomial of self. For links with one component,
        can evaluate the alexander polynomial at v.
       
        >>> K = Link('4a1')
        >>> K.alexander_poly()
        -t - 1/t + 3
        >>> K.alexander_poly(4)
        -5/4
        """
        if not _within_sage:
            raise RuntimeError('alexander_poly '+not_in_sage_msg)

        if method == 'snappy':
            try:
                import snappy
                return snappy.snap.alexander_polynomial(self.exterior())
            except ImportError:
                raise RunTimeError('this method for alexander_poly '+no_snappy_msg)
        else:
            t = var('t')
            C = self.alexander_matrix()
            m = C.nrows()
            n = C.ncols()
            if n>m:
                k = m-1
            else:
                k = n-1
            minor = C[0:k,0:k]
            p = minor.determinant().expand()
            exps = [x[1] for x in p.coeffs()]
            a = max(exps)
            b = min(exps)
            c = -1*(b+a)/2
            q = p*(t**c)
            if q.subs(t==1) == -1:
                q = -1*q
            if v=='no':
                return q.expand()
            else:
                return q.subs(t==v)

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
        return Link(first.crossings)

    def black_graph(self, return_signs=False):
        """
        Finds the black graph for a knot, and returns just
        one component of the graph.                
       
        >>> K=Link('5_1')                                                                                
        >>> K.black_graph()
        Subgraph of (): Multi-graph on 2 vertices
        """
        # this is a bit hacky, could stand to be re-written

        faces = []
        for x in self.faces():
            l = []
            for y in x:
                l.append((y[0],y[1]))
                l.append((y[0],(y[1]+1)%4))
            faces.append(l)

        coords=list()
        signs=list()
        for i in range(len(faces)-1):
            for j in range (i+1, len(faces)):
                a=set(faces[i])
                b=set(faces[j])
                s=a.union(b)
                for x in range(len(self.crossings)):
                    crossings=[self.crossings[x][0],self.crossings[x][1],self.crossings[x][2],self.crossings[x][3]]
                    total=set(crossings)
                    if total.issubset(s):
                        coords.append((tuple(faces[i]),tuple(faces[j])))
                        if set([self.crossings[x][1], self.crossings[x][2]]).issubset(set(faces[i])) or set([self.crossings[x][3], self.crossings[x][0]]).issubset(set(faces[i])):
                                signs.append(-1)
                        elif set([self.crossings[x][2], self.crossings[x][3]]).issubset(set(faces[i])) or set([self.crossings[x][0], self.crossings[x][1]]).issubset(set(faces[i])):
                                signs.append(1)

        G=graph.Graph(coords)
        component=G.connected_components()[1]
        G=G.subgraph(component)
        #Built shorter versions of coords and signs corresponding just to those edges in the subgraph:   
        new_coords = list()
        new_signs = list()
        edges = G.edges()
        for n in range(len(coords)):
            if coords[n]+(None,) in edges:
                new_coords.append(coords[n])
                new_signs.append(signs[n])
            if (coords[n][1],coords[n][0],None) in edges:
                new_coords.append((coords[n][1],coords[n][0]))
                new_signs.append(signs[n])
        if return_signs==True:
            return (new_signs,new_coords)
        else:
            return G

    def goeritz_matrix(self):
        """
        Finds the black graph of a knot, and from that, returns the Goeritz matrix of that knot.
       
        >>> K=Link('4_1')
        >>> K.goeritz_matrix()   #doctest: +SKIP
        [-3  2]
        [ 2 -3]
        """
        if not _within_sage:
            raise RuntimeError('goeritz_matrix '+not_in_sage_msg)

        (all_signs,all_edges)=self.black_graph(True)
        g=self.black_graph()
        l=g.vertices()
        m=matrix(len(g.vertices()),len(g.vertices()))
        for n in range(len(all_edges)):
            e = all_edges[n]
            i=l.index(e[0])
            j=l.index(e[1])
            m[(i,j)]=m[(i,j)]+all_signs[n]
            m[(j,i)]=m[(j,i)]+all_signs[n]
        for i in range(len(g.vertices())):
            m[(i,i)]=sum(m.column(i))*(-1)
        m=m.delete_rows([0])
        m=m.delete_columns([0])
        return m

    def _find_crossing(self, e):
        """
        Auxiliary function used by signature to find which
        crossing corresponds to an edge in the black graph.
        """
        a=set(e[0])
        b=set(e[1])
        s=a.union(b)
        for x in range(len(self.crossings)):
            crossings=[self.crossings[x][0],self.crossings[x][1],self.crossings[x][2],self.crossings[x][3]]
            total=set(crossings)
            if total.issubset(s):
                return self.crossings[x]

    def signature(self):
        """
        Returns the signature of the link.       
        >>> K = Link('4a1')                                                                               
        >>> K.signature()                                                                              
        0
        """
        if not _within_sage:
            raise RuntimeError('signature '+not_in_sage_msg)

        answer=0
        (signs,edges)=self.black_graph(True)
        for i in range(len(signs)):
            v=self._find_crossing(edges[i])
            if signs[i]*v.sign==1:
                answer=answer+signs[i]
        m=self.goeritz_matrix()
        vals=m.eigenvalues()
        pos=0
        neg=0
        for v in vals:
            if v>0:
                pos+=1
            elif v<0:
                neg+=1
        sig=pos-neg+answer
        return sig

    def copy(self):
        """
        Returns a copy of the knot.       
        >>> K=Link('4_1')
        >>> copy=K.copy()
        >>> K.PD_code() == copy.PD_code()
        True
        """
        return copy.deepcopy(self)

    def mirror(self):
        """
        Returns the mirror of a knot. Is not consistent about orientations.
        >>> K = Link('3_1')
        >>> K.crossings[0].sign                                                               
        1
        >>> mirror =K.mirror()                                                                   
        >>> mirror.crossings[0].sign                                                            
        -1
        """
        pd = self.PD_code()
        new_pd = list()
        for x in pd:
            new_pd.append (x[1:]+(x[0],))
        return Link(new_pd)

    def colorability_matrix(self):
        """Auxiliary function used by determinant. Returns 'colorability matrix'."""
        edges = self.pieces()
        m=matrix(len(self.crossings), len(edges))
        for c in self.crossings:
            for i in range(4):
                for s in edges:
                    if (c,i) in s:
                        ind=edges.index(s)
                        if i==1 or i==3:
                            m[(self.crossings.index(c),ind)]+=1
                        if i==2 or i==0:
                            m[(self.crossings.index(c),ind)]-=1
                        break
        return m

    def determinant(self, method='wirt'):
        """Returns the determinant of the knot K, a non-negative integer.                
       
        >>> K = Link( [(4,1,5,2),(6,4,7,3),(8,5,1,6),(2,8,3,7)] )  # Figure 8 knot
        >>> K.determinant()                                                        
        5
        """
        if not _within_sage:
            raise RuntimeError('determinant '+not_in_sage_msg)
        if method=='color':
            M = self.colorability_matrix()
            size = len(self.crossings)-1
            N = matrix(size, size)
            for i in range(size):
                for j in range(size):
                    N[(i,j)]=M[(i+1,j+1)]
            return abs(N.determinant())
        elif method=='goeritz':
            return abs(self.goeritz_matrix().determinant())
        else:
            t = var('t')
            return abs(self.alexander_poly(v=-1))
            


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

# ---- Drawing the link --------
#
#   This code is no longer used. 
#
#----------------------------

from . import draw
import tempfile 

def save_link_pdf(self, filename):
    """
    Requires Bartholomew's "draw programme",
    MetaPost, and the Ghostscript based "epstopdf". 
    """
    file = open(filename, 'wb')
    file.write(draw.link_pdf(self.peer_code()))
    file.close()

def show(self):
    # this leaves crud in your tmp directory
    file = tempfile.NamedTemporaryFile(delete=False)
    file.write(draw.link_pdf(self.peer_code()))
    file.close()
    if sys.platform == 'darwin':
        os.system('open ' + file.name)
    # this is not very portable.
    elif sys.platform == 'linux2':
        os.system('okular ' + file.name)

# Link.save_link_pdf, Link.show = save_link_pdf, show
