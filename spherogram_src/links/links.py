from __future__ import print_function
"""
Links are made from Crossings.  The general model is that of
a PD diagram as described in 

  http://katlas.org/wiki/Planar_Diagrams

See the file "doc.pdf" for the conventions, and the file
"test.py" for some examples of creating links.

"""
from .. import graphs
CyclicList = graphs.CyclicList
import  string, os, sys, re, collections
try:
    import cPickle as pickle
except ImportError: # Python 3
    import pickle
    
class Crossing(object):
    """
    See crossings.pdf for the conventions.  The sign of the
    crossing can be in {0, +1, -1}.  In the first case, the
    strands at the crossings can have any orientations.
    """
    def __init__(self, label=None):
        self.label, self.sign, self.directions = label, 0, set()
        self.adjacent = CyclicList([None, None, None, None])
        self.strand_labels = CyclicList([None, None, None, None])
        self.strand_components = CyclicList([None, None, None, None])

    def make_tail(self, a):
        b = (a, (a + 2) % 4)
        self.directions.add(b)

    def rotate(self, s):
        """
        Rotate the incoming connections by 90*s degrees anticlockwise.  
        """
        rotate = lambda v : (v + s) % 4
        self.adjacent = CyclicList(self.adjacent[s:] + self.adjacent[:s])
        for i, (o, j) in enumerate(self.adjacent):
            if o != self:
                o.adjacent[j] = (self, i)
            else:
                self.adjacent[i] = (self, rotate(j))

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
        other[0].adjacent[other[1]] = (self, i)
        

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
        return [o for o in self.crossing.entry_points() if o != self][0]

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
    
    def __init__(self, crossings, check_planarity=True):
        # If crossings is just a PD code rather than a list of Crossings,
        # we create the corresponding Crossings.
        if len(crossings) > 0 and not isinstance(crossings[0], (Strand, Crossing)):
            gluings = collections.defaultdict(list)
            for c, X in enumerate(crossings):
                for i, x in enumerate(X):
                    gluings[x].append( (c,i) )

            if {len(v) for v in gluings.values()} != {2}:
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
        self._crossing_entries = set()
        self._orient_crossings()
        self._build_components()

        if check_planarity and not self.is_planar():
            raise ValueError("Link isn't planar")

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
                remaining = remaining - set( [ (c, i), (d, j) ])
                c, i = d, (j+2) % 4
                finished = (c, i) == start

        for c in self.crossings:
            c.orient()

    def crossing_entries(self):
        return sum([C.entry_points() for C in self.crossings], [])

    def crossing_strands(self):
        return sum([C.crossing_strands() for C in self.crossings], [])

    def _build_components(self):
        """
        Each component is stored as a list of *entry points* to
        crossings.  The labeling of the entry points (equivalently
        oriented edges) is compatible with the DT convention that each
        crossing has both an odd and and even incoming strand.
        """
        remaining, components = set( self.crossing_entries() ), LinkComponents()
        self.labels = labels = Labels()
        while len(remaining):
            if len(components) == 0:
                d = remaining.pop()
            else:
                found = False
                for c in sum(components, []):
                    d = c.other()
                    if d in remaining:
                        if labels[c]  % 2 == 0:
                            d = d.next()
                        found = True
                        break
                if not found:
                    d = remaining.pop()
        
            component = components.add(d)
            for c in component:
                labels.add( c )
            for c in component:
                c.label_crossing(len(components) - 1, labels)
            remaining = remaining - set(component)

        # Build the underlying graph
        self.digraph = graphs.Digraph()
        for component in components:
            for c in component:
                cs0 = CrossingStrand(c.crossing, c.entry_point)
                cs1 = cs0.opposite()
                e = self.digraph.add_edge(cs0.crossing, cs1.crossing)

        self.link_components = components

    def copy(self):
        return pickle.loads(pickle.dumps(self))
    
    def is_planar(self):
        if not self.digraph.is_weakly_connected():
            return False
        v = len(self.crossings)
        assert 2*v == len(self.digraph.edges)
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
        DT = [ [DT_dict[x] for x in component] for component in odd_labels]

        if DT_alpha:
            if len(self) > 52:
                raise ValueError("Too many crossing for alphabetic DT code")
            DT_alphabet = '_abcdefghijklmnopqrstuvwxyzZYXWVUTSRQPONMLKJIHGFEDCBA'
            init_data = [len(self), len(DT)] +  [len(c) for c in DT]
            DT = ''.join([DT_alphabet[x] for x in init_data] + [DT_alphabet[x>>1] for x in sum(DT, [])])
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

    def exterior(L):
        raise RuntimeError("SnapPy doesn't seem to be available.  Try: from snappy import *")

    def __repr__(self):
        return "<Link: %d comp; %d cross>" % (len(self.link_components), len(self.crossings))

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

Link.save_link_pdf, Link.show = save_link_pdf, show
