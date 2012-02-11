from spherogram.graphs import Graph, Digraph, DirectedEdge
import link_exterior, copy

class Crossing:
    """
    See crossings.pdf for the conventions.  The sign of the
    crossing can be in {0, +1, -1}.  In the first case, the
    strands at the crossings can have any orientations.
    """
    def __init__(self, label=None):
        self.label, self.sign, self.directions = label, 0, set()
        self.adjacent = [None, None, None, None]
        self.strand_labels = {'over':None, 'under':None}

    def make_tail(self, a):
        b = (a, (a + 2) % 4)
        self.directions.add(b)

    def rotate(self, s):
        """
        Rotate the incoming connections by 90*s degrees anticlockwise.  
        """
        rotate = lambda v : (v + s) % 4
        self.adjacent = self.adjacent[s:] + self.adjacent[:s]
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

    def slot_is_empty(self, i):
        return self.adjacent[i % 4] == None

    def entry_vertex(self, s):
        if s == 'under':
            return 0
        else:
            return 3 if self.sign == 1 else 1
    
    def __setitem__(self, i, other):
        o, j = other
        self.adjacent[i % 4] = other
        other[0].adjacent[other[1]] = (self, i)
        

    def __repr__(self):
        return self.label

    def info(self):
        def format_adjacent(a):
            return (a[0].label, a[1]) if a else None
        print "<%s : %s : %s : %s : %s>" % (self.label, self.sign, [format_adjacent(a) for a in self.adjacent], self.directions, self.strand_labels)

class Strand:
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
        return self.label

    def info(self):
        def format_adjacent(a):
            return (a[0].label, a[1]) if a else None
        print "<%s : %s>" % (self.label, [format_adjacent(a) for a in self.adjacent])
    

class Link(Digraph):
    def __init__(self, crossings, check_planarity=True):
        if True in [ None in c.adjacent for c in crossings]:
            raise ValueError("No loose strands allowed")

        # Fuse the strands.  If there any components made up
        # only of strands, these thrown out here.
        [s.fuse() for s in crossings if isinstance(s, Strand)]
        self.crossings = [c for c in crossings if not isinstance(c, Strand)]
        Digraph.__init__(self, [], [])
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

    def _build_components(self):
        """
        Each component is stored as a list of *entry points*
        to crossings.
        """
        remaining = set( [ (c, s) for c in self.crossings for s in ('over', 'under')] )
        components = []
        while len(remaining):
            component, n = [], 0
            c, s = remaining.pop()
            i = c.entry_vertex(s)
            start, finished = (c, i), False
            while not finished:
                component.append( (c, i) )
                remaining.discard( (c, s) )
                c.strand_labels[s] = (len(components), n)
                d, i = c.adjacent[(i+2)%4]
                self.add_edge(c, d)
                s = 'under' if i == 0 else 'over'
                c, n = d, n+1
                finished =  (c, i) == start

            components.append(component)

        self.link_components = components

    def copy(self):
        return copy.deepcopy(self)
    
    def is_planar(self):
        if not self.is_connected():
            return False
        
        G = Graph()
        for c in self.vertices:
            for g in (0.5, 1.5, 2.5, 3.5):  # corners of complementary regions at c
                d, j = c.adjacent[ int(g+0.5) % 4 ]
                G.add_edge( (c, g), (d, j + 0.5) )

        euler = len(self.vertices) - len(self.edges) + len(G.components())
        return euler == 2

    def PD_code(self, KnotTheory=False):
        components = self.link_components
        comp_lens = [len(component) for component in components]

        def label( (n, i) ):
            return sum(comp_lens[:n]) + (i % comp_lens[n]) + 1
        def next_label( (n, i) ):
            return label( (n, i+1) )
        PD = []
        for c in self.vertices:
            over, under = c.strand_labels['over'], c.strand_labels['under']
            if c.sign == 1:
                PD.append( [label(under), next_label(over), next_label(under), label(over)] )
            else:
                PD.append([label(under), label(over), next_label(under), next_label(over)])

        if KnotTheory:
            PD = "PD" + repr(PD).replace('[', 'X[')[1:]
        return PD

Link.exterior = link_exterior.link_to_complement
