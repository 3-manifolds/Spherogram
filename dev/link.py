from spherogram import *

class Crossing:
    """
    See crossing.pdf for the convention, which matches
    that of KnotTheory.
    """
    def __init__(self, label=None, sign=0, adjacent=None):
        if adjacent == None:
            adjacent = [None, None, None, None]

        self.label, self.sign, self.adjacent = label, sign, adjacent
        self.strand_labels = [None, None, None, None]
        
    def __getitem__(self, i):
        return (self, i%4)

    def slot_is_empty(self, i):
        return self.adjacent[i % 4] == None
    
    def __setitem__(self, i, other):
        o, j = other
        if not (self.slot_is_empty(i) and o.slot_is_empty(j)):
            raise ValueError('Gluing two strands in the same place')
        self.adjacent[i % 4] = other
        other[0].adjacent[other[1]] = (self, i)
        

    def __repr__(self):
        return self.label

    def info(self):
        def format_adjacent(a):
            return (a[0].label, a[1]) if a else None
        print "<%s : %s>" % (self.label, [format_adjacent(a) for a in self.adjacent])

class CrossingStrand(DirectedEdge):
    def __mul__(self, other):
        a, b = self
        c, d = other
        assert set([a,b,c,d]) == set([0,1,2,3])
        return 1 if (a < b and c > d) or (a > b and c < d) else -1

    def kind(self):
        return "under" if self.tail() % 2 == 0 else "over"

over_right =CrossingStrand(3, 1)
over_left = CrossingStrand(1, 3)
under_right = CrossingStrand(2, 0)
under_left = CrossingStrand(0, 2)
strand_from_endpoint = (under_right, over_right, under_left, over_left)

class Link(Digraph):
    def __init__(self, crossings, check_planarity=True):
        if True in [ None in c.adjacent for c in crossings]:
            raise ValueError("No loose strands allowed")

        Digraph.__init__(self, [], [])

        remaining = set( [ (c, i) for c in crossings for i in range(4)] )
        components = []
        while len(remaining):
            start = remaining.pop()
            c, i = start
            component = []
            while True:
                component.append((c, strand_from_endpoint[i]))
                remaining.discard( (c,i) )
                d, j = c.adjacent[i]
                remaining.remove( (d,j) )
                self.add_edge(c, d)
                c, i = d, (j + 2) % 4
                if (c, i) == start:
                    break

            components.append(component)

        self.link_components = components
        if check_planarity and not self.is_planar():
            raise ValueError("Link isn't planar")

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
        for n, component in enumerate(components):
            for i, (c, s) in enumerate(component):
                c.strand_labels[s.tail()] = (n, i)
                c.strand_labels[s.head()] = (n, i + 1)

        component_lens = [len(component) for component in components]

        def strand_label( (n, i) ):
            return sum(component_lens[:n]) + (i % component_lens[n]) + 1

        PD = []
        for c in self.vertices:
            labels = [strand_label(s) for s in c.strand_labels]
            if c.strand_labels[0] > c.strand_labels[2]:
                labels = labels[2:] + labels[:2]
            PD.append(labels)

        if KnotTheory:
            PD = "PD" + repr(PD).replace('[', 'X[')[1:]
        return PD


def figure8():
    a, b, c, d = [Crossing(x) for x in 'abcd']
    a[0] = c[1]
    a[1] = d[0]
    a[2] = b[1]
    a[3] = b[0]
    b[2] = d[3]
    b[3] = c[2]
    c[3] = d[2]
    c[0] = d[1]
    return Link([a,b,c,d])

def punct_torus():
    a = Crossing('a')
    a[0] = a[2]
    a[1] = a[3]
    return Link([a], check_planarity=False)

def whitehead():
    a, b, c, d, e =  crossings = [Crossing(x) for x in 'abcde']
    a[0] = b[3]
    a[1] = d[0]
    a[2] = d[3]
    a[3] = c[0]
    b[0] = c[3]
    b[1] = e[2]
    b[2] = e[1]
    c[1] = d[2]
    c[2] = e[3]
    d[1] = e[0]
    return Link(crossings)

K, W, T = figure8(), whitehead(), punct_torus()
print K.is_planar(), W.is_planar(), punct_torus().is_planar()
print K.PD_code(True)
print W.PD_code(True)
print W.link_components

