def vertex_to_KLP(c, v):
    i = v if c.sign == 1 else (v - 1) % 4
    return ['Ybackward', 'Xforward', 'Yforward', 'Xbackward'][i]

class KLPCrossing():
    """
    SnapPea uses a convention where the orientation
    of the strands is fixed in the master picture but
    which strand is on top varies.
    """
    def __init__(self, c):
        self.adjacent = adjacent = 4*[None]
        self.index = c._KLP_index
        if c.sign == 1:
            strands, self.sign = ['over','under'], 'R'
        else:
            strands, self.sign = ['under','over'], 'L'
            
        components = [ c.strand_labels[s][0] for s in strands]
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
    vertices = list(L.vertices)
    for i, v in enumerate(vertices):
        v._KLP_index = i

    return [len(vertices), 0, len(L.link_components), [KLPCrossing(c) for c in vertices]]

try:
    import snappy
    def link_to_complement(L):
        P = python_KLP(L)
        return snappy.SnapPy.triangulate_link_complement_from_data(P)

except ImportError:
    def link_to_complement(L):
        raise RuntimeError("SnapPy doesn't seem to be available")
            

            
