"""
Orthogonal layouts for link diagrams

Finds a layout of the given link diagram where the strands are unions of edges in the standard integer grid, following:

Tamassia, On embedding a graph in the grid with the minimum number of bends.
Siam J. Comput. 16 (1987)
http://dx.doi.org/10.1137/0216030

As all vertices (=crossings) of the underlying graph are 4-valent, things simpify;
the associated network N(P) has A_V empty and A_F has no self-loops.  
"""

import spherogram, snappy, collections, networkx
from spherogram.links.links import CrossingEntryPoint


    
class Face:
    def __init__(self, edges, exterior=False, bends=None, angles=None):
        self.edges, self.exterior = edges, exterior
        if bends is None:
            bends = {e:'' for e in edges }
        if angles is None:
            angles = {e:None for e in edges}
        self.bends, self.angles = bends, angles

    def intersection(self, other):
        return set(self.edges).intersection(set(other.edges))

    def source_capacity(self):
        if self.exterior:
            return 0
        return max(4 - len(self), 0)

    def sink_capacity(self):
        if self.exterior:
            return len(self) + 4
        return max(len(self) - 4, 0)

    def __len__(self):
        return len(self.edges)

def all_faces(link):
    ans = [Face(F) for F in link.faces()]
    F = max(ans, key=len)
    F.exterior = True
    return ans

def flow_network(link):
    """
    Tamassia's associated graph N(P) where the flow
    problem resides.  
    """
    faces = all_faces(link)
    vertices = {'s', 't'}.union(set(faces))
    G = spherogram.Digraph(singles=vertices)
    for A in faces:
        for B in faces:
            if A != B and A.intersection(B):
                G.add_edge(A, B)

    for F in faces:
        c = F.source_capacity()
        if c:
            G.add_edge('s', F)

        C = F.sink_capacity()
        if C:
            G.add_edge(F, 't')
    
    return G

def flow_networkx(link):
    """
    Tamassia's associated graph N(P) where the flow
    problem resides.  
    """
    faces = all_faces(link)
    G = networkx.DiGraph()
    
    # Source
    source_demand = sum(F.source_capacity() for F in faces)
    G.add_node('s', demand = -source_demand)
    for F in faces:
        if F.source_capacity():
            G.add_edge('s', faces.index(F), weight=0, capacity=F.source_capacity())

    # Sink
    sink_demand = sum(F.sink_capacity() for F in faces)
    assert sink_demand == source_demand
    G.add_node('t', demand = sink_demand)
    for F in faces:
        if F.sink_capacity():
            G.add_edge( faces.index(F), 't', weight=0, capacity=F.sink_capacity())

    # Rest of edges
    for F0 in faces:
        for F1 in faces:
            if F0 != F1:
                G.add_edge(faces.index(F0), faces.index(F1), weight=1)  # infinite capacity

    return G
    



def random_link():
    return spherogram.DTcodec(snappy.HTLinkExteriors.random().DT_code()).link()

def check_faces(link):
    faces = link.faces()
    assert len(link.vertices) - len(link.edges) + len(faces) == 2
    assert set(collections.Counter(sum( faces, [] )).values()) == {2}
    assert link.is_planar()

def test_face_method(N):
    for i in xrange(N):
        check_faces(random_link())
    
unknot = spherogram.RationalTangle(1).numerator_closure()
hopf = spherogram.RationalTangle(2).numerator_closure()
trefoil = spherogram.DTcodec([(4,6,2)]).link()
big_knot = spherogram.DTcodec([(4, 12, 14, 22, 20, 2, 28, 24, 6, 10, 26, 16, 8, 18)]).link()
big_link = spherogram.DTcodec([(8, 12, 16), (18, 22, 24, 20), (4, 26, 14, 2, 10, 6)]).link()
