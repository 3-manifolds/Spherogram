"""
Orthogonal layouts for link diagrams

Finds a layout of the given link diagram where the strands are unions of edges in the standard integer grid, following:

Tamassia, On embedding a graph in the grid with the minimum number of bends.
Siam J. Comput. 16 (1987)
http://dx.doi.org/10.1137/0216030

and

Bridgeman et. al., Turn-Regularity and Planar Orthogonal Drawings
ftp://ftp.cs.brown.edu/pub/techreports/99/cs99-04.pdf

As all vertices (=crossings) of the underlying graph are 4-valent, things simpify;
the associated network N(P) has A_V empty and A_F has no self-loops.  
"""

import spherogram, snappy, collections, networkx, random
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

def element_map(partition):
    ans = dict()
    for P in partition:
        for x in P:
            ans[x] = P
    return ans

def partial_sums(L):
    ans, sum = [], 0
    for x in L:
        sum += x
        ans.append(sum)
    return ans
        
def topological_numbering(G):
    """
    Finds an optimal weighted topological
    numbering a directed acyclic graph
    """
    in_valences = { v:G.in_valence(v) for v in G.vertices }
    numbering = {}
    curr_sources = [v for v,i in in_valences.iteritems() if i == 0]
    curr_number = 0
    while len(in_valences):
        new_sources = []
        for v in curr_sources: 
            in_valences.pop(v)
            numbering[v] = curr_number
            for e in G.incident(v):
                w = e.head
                in_valences[w] -= 1
                if in_valences[w] == 0:
                    new_sources.append(w)
            curr_sources = new_sources
        curr_number += 1

    return numbering

   
class OrthogonalFace(list):
    """
    A face of an OrthogonalRep oriented clockwise
    and stored as a list of pairs (edge, vertex) where
    vertex is gives the clockwise orientation of edge.
    """
    def __init__(self, graph, (edge, vertex)):
        self.append( (edge, vertex) )
        while True:
            edge, vertex = self[-1]
            edge = graph.next_edge_at_vertex(edge, vertex)
            vertex = edge(vertex)
            if self[0] == (edge, vertex):
                break
            else:
                self.append( (edge, vertex) )
                
        self._add_turns()
        
    def _add_turns(self):
        self.turns = turns = []
        for i, (e0, v0) in enumerate(self):
            (e1, v1) = self[(i + 1) % len(self)]
            if e0.kind == e1.kind:
                turns.append(0)
            else:
                t = (e0.tail == v0)^(e1.head == v0)^(e0.kind == 'horizontal')
                turns.append(-1 if t else 1)

        rotation = sum(turns)
        assert abs(rotation) == 4
        self.exterior = (rotation == -4)

    def kitty_corner(self):
        if self.exterior:
            return None
        turns = self.turns
        rotations = partial_sums(turns)
        reflex_corners = [i for i, t in enumerate(turns) if t == -1]
        for r0 in reflex_corners:
            for r1 in reflex_corners:
                if rotations[r1] - rotations[r0] == 2:
                    return (r0, r1)

    def is_turn_regular(self):
        return self.kitty_corner() is None
            
        
class OrthogonalRep(spherogram.Digraph):
    """
    An orthogonal representation is an equivalence class
    of planar embeddings of a graph where all edges are
    either vertical or horizontal.   We assume there are
    no degree 1 vertices.  

    Horizontal edges are oriented to the right, and vertical
    edges oriented upwards.

    >>> square = OrthogonalRep([(0, 1), (3, 2)], [(0, 3), (1,2)])
    """
    def __init__(self, horizonal_pairs, vertical_pairs):
        spherogram.Digraph.__init__(self)
        for (a,b) in horizonal_pairs:
            e = self.add_edge(a, b)
            e.kind = 'horizontal'
        for (a,b) in vertical_pairs:
            e = self.add_edge(a, b)
            e.kind = 'vertical'

        self._build_faces()
        self._make_turn_regular()

    def link(self, vertex):
        """
        The edges in the link oriented counterclockwise.
        """
        link = self.incidence_dict[vertex]
        link.sort(key=lambda e: (e.tail == vertex, e.kind == 'vertical'))
        return link

    def next_edge_at_vertex(self, edge, vertex):
        """
        The next edge at the vertex
        """
        link = self.link(vertex)
        return link[ (link.index(edge) + 1) % len(link) ]
        
    def _build_faces(self):
        self.faces = []
        edge_sides = { (e, e.head) for e in self.edges } |  { (e, e.tail) for e in self.edges }
        while len(edge_sides):
            face = OrthogonalFace(self, edge_sides.pop())
            edge_sides.difference_update(face)
            self.faces.append(face)

    def _make_turn_regular(self):
        dummy = set()
        regular = [F for F in self.faces if F.is_turn_regular()]
        irregular =[F for F in self.faces if not F.is_turn_regular()]
        while len(irregular):
            F = irregular.pop()
            i,j = F.kitty_corner()
            v0, v1 = F[i][1], F[j][1]
            kind = random.choice( ('vertical', 'horizontal'))
            if len([e for e in self.incident_to(v0) if e.kind == kind]):
                e = self.add_edge(v0, v1)
            else:
                e = self.add_edge(v1, v0)
            e.kind = kind
            dummy.add(e)
            for v in [v0, v1]:
                F = OrthogonalFace(self, (e, v))
                if F.is_turn_regular():
                    regular.append(F)
                else:
                    irregular.append(F)

        self.faces, self.dummy = regular, dummy

    def DAG_from_direction(self, kind):
        H = spherogram.Digraph(
            pairs = [e for e in self.edges if e.kind == kind],
            singles = self.vertices)
        maximal_chains = H.components()
        vertex_to_chain = element_map(maximal_chains)
        D = spherogram.Digraph(singles=maximal_chains)
        for e in [e for e in self.edges if e.kind != kind]:
            D.add_edge(vertex_to_chain[e.tail],
                       vertex_to_chain[e.head])
        D.vertex_to_chain = vertex_to_chain
        return D

    def chain_coordinates(self, kind):
        D = self.DAG_from_direction(kind)
        chain_coors = topological_numbering(D)
        return { v:chain_coors[D.vertex_to_chain[v]] for v in self.vertices }

    def basic_grid_embedding(self):
        V = self.chain_coordinates('horizontal')
        H = self.chain_coordinates('vertical')
        return {v:(H[v], V[v]) for v in self.vertices}

    def show(self, unit=10):
        pos = self.basic_grid_embedding()
        for v, (a,b) in pos.iteritems():
            pos[v] = (unit*a, unit*b)
        verts = [ circle(p, 1, fill=True) for p in pos.itervalues() ]
        edges = [ line( [pos[e.tail], pos[e.head]] ) for e in
                  self.edges if not e in self.dummy]
        return sum(verts + edges, Graphics())
        

#----- testing code --------

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

square = OrthogonalRep([(0, 1), (3, 2)], [(0, 3), (1,2)])
OR = OrthogonalRep([ (0,1), (2, 3), (3, 4), (5, 6), (6, 7), (7, 8), (9, 10)],
                     [(0, 2), (1, 3), (2, 6), (3, 7), (4,8),(5,9),(6, 10)])
kinked = OrthogonalRep([ (0, 1), (1, 2), (3, 4), (6, 7) ],
                       [(0,3),(4,6), (2,5), (5,7)])
kinked2 = OrthogonalRep([ (0,1), (2,3), (4,5), (6,7)], [(1,2), (0, 4), (3,7), (5,6)])
kinked3 = OrthogonalRep([ (0,1), (2,3), (4,5), (6,7)], [(1,2), (0, 8),
                                                        (8,4), (3,7), (5,6)])
