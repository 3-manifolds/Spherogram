"""
Orthogonal layouts for link diagrams

Finds a layout of the given link diagram where the strands are unions of
edges in the standard integer grid, following:

Tamassia, On embedding a graph in the grid with the minimum number of
bends. Siam J. Comput. 16 (1987) http://dx.doi.org/10.1137/0216030

and

Bridgeman et. al., Turn-Regularity and Planar Orthogonal Drawings
ftp://ftp.cs.brown.edu/pub/techreports/99/cs99-04.pdf

As all vertices (=crossings) of the underlying graph are 4-valent,
things simpify; the associated network N(P) has A_V empty and A_F has no
self-loops.
"""

import spherogram, snappy, collections, networkx, random, plink
from spherogram.links.links import CrossingStrand, CrossingEntryPoint, Strand


#----- Utility code ------

def element_map(partition):
    ans = dict()
    for P in partition:
        for x in P:
            ans[x] = P
    return ans

def appearances_dict(list_of_lists):
    ans = dict()
    for L in list_of_lists:
        for x in L:
            if not ans.has_key(x):
                ans[x] = [L]
            else:
                ans[x].append(L)

    return ans

def partial_sums(L):
    ans, sum = [], 0
    for x in L:
        sum += x
        ans.append(sum)
    return ans

def topological_numbering(G):
    """
    Finds an optimal weighted topological numbering a directed acyclic graph
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

#----- Faces of link diagrams ------

class DummyVertex:
    def __repr__(self):
        return "<DummyVertex>"

def _subdivide_edge(self, crossing_strand, n):
    """
    Given a CrossingStrand, subdivides the edge for which it's the *head*
    into (n + 1) pieces.
    
    WARNING: this breaks several of the link's internal data structures.
    """
    head = crossing_strand
    backwards = not (head in head.crossing.entry_points())
    if backwards:
        head = head.opposite()
    tail = head.opposite()
    strands = [Strand() for i in range(n)]
    strands[0][0] = tail.crossing[tail.entry_point]
    for i in range(n - 1):
        strands[i][1] = strands[i+1][0]
    strands[-1][1] = head.crossing[head.entry_point]
        
spherogram.Link._subdivide_edge = _subdivide_edge
    
class Face(list):
    def __init__(self, link, crossing_strands, exterior=False):
        list.__init__(self, crossing_strands)
        self.edges = {link.CS_to_edge[c]:c for c in crossing_strands}
        self.exterior = exterior
        self.angles =[90 for e in self] 

    def edge_of_intersection(self, other):
        """
        Returns edge as the pair of CrossingStrands for self and other which
        defines the common edge. Returns None if the faces are disjoint.
        """
        common_edges = set(self.edges) & set(other.edges)
        if common_edges:
            e = common_edges.pop()
            return (self.edges[e], other.edges[e])
            
    def source_capacity(self):
        if self.exterior:
            return 0
        return max(4 - len(self), 0)

    def sink_capacity(self):
        if self.exterior:
            return len(self) + 4
        return max(len(self) - 4, 0)

    def iterate_from(self, edge):
        i = self.index(edge)
        return zip(self[i:] + self[:i], self.angles[i:] + self.angles[:i])

    def bend(self, edge, angles):
        """
        Assumes edge has already been subdivided by adding in len(angles) new
        vertices. 
        """
        i = self.index(edge)
        for a in reversed(angles):
            edge = edge.previous_corner()
            self.insert(i, edge)
            self.angles.insert(i, a)

    def orient_edges(self, edge, orientation):
        """
        Returns a dictionary of the induced directions of the edges (and their
        opposites) around the face.
        """
        dirs = ["left", "up", "right", "down"]
        dir = dirs.index(orientation)
        ans = dict()
        for e, a in self.iterate_from(edge):
            ans[e] = dirs[dir]
            ans[e.opposite()] = dirs[(dir + 2)%4]
            dir = (dir + a // 90) % 4
        return ans
    
    def __repr__(self):
        ext = '*' if self.exterior else ''
        return list.__repr__(self) + ext


class Faces(list):
    def __init__(self, link):
        self.link = link = link.copy()
        list.__init__(self, [Face(link, F) for F in link.faces()])
        F = max(self, key=len)
        F.exterior = True
        self.face_network = self.flow_networkx()
        self.bend()
        self.orient_edges()
        self.edges = sum([F for F in self], [])
        strands = {e.crossing for e in self.edges
                      if isinstance(e.crossing, Strand)}
        self.strand_CEPs =  [CrossingEntryPoint(s, 0) for s in strands]
        
    
    def flow_networkx(faces):
        """
        Tamassia's associated graph N(P) where the flow problem resides.
        """
        G = networkx.DiGraph()
        
        # Source
        source_demand = sum(F.source_capacity() for F in faces)
        G.add_node('s', demand = -source_demand)
        for i, F in enumerate(faces):
            if F.source_capacity():
                G.add_edge('s', i, weight=0, capacity=F.source_capacity())

        # Sink
        sink_demand = sum(F.sink_capacity() for F in faces)
        assert sink_demand == source_demand
        G.add_node('t', demand = sink_demand)
        for i, F in enumerate(faces):
            if F.sink_capacity():
                G.add_edge(i, 't', weight=0, capacity=F.sink_capacity())

        # Rest of edges
        for A in faces:
            for B in faces:
                if A != B and A.edge_of_intersection(B):
                    G.add_edge(faces.index(A), faces.index(B), weight=1)  # infinite capacity

        return G

    def bend(self):
        N = self.face_network
        flow = networkx.min_cost_flow(N)
        for a, flows in flow.iteritems():
            for b, w_a in flows.iteritems():
                if w_a and {'s', 't'}.isdisjoint({a, b}):
                    w_b = flow[b][a]
                    A, B = self[a], self[b]
                    e_a, e_b = A.edge_of_intersection(B)
                    angles_a = w_a*[90] + w_b*[270]
                    angles_b = w_b*[90] + w_a*[270]
                    self.link._subdivide_edge(e_a, len(angles_a))
                    A.bend(e_a, angles_a)
                    B.bend(e_b, angles_b)

    def orient_edges(self):
        """
        For each edge in a face, assign it one of four possible orientations:
        "left", "right", "up", "down".
        """
        orientations = {self[0][0]:'right'}
        N = self.face_network.to_undirected()
        N.remove_node('s'), N.remove_node('t')
        for i in networkx.traversal.dfs_preorder_nodes(N, 0):
            F = self[i]
            for edge in F:
                if orientations.has_key(edge):
                    new_orientations = F.orient_edges(edge, orientations[edge])
                    for e, dir in new_orientations.iteritems():
                        if orientations.has_key(e):
                            assert orientations[e] == dir
                        else:
                            orientations[e] = dir
                    break

        assert len(orientations) == sum(len(F) for F in self)
        self.orientations = orientations

    def orthogonal_rep(self):
        orientations = self.orientations
        spec = [[ (e.crossing, e.opposite().crossing) for e in self.edges if orientations[e] == dir] for dir in ['right', 'up']]
        ans = OrthogonalRep(*spec)

        face_sizes = sorted(len(F) for F in self)
        new_face_sizes = sorted(len(F) for F in ans.faces)
        assert face_sizes == new_face_sizes
        return ans


                
    def break_into_arrows(self):
        arrows = []
        for s in self.strand_CEPs:
            arrow = [s, s.next()]
            while not isinstance(arrow[-1].crossing, Strand):
                arrow.append(arrow[-1].next())
            arrows.append(arrow)

        undercrossings = dict()
        for i, arrow in enumerate(arrows):
            for a in arrow[1:-1]:
                if a.is_under_crossing():
                    undercrossings[a] = i

        crossings = []    # (under arrow, over arrow)
        for i, arrow in enumerate(arrows):
            for a in arrow[1:-1]:
                if a.is_over_crossing():
                    crossings.append( (undercrossings[a.other()], i) )
            
        return arrows, crossings

    def plink_data(self, spacing = 25):
        """
        Returns:
        * a list of vertex positions
        * a list of arrows join vertices
        * a list of crossings in the format (arrow over, arrow under)
        """
        emb = self.orthogonal_rep().basic_grid_embedding()
        vertex_positions = []
        for v in self.strand_CEPs:
            a, b = emb[v.crossing]
            vertex_positions.append( (spacing*(a+1), spacing*(b+1)) )

        vert_indices = {v:i for i, v in enumerate(self.strand_CEPs)}
        arrows, crossings = self.break_into_arrows()
        arrows = [ (vert_indices[a[0]], vert_indices[a[-1]]) for a in arrows]
        return vertex_positions, arrows, crossings
        

#----- Orhogonal representations ------
   
class OrthogonalFace(list):
    """
    A face of an OrthogonalRep oriented *clockwise* and stored as a list of
    pairs (edge, vertex) where vertex is gives the *clockwise* orientation
    of edge.
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
    An orthogonal representation is an equivalence class of planar
    embeddings of a graph where all edges are either vertical or
    horizontal. We assume there are no degree 1 vertices.
    
    Horizontal edges are oriented to the right, and vertical edges
    oriented upwards.
    
    >>> square = OrthogonalRep([(0, 1), (3, 2)], [(0, 3), (1,2)])
    """
    def __init__(self, horizonal_pairs=[], vertical_pairs=[]):
        spherogram.Digraph.__init__(self)
        for (a,b) in horizonal_pairs:
            self.add_edge(a, b, 'horizontal')
        for (a,b) in vertical_pairs:
            self.add_edge(a, b, 'vertical')

        self._build_faces()
        self._make_turn_regular()

    def add_edge(self, a, b, kind):
        e = spherogram.Digraph.add_edge(self, a, b)
        e.kind = kind
        return e

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
                e = self.add_edge(v0, v1, kind)
            else:
                e = self.add_edge(v1, v0, kind)
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

    def basic_grid_embedding(self, rotate=False):
        """
        Returns the positions of vertices under the grid embedding. 
        """
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
        
#----- plink code --------

from plink import Vertex, Arrow, Crossing

def load_from_spherogram(self, link, spacing=25):
    vertices, arrows, crossings = Faces(link).plink_data(spacing)
    self.clear()
    self.clear_text()
    for (x, y) in vertices:
        self.Vertices.append(Vertex(x, y, self.canvas))
    for u, v in arrows:
        U, V = self.Vertices[u], self.Vertices[v]
        self.Arrows.append(Arrow(U, V, self.canvas))
    for u, o in crossings:
        U, O = self.Arrows[u], self.Arrows[o]
        self.Crossings.append(Crossing(O, U))
    
    self.create_colors()
    self.goto_start_state()

plink.LinkEditor.load_from_spherogram = load_from_spherogram

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

def element(S):
    return list(S)[0]
    
    
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

cs = CrossingStrand(element(trefoil.vertices), 0)

BOR = OrthogonalRep([(0, 1), (2, 3), (3, 4), (5,6), (6,7), (8, 9)],
                    [(0,3), (2,5), (3, 6), (4, 7), (6, 8), (1, 9)])
faces = Faces(trefoil)
