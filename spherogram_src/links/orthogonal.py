"""
Orthogonal layouts for link diagrams

Finds a layout of the given link diagram where the strands are unions of
edges in the standard integer grid, following:

Tamassia, On embedding a graph in the grid with the minimum number of
bends. Siam J. Comput. 16 (1987) http://dx.doi.org/10.1137/0216030

and

Bridgeman et. al., Turn-Regularity and Planar Orthogonal Drawings
ftp://ftp.cs.brown.edu/pub/techreports/99/cs99-04.pdf

A more concise summary of the algorithm is contained in 

Hashemi and Tahmasbi, A better heuristic for area-compaction of orthogonal 
representations.  http://dx.doi.org/10.1016/j.amc.2005.03.007

As all vertices (=crossings) of the underlying graph are 4-valent, things simpify; 
the associated network N(P) has A_V empty and A_F has no self-loops.
"""
from future.utils import iteritems, itervalues
import networkx, random, string
from .links import CrossingStrand, CrossingEntryPoint, Strand
from ..graphs import CyclicList, Graph, Digraph
from collections import *

#---------------------------------------------------
#
# Utility code
#
#---------------------------------------------------

def element_map(partition):
    ans = dict()
    for P in partition:
        for x in P:
            ans[x] = P
    return ans

def partial_sums(L):
    ans = [0]
    for x in L:
        ans.append(ans[-1] + x)
    return ans

def basic_topological_numbering(G):
    """
    Finds an optimal weighted topological numbering a directed acyclic graph
    """
    in_valences = dict( (v,G.indegree(v)) for v in G.vertices  )
    numbering = {}
    curr_sources = [v for v,i in iteritems(in_valences) if i == 0]
    curr_number = 0
    while len(in_valences):
        new_sources = []
        for v in curr_sources: 
            in_valences.pop(v)
            numbering[v] = curr_number
            for e in G.outgoing(v):
                w = e.head
                in_valences[w] -= 1
                if in_valences[w] == 0:
                    new_sources.append(w)
            curr_sources = new_sources
        curr_number += 1

    return numbering

def topological_numbering(G):
    """
    Finds an optimal weighted topological numbering a directed acyclic graph
    which doesn't have any local moves which decrease the lengths of the
    (non-dummy) edges.
    """
    n = basic_topological_numbering(G)
    success = True
    while success:
        success = False
        for v in G.vertices:
            below = len([e for e in G.incoming(v) if e.dummy == False])
            above = len([e for e in G.outgoing(v) if e.dummy == False])
            if above != below:
                if above > below:
                    new_pos = min( n[e.head] for e in G.outgoing(v) ) - 1
                else: 
                    new_pos = max( n[e.tail] for e in G.incoming(v) ) + 1
                if new_pos != n[v]:
                    n[v] = new_pos
                    success = True
    return n


def kitty_corner(turns):
    rotations = partial_sums(turns)
    reflex_corners = [i for i, t in enumerate(turns) if t == -1]
    for r0 in reflex_corners:
            for r1 in [r for r in reflex_corners if r > r0]:
                if rotations[r1] - rotations[r0] == 2:
                    return (r0, r1)

#---------------------------------------------------
#
# Orthogonal Representations 
#
#---------------------------------------------------

LabeledFaceVertex = namedtuple('LabeledFaceVertex', ['index', 'kind', 'turn'])

class OrthogonalFace(CyclicList):
    """
    A face of an OrthogonalRep oriented *clockwise* and stored as a list of
    pairs (edge, vertex) where vertex is gives the *clockwise* orientation
    of edge.
    """
    def __init__(self, graph, edge_and_vertex):
        edge, vertex = edge_and_vertex
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
            (e1, v1) = self[i + 1]
            if e0.kind == e1.kind:
                turns.append(0)
            else:
                t = (e0.tail == v0)^(e1.head == v0)^(e0.kind == 'horizontal')
                turns.append(-1 if t else 1)

        rotation = sum(turns)
        assert abs(rotation) == 4
        self.exterior = (rotation == -4)

    def kitty_corner(self):
        if not self.exterior:
            return kitty_corner(self.turns)

    def is_turn_regular(self):
        return self.kitty_corner() is None

    def switches(self, swap_hor_edges):
        """
        Returns a list of (index, source|sink, -1|1). 
        """
        def edge_to_endpoints(e):
            if swap_hor_edges and e.kind == 'horizontal':
                return e.head, e.tail
            return e.tail, e.head 
        
        ans = []
        for i, (e0, v0) in enumerate(self):
            t0, h0 = edge_to_endpoints(e0)
            t1, h1 = edge_to_endpoints(self[i+1][0])
            if t0 == t1 == v0:
                ans.append( LabeledFaceVertex(i, 'source', self.turns[i]) )
            elif h0 == h1 == v0: 
                ans.append( LabeledFaceVertex(i, 'sink', self.turns[i]) )


        return ans

    def saturation_edges(self, swap_hor_edges):
        def saturate_face( face_info ):
            # Normalize so it starts with a -1 turn, if any
            for i, a in enumerate(face_info):
                if a.turn == -1:
                    face_info = face_info[i:] + face_info[:i]
                    break
            for i in range(len(face_info) - 2):
                x, y, z = face_info[i:i+3]
                if x.turn == -1 and y.turn== z.turn == 1:
                    a,b = (x, z) if  x.kind == 'sink' else (z, x)
                    remaining = face_info[:i] + [LabeledFaceVertex(z.index, z.kind, 1)] + face_info[i+3:]
                    return [ (a.index, b.index)  ] + saturate_face(remaining)
            return []
        
        new_edges = saturate_face(self.switches(swap_hor_edges))
        return [ (self[a][1], self[b][1]) for a, b in new_edges]

    def __repr__(self):
        ext = '*' if self.exterior else ''
        return list.__repr__(self) + ext
            
def saturate_face( face_info ):
            # Normalize so it starts with a -1 turn, if any
            for i, a in enumerate(face_info):
                if a.turn == -1:
                    face_info = face_info[i:] + face_info[:i]
                    break
            for i in range(len(face_info) - 2):
                x, y, z = face_info[i:i+3]
                if x.turn == -1 and y.turn== z.turn == 1:
                    a,b = (x, z) if  x.kind == 'sink' else (z, x)
                    remaining = face_info[:i] + [LabeledFaceVertex(z.index, z.kind, 1)] + face_info[i+3:]
                    return [ (a.index, b.index)  ] + saturate_face(remaining)
            return []
        
class OrthogonalRep(Digraph):
    """
    An orthogonal representation is an equivalence class of planar
    embeddings of a graph where all edges are either vertical or
    horizontal. We assume there are no degree 1 vertices.
    
    Horizontal edges are oriented to the right, and vertical edges
    oriented upwards.
    
    >>> square = OrthogonalRep([(0, 1), (3, 2)], [(0, 3), (1,2)])
    """
    def __init__(self, horizonal_pairs=[], vertical_pairs=[]):
        Digraph.__init__(self)
        for (a,b) in horizonal_pairs:
            self.add_edge(a, b, 'horizontal')
        for (a,b) in vertical_pairs:
            self.add_edge(a, b, 'vertical')

        self._build_faces()
        self._make_turn_regular()

    def add_edge(self, a, b, kind):
        e = Digraph.add_edge(self, a, b)
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
        edge_sides = set([ (e, e.head) for e in self.edges] + [ (e, e.tail) for e in self.edges ])
        while len(edge_sides):
            es = edge_sides.pop()
            face = OrthogonalFace(self, es)
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
            if len([e for e in self.incoming(v0) if e.kind == kind]):
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

    def saturation_edges(self, swap_hor_edges):
        return sum( [face.saturation_edges(swap_hor_edges) for face in self.faces], [])
        
    def DAG_from_direction(self, kind):
        H = Digraph(
            pairs = [e for e in self.edges if e.kind == kind],
            singles = self.vertices)

        maximal_chains = H.weak_components()
        vertex_to_chain = element_map(maximal_chains)
        D = Digraph(singles=maximal_chains)
        for e in [e for e in self.edges if e.kind != kind]:
            d = D.add_edge(vertex_to_chain[e.tail],
                       vertex_to_chain[e.head])
            d.dummy = e in self.dummy

        for u, v in self.saturation_edges(False):
            d = D.add_edge(vertex_to_chain[u], vertex_to_chain[v])
            d.dummy = True
        for u, v in self.saturation_edges(True):
            if kind == 'vertical':
                u, v = v, u
            d = D.add_edge(vertex_to_chain[u], vertex_to_chain[v])
            d.dummy = True

        D.vertex_to_chain = vertex_to_chain
        return D

    def chain_coordinates(self, kind):
        D = self.DAG_from_direction(kind)
        chain_coors = topological_numbering(D)
        return dict( (v,chain_coors[D.vertex_to_chain[v]]) for v in self.vertices )

    def basic_grid_embedding(self, rotate=False):
        """
        Returns the positions of vertices under the grid embedding. 
        """
        V = self.chain_coordinates('horizontal')
        H = self.chain_coordinates('vertical')
        return dict( (v,(H[v], V[v])) for v in self.vertices)

    def show(self, unit=10, labels=True):
        from sage.all import circle, text, line,  Graphics
        pos = self.basic_grid_embedding()
        for v, (a,b) in iteritems(pos):
            pos[v] = (unit*a, unit*b)
        if not labels:
            verts = [ circle(p, 1, fill=True) for p in itervalues(pos) ]
        else:
            verts =  [ text(repr(v), p, fontsize=20, color='black') for v, p in iteritems(pos) ]
            verts += [ circle(p, 1.5, fill=False) for p in itervalues(pos) ]
        edges = [ line( [pos[e.tail], pos[e.head]] ) for e in
                  self.edges if not e in self.dummy]
        G = sum(verts + edges, Graphics())
        G.axes(False)
        return G

#---------------------------------------------------
#
# Orthogonal link diagrams
#
#---------------------------------------------------

class Face(CyclicList):
    def __init__(self, link, crossing_strands, exterior=False):
        list.__init__(self, crossing_strands)
        self.edges = dict( (c.oriented(),c) for c in crossing_strands)
        self.exterior = exterior
        self.turns =[1 for e in self] 
    
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
        return zip(self[i:] + self[:i], self.turns[i:] + self.turns[:i])

    def bend(self, edge, turns):
        """
        Assumes edge has already been subdivided by adding in len(turns) new
        vertices. 
        """
        i = self.index(edge)
        for t in reversed(turns):
            edge = edge.previous_corner()
            self.insert(i, edge)
            self.turns.insert(i, t)

    def orient_edges(self, edge, orientation):
        """
        Returns a dictionary of the induced directions of the edges (and their
        opposites) around the face.
        """
        dirs = CyclicList(["left", "up", "right", "down"])
        dir = dirs.index(orientation)
        ans = dict()
        for e, t in self.iterate_from(edge):
            ans[e] = dirs[dir]
            ans[e.opposite()] = dirs[dir + 2]
            dir = (dir + t) % 4
        return ans

    def is_turn_regular(self):
        if self.exterior:
            return True
        else:
            return kitty_corner(self.turns) is None
    
    def __repr__(self):
        ext = '*' if self.exterior else ''
        return list.__repr__(self) + ext

def subdivide_edge(crossing_strand, n):
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
    strands[0][0] = tail.crossing[tail.strand_index]
    for i in range(n - 1):
        strands[i][1] = strands[i+1][0]
    strands[-1][1] = head.crossing[head.strand_index]

        
class OrthogonalLinkDiagram(list):
    """
    A link diagram where all edges are made up of horizontal and vertical
    segments.
    """
    
    def __init__(self, link):
        self.link = link = link.copy()
        list.__init__(self, [Face(link, F) for F in link.faces()])
        F = max(self, key=len)
        F.exterior = True
        self.face_network = self.flow_networkx()
        self.bend()
        self.orient_edges()
        self.edges = sum([F for F in self], [])
        self.repair_components()
        
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
        """
        Computes a minimal size set of edge bends that allows the link diagram
        to be embedded orthogonally. This follows directly Tamassia's first
        paper.
        """
        N = self.face_network
        flow = networkx.min_cost_flow(N)
        for a, flows in iteritems(flow):
            for b, w_a in iteritems(flows):
                if w_a and set(['s', 't']).isdisjoint(set([a, b])):
                    w_b = flow[b][a]
                    A, B = self[a], self[b]
                    e_a, e_b = A.edge_of_intersection(B)
                    turns_a = w_a*[1] + w_b*[-1]
                    turns_b = w_b*[1] + w_a*[-1]
                    subdivide_edge(e_a, len(turns_a))
                    A.bend(e_a, turns_a)
                    B.bend(e_b, turns_b)

    def repair_components(self):
        # Repair the link components and store their numbering for
        # easy reference.  Also order the strand_CEPs so they go component
        # by component.
        self.link.link_components = [component[0].component()
                                     for component in self.link.link_components]
        self.strand_CEP_to_component = stc = dict()
        self.strand_CEPs = []
        for n, component in enumerate(self.link.link_components):
            for ce in component:
                if isinstance(ce.crossing, Strand):
                    stc[ce] = n
                    self.strand_CEPs.append(ce)


    def orient_edges(self):
        """
        For each edge in a face, assign it one of four possible orientations:
        "left", "right", "up", "down".
        """
        orientations = {self[0][0]:'right'}
        N = self.face_network.to_undirected()
        G = Graph(N.edges())
        G.remove_vertex('s'), G.remove_vertex('t')
        for i in G.depth_first_search(0):
            F = self[i]
            for edge in F:
                if edge in orientations:
                    new_orientations = F.orient_edges(edge, orientations[edge])
                    for e, dir in iteritems(new_orientations):
                        if e in orientations:
                            assert orientations[e] == dir
                        else:
                            orientations[e] = dir
                    break

        assert len(orientations) == sum(len(F) for F in self)
        self.orientations = orientations

    def orthogonal_rep(self):
        orientations = self.orientations
        spec = [[ (e.crossing, e.opposite().crossing) for e in self.edges if orientations[e] == dir] for dir in ['right', 'up']]
        return OrthogonalRep(*spec)

    def orthogonal_spec(self):
        orientations = self.orientations
        return [[ (e.crossing.label, e.opposite().crossing.label) for e in self.edges if orientations[e] == dir] for dir in ['right', 'up']]

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
                    crossings.append( (undercrossings[a.other()], i, False, a.crossing.label) )
            
        return arrows, crossings

    def plink_data(self):
        """
        Returns:
        * a list of vertex positions
        * a list of arrows joining vertices
        * a list of crossings in the format (arrow over, arrow under)
        """
        emb = self.orthogonal_rep().basic_grid_embedding()
        x_max = max(a for a,b in emb.values())
        y_max = max(b for a,b in emb.values())

        # We rotate things so the long direction is horizontal.  The
        # Plink canvas coordinate system forces us to flip things to
        # preserve crossing type.
        vertex_positions = []
        for v in self.strand_CEPs:
            if x_max >= y_max:
                a, b = emb[v.crossing]
                b = y_max - b
            else:
                b, a = emb[v.crossing]
            vertex_positions.append( (10*(a+1), 10*(b+1)) )

        vert_indices = dict( (v,i) for i, v in enumerate(self.strand_CEPs))
        arrows, crossings = self.break_into_arrows()
        arrows = [ (vert_indices[a[0]], vert_indices[a[-1]]) for a in arrows]

        return vertex_positions, arrows, crossings
        
#---------------------------------------------------
#
#  Opening in plink.  
#
#---------------------------------------------------

def orthogonal_draw(self, link_editor=None):
    """
    Opens a Plink link editor window with displaying the current link.
    The strands of the links are unions of edges in the standard
    integer grid, following the work of `Tamassia
    <http://dx.doi.org/10.1137/0216030>`_ and `Bridgeman
    et. al. <ftp://ftp.cs.brown.edu/pub/techreports/99/cs99-04.pdf>`_
    """
    if link_editor is None:
        import plink
        link_editor = plink.LinkEditor()
    diagram = OrthogonalLinkDiagram(self)
    link_editor.unpickle(*diagram.plink_data())
    try:
        link_editor.zoom_to_fit()
        link_editor.goto_start_state()
    except AttributeError:  # Have been passed just a LinkManager which has no Tk canvas
        pass
    return link_editor

#---------------------------------------------------
#
#  Testing code
#
#---------------------------------------------------



if __name__ == '__main__':
    import snappy
    import plink
    from spherogram import DTcodec
    
    def link_from_manifold(manifold):
        return DTcodec(manifold.DT_code()).link()

    def random_link():
        return link_from_manifold(HTLinkExteriors.random())

    def check_faces(link):
        faces = link.faces()
        assert len(link.vertices) - len(link.edges) + len(faces) == 2
        assert set(Counter(sum( faces, [] )).values()) == set([1])
        assert link.is_planar()

    def test_face_method(N):
        for i in xrange(N):
            check_faces(random_link())

    def element(S):
        return list(S)[0]

    def appears_hyperbolic(M):
        acceptable = ['all tetrahedra positively oriented',
                      'contains negatively oriented tetrahedra']
        return M.solution_type() in acceptable and M.volume() > 1.0 

    def test(manifold_with_DT, plink_manifold=None):
        L = link_from_manifold(manifold_with_DT)
        PM = plink_manifold
        if PM is None:
            PM = snappy.Manifold()
        PM.LE.load_from_spherogram(L, None, False)
        PM.LE.callback()
        if appears_hyperbolic(PM):
            assert abs(manifold_with_DT.volume() - PM.volume())  < 0.000001
            assert manifold_with_DT.is_isometric_to(PM)

    def big_test():
        PM = snappy.Manifold()
        while 1:
            M = snappy.HTLinkExteriors.random()
            print('Testing Manifold: ' + M.name())
            ans = test(M, PM)

    from .tangles import RationalTangle
    unknot = RationalTangle(1).numerator_closure()
    hopf = RationalTangle(2).numerator_closure()
    trefoil = DTcodec([(4,6,2)]).link()
    big_knot = DTcodec([(4, 12, 14, 22, 20, 2, 28, 24, 6, 10, 26, 16, 8, 18)]).link()
    big_link = DTcodec([(8, 12, 16), (18, 22, 24, 20), (4, 26, 14, 2, 10, 6)]).link()

    square = OrthogonalRep([(0, 1), (3, 2)], [(0, 3), (1,2)])
    OR = OrthogonalRep([ (0,1), (2, 3), (3, 4), (5, 6), (6, 7), (7, 8), (9, 10)],
                         [(0, 2), (1, 3), (2, 6), (3, 7), (4,8),(5,9),(6, 10)])
    kinked = OrthogonalRep([ (0, 1), (1, 2), (3, 4), (6, 7) ],
                           [(0,3),(4,6), (2,5), (5,7)])
    kinked2 = OrthogonalRep([ (0,1), (2,3), (4,5), (6,7)], [(1,2), (0, 4), (3,7), (5,6)])
    kinked3 = OrthogonalRep([ (0,1), (2,3), (4,5), (6,7)], [(1,2), (0, 8),
                                                            (8,4), (3,7), (5,6)])
    BOR = OrthogonalRep([(0, 1), (2, 3), (3, 4), (5,6), (6,7), (8, 9)],
                        [(0,3), (2,5), (3, 6), (4, 7), (6, 8), (1, 9)])

    BOR2 = OrthogonalRep( [(0,3), (2,5), (3, 6), (4, 7), (6, 8), (1, 9)], [(0, 1), (2, 3), (3, 4), (5,6), (6,7), (8, 9)])

    PDG = Digraph([(0, 4), (3, 0), (3, 0), (1, 2), (0, 3), (1, 3), (0, 4), (0, 4), (1, 0), (4, 2), (3,0)])

    FOR2 = OrthogonalRep(  [(0, 8), (1, 2), (3,4), (6,5), (7, 9)], [ (0,1), (2, 3), (4, 5), (6, 7), (8, 9) ])
    FOR1 = OrthogonalRep(  [ (0,1), (2, 3), (4, 5), (6, 7), (8, 9) ], [(0, 8), (1, 2), (3,4), (6,5), (7, 9)])
    BDG = Digraph([(1, 2), (5, 4), (1, 3), (1, 5), (2, 0), (0, 4), (5, 3), (5, 2), (3, 5)])
    BOR3 = OrthogonalRep( [ (0, 1), (2, 3), (4,5), (6,7)], [(0,2), (3,4), (5,6), (1, 7)] )
    BOR4 = OrthogonalRep( [ ('a', 'b'), ('c', 'd'), ('e', 'f'), ('g', 'h') ], [('a','g'), ('b','d'), ('c','f'), ('e', 'h')])
    BOR5 = OrthogonalRep([(8,9), (3,2), (1, 0), (5, 4), (7,6)],[ (7, 8), (6, 5), (4,3), (1, 2), (0, 9)] )
    face = OrthogonalFace(BOR5, (BOR5.outgoing(0).pop(), 0) )
    spec = [[(3, 5), (2, 'b'), ('j', 'd'), (5, 2), ('i', 'n'), ('h', 'f'), ('k', 0), ('m', 'g'), ('a', 3), (4, 'e'), (0, 1), (6, 4), ('l', 6), (1, 'c')], [('b', 'g'), (5, 'j'), (2, 'd'), ('n', 4), ('i', 6), ('a', 'h'), (0, 3), ('l', 'k'), (3, 'm'), ('e', 'f'), (1, 5), (4, 1), (6, 0), ('c', 2)]]

    #big_test()
