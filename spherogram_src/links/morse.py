"""
The *Morse number* of a planar link diagram D is

m(D) = min { # of maxima of h on D }

where h is a height function on R^2 which is generic on D; alternatively,
this is the minimum number of cups/caps in a "MorseLink" presentation:

   http://katlas.math.toronto.edu/wiki/MorseLink_Presentations

of the diagram D.  The Morse number is very closely related to the more
traditional bridge number.  In this submodule, we compute m(D) via
integer linear programming, using ideas from

[DP] Didimo and Pizzonia, Upward Embeddings and Orientations of
Undirected Planar Graphs. <http://jgaa.info/getPaper?id=68>

To each corner of a face of D is classified into "large", "flat", or
"small" (which is the default) as per [DP], which correspond to integers
2, 1, and 0 respectively.  The geometric conditions are

* Every vertex has either a single "large" angle or two "flat" ones which
are opposite.

* The sum of the corner types around every face is degree - 2, with the
exception of the exterior face (which is degree + 2).

This code written by Nathan Dunfield who claims (not very plausibly) that
he will write a paper about this algorithm at some point.  The bit that
annoys him is that [DP] is polynomial time, but the algorithm here is not
known to be.  The issue is that [DP] creates a very special kind of ILP
(a network flow) which can be solved in polynomial time, but below we're
reduced to using a generic ILP solver.  
"""
from future.utils import iteritems
from ..sage_helper import _within_sage
from ..graphs import CyclicList, Digraph
from .links import CrossingStrand, Crossing, Strand, Link
from .orthogonal import basic_topological_numbering
from .tangles import RationalTangle
if _within_sage:
    from sage.numerical.mip import MixedIntegerLinearProgram

def morse_via_LP(link, solver='GLPK'):
    """
    An integer linear program which computes the Morse number of the given
    link diagram.

    sage: K = RationalTangle(23, 43).denominator_closure()
    sage: morse, details = morse_via_LP(K)
    sage: morse
    2
    sage: morse_via_LP(Link('8_20'))[0]
    3
    """
    LP = MixedIntegerLinearProgram(maximization=False, solver=solver)

    # Near a crossing, the level sets of the height function are either
    # horizontal or vertical, and so there are two flat corners which
    # must be opposite.  
    hor_cross = LP.new_variable(binary=True)
    vert_cross = LP.new_variable(binary=True)

    # We add a dummy vertex in the middle of each edge so it can bend.
    # The corners here must either both be flat or just one is large.  For
    # the corner to large, *both* flat_edge and large_edge must be 1.  
    flat_edge = LP.new_variable(binary=True)
    large_edge = LP.new_variable(binary=True)

    # Exactly one complementary region is the exterior one.
    exterior = LP.new_variable(binary=True)
    faces = link.faces()
    LP.add_constraint(sum(exterior[i] for i, F in enumerate(faces)) == 1)

    # Now for the main constraints 
    for c in link.crossings:
        LP.add_constraint(hor_cross[c] + vert_cross[c] == 1)
        for ce in c.entry_points():
            s = CrossingStrand(c, ce.strand_index)
            t = s.opposite()
            LP.add_constraint( flat_edge[s] == flat_edge[t] )
            LP.add_constraint( flat_edge[s] + large_edge[s] + large_edge[t] == 1 )

    for i, face in enumerate(faces):
        eqn = 0
        for cs in face:
            flat = hor_cross if cs.strand_index % 2 == 0 else vert_cross
            eqn += flat[cs.crossing] + flat_edge[cs] + 2*large_edge[cs]
        LP.add_constraint(eqn == (2*len(face) - 2 + 4*exterior[i]))

    LP.set_objective(sum(large_edge.values()))
    morse = int(LP.solve())
    assert morse % 2 == 0
    return morse//2, LP.get_values([hor_cross, vert_cross, flat_edge, large_edge, exterior])


def have_positive_value(D):
    return [k for k, v in iteritems(D) if v > 0]

class ImmutableValueDict(dict):
    def __setitem__(self, index, value):
        if index in self:
            if self[index] != value:
                raise ValueError("Can't change an assigned value")
        else:
            dict.__setitem__(self, index, value)

def pairing_to_permuation(pairing):
    points = sorted(sum(pairing, tuple()))
    assert points == range(len(points))
    ans = len(points)*[None]
    for x, y in pairing:
        ans[x], ans[y] = y, x
    return ans

class UpwardSnake(tuple):
    """
    Start with an MorseLinkDiagram, resolve all the crossings vertically,
    and snip all the mins/maxes.  The resulting pieces are the UpwardSnakes.
    """
    def __new__(self, crossing_strand, link):
        cs, kind = crossing_strand, link.orientations
        assert kind[cs] is 'min'
        snake = [cs]
        while True:
            ca, cb = cs.rotate(1).opposite(), cs.rotate(-1).opposite()
            if kind[ca] is 'max' or kind[cb] is 'max':
                break
            if kind[ca] is 'up':
                cs = ca
            else:
                assert kind[cb] is 'up'
                cs = cb
            snake.append(cs)

        ans = tuple.__new__(UpwardSnake, snake)
        ans.final = ca if kind[ca] is 'max' else cb

        heights = [link.heights[cs.crossing] for cs in ans]
        assert heights == sorted(heights)
        #assert heights[-1] == link.heights[ans.final.crossing]
        ans.heights = heights
        return ans
        
class MorseLinkDiagram(object):
    """
    A planar link diagram with a height function on R^2 which
    is Morse on the link. 
    """
    def __init__(self, link, solver='GLPK'):
        self.link = link = link.copy()
        morse, values = morse_via_LP(link, solver)
        self.morse_number = morse
        self.bends = set(have_positive_value(values[3]))
        self.faces = faces = link.faces()
        self.exterior = self.faces[have_positive_value(values[4])[0]]
        for c in have_positive_value(values[0]):
            c.kind = 'horizontal'
        for c in have_positive_value(values[1]):
            c.kind = 'vertical'
        self.orient_edges()
        self.set_heights()
        self.upward_snakes()

    def orient_edges(self):
        """
        Orients the edges of the link (that is, its CrossingStrands) with
        respect to the height function.

        sage: L = Link('K3a1')
        sage: D = MorseLinkDiagram(L)
        sage: orients = D.orientations.values()
        sage: len(orients) == 4*len(L.crossings)
        True
        sage: sorted(orients)
        ['down', 'down', 'max', 'max', 'max', 'max', 'min', 'min', 'min', 'min', 'up', 'up']
        sage: orients.count('max') == 2*D.morse_number
        True
        """
        def expand_orientation(cs, kind):
            c, i = cs
            kinds = CyclicList(['up', 'down', 'down', 'up'])
            if c.kind in 'horizontal':
                s = 0 if i in [0, 3] else 2
            elif c.kind == 'vertical':
                s = 1 if i in [2, 3] else 3
            if kind in ['down', 'max']:
                s += 2
            return [ (CrossingStrand(c,i), kinds[i+s]) for i in range(4) ]

        orientations = ImmutableValueDict()
        cs = list(self.bends)[0]
        co = cs.opposite()
        orientations[cs] = 'max'
        orientations[co] = 'max'
        current = [cs, co]
        while len(current):
            new = []
            for cs in current:
                for cn, kind in expand_orientation(cs, orientations[cs]):
                    co = cn.opposite()
                    if cn in self.bends or co in self.bends:
                        kind = {'up':'min', 'down':'max'}[kind]
                    if co not in orientations:
                        new.append(co)
                    orientations[cn] = kind
                    orientations[co] = {'up':'down', 'down':'up', 'max':'max', 'min':'min'}[kind]

            current = new
                    
        self.orientations = orientations
        
    def strands_below(self, crossing):
        """
        The two upward strands below the crossing.  
        """
        kinds = self.orientations
        a = CrossingStrand(crossing, 0)
        b = a.rotate()
        while True:
            if set([kinds[a], kinds[b]]).issubset(set(['up', 'min'])):
                return a, b
            a, b = b, b.rotate()

    def adjacent_upwards(self, crossing_strand):
        a, b = crossing_strand.rotate(), crossing_strand.rotate(-1)
        if self.orientations[a] in ['up', 'min']:
            return a
        else:
            assert self.orientations[b] in ['up', 'min']
            return b
                    
    def digraph(self):
        """
        The directed graph whose vertices are the mins/maxes of the height
        function together with the crossings, and where the edges come from
        the link and are directed upwards with respect to the height function.

        sage: L = Link('K4a1')
        sage: D = MorseLinkDiagram(L)
        sage: G = D.digraph()
        sage: len(G.vertices)
        8
        """
        G = Digraph()
        kinds = self.orientations
        for cs in self.bends:
            c, d = cs.crossing, cs.opposite().crossing
            if kinds[cs] == 'min':
                G.add_edge(cs, c), G.add_edge(cs, d)
            elif kinds[cs] == 'max':
                G.add_edge(c, cs), G.add_edge(d, cs)
                
        for cs, kind in iteritems(kinds):
            if kind == 'up':
                c, d  = cs.crossing, cs.opposite().crossing
                G.add_edge(d, c)

        return G

    def set_heights(self):
        """
        Assigns a height to each min/max and crossing of the
        diagram.
        """
        D = self.digraph()
        self.heights = basic_topological_numbering(D)

    def upward_snakes(self):
        """
        Resolve all the crossings vertically and snip all the mins/maxes. The
        resulting pieces are the UpwardSnakes.  For a diagram in bridge position,
        the number of snakes is just twice the bridge number.

        sage: D = MorseLinkDiagram(Link('8a1'))
        sage: len(D.snakes)
        4
        """
        kinds = self.orientations
        self.snakes = snakes = []
        for cs in self.bends:
            if kinds[cs] is 'min':
                snakes += [UpwardSnake(cs, self), UpwardSnake(cs.opposite(), self)]

        self.strand_to_snake = dict()
        for snake in snakes:
            for s in snake:
                self.strand_to_snake[s] = snake
            self.strand_to_snake[snake.final] = snake

        self.pack_snakes()
        

    def pack_snakes(self):
        """
        Give the snakes horizonal positions.
        """
        snakes, to_snake = self.snakes, self.strand_to_snake
        S = Digraph(singles=snakes)
        for c in self.link.crossings:
            a, b = self.strands_below(c)
            S.add_edge(to_snake[a], to_snake[b])

        for b in self.bends:
            a = b.opposite()
            S.add_edge(to_snake[a], to_snake[b])
            
        snake_pos = basic_topological_numbering(S)
        self.S, self.snake_pos = S, snake_pos
        heights = self.heights
        max_height = max(heights.values())
        snakes_at_height = dict()
        for h in range(0, max_height+1):
            at_this_height = []
            for snake in snakes:
                if heights[snake[0].crossing] <= h <= heights[snake[-1].crossing]:
                    at_this_height.append(snake)

            at_this_height.sort(key=lambda s:snake_pos[s])
            for i, s in enumerate(at_this_height):
                snakes_at_height[ (s, h) ] = i

        self.snakes_at_height = snakes_at_height

    def is_bridge(self):
        """
        Returns whether the link is in bridge position with respect to this
        height function.
        """
        return sorted(self.snake_pos.values()) == range(len(self.snakes))

    def bridge(self):
        if not self.is_bridge():
            raise ValueError("Morse function doesn't give a bridge diagram")

        def to_index(cs):
            return self.snake_pos[self.strand_to_snake[cs]]

        cross_data = []
        for c in self.link.crossings:
            a, b = self.strands_below(c)
            i, j = to_index(a), to_index(b)
            assert i < j
            cross = (i, j) if a.strand_index % 2 == 1 else (j, i)
            cross_data.append( (self.heights[c], cross) )
        cross_data.sort()
        
        def bottom_pairing(snake):
            cs = snake[0]
            return tuple(sorted([to_index(cs), to_index(cs.opposite())]))

        bottom = set(bottom_pairing(snake) for snake in self.snakes)

        def top_pairing(snake):
            cs = snake[-1]
            cn = self.adjacent_upwards(snake.final)
            return tuple(sorted([to_index(cs), to_index(cn)]))

        top = set(top_pairing(snake) for snake in self.snakes)
            
        return BridgeDiagram(bottom, [cd[1] for cd in cross_data], top)
                

class BridgeDiagram(object):
    """
    A proper bridge diagram of a link, that is, a height function
    where all the mins are below all the maxes.  
    """
    def __init__(self, bottom, crossings, top):
        self.bottom, self.crossings, self.top = bottom, crossings, top
        self.width = 2*len(bottom)
        self.name = 'None'
        
    def link(self):
        crossings = []
        curr_endpoints = self.width*[None]

        for x, y in self.bottom:
            s = Strand()
            crossings.append(s)
            curr_endpoints[x] = (s, 0)
            curr_endpoints[y] = (s, 1)

        for a, b in self.crossings:
            c = Crossing()
            crossings.append(c)
            if a < b:
                ins, outs = (3, 0), (2,1)
            else:
                ins, outs = (0, 1), (3,2)
                a, b = b, a
                        
            c[ins[0]] = curr_endpoints[a]
            c[ins[1]] = curr_endpoints[b]
            curr_endpoints[a] = (c, outs[0])
            curr_endpoints[b] = (c, outs[1])

        for x, y in self.top:
            join_strands(curr_endpoints[x], curr_endpoints[y])

        return Link(crossings)

    def bohua_code(self):
        b = self.width//2
        ans = [b]
        ans += pairing_to_permuation(self.bottom)
        ans += [len(self.crossings)] + list(sum(self.crossings, tuple()))
        ans += pairing_to_permuation(self.top)
        return self.name + '\t' + ' '.join(repr(a) for a in ans)

    def HF(self):
        import bohua_HF
        return bohua_HF.compute_HF(self.bohua_code())

    def is_proper(self):
        return max( abs(a-b) for a, b in self.crossings) < 2
