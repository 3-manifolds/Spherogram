from spherogram import FatGraph

def sign(x):
    return 1 if x > 0 else -1 if x < 0 else 0

# To construct a knot projection from a DT code we first construct
# a fat graph, which may not be a planar surface.  There are only
# two possible orderings of edges at each vertex, since we know which
# pairs of edges are opposites.  Thus the process of finding the
# projection consists of reversing the orderings of some vertices until
# we get a planar surface.
#
# Each crossing in the diagram is traversed twice in building the DT
# code.  Our convention is to label the four edges at the crossing by
# 0, 1, 2, 3 so that the first component to pass through the vertex
# enters at 0 and leaves at 2, and the second component enters at 1 and
# leaves at 3.
        
class DTvertex:
    """
    A vertex of the 4-valent graph which is described by a DT code.
    Instantiate with an even-odd pair, in either order.
    """
    def __init__(self, pair, overcrossing=1):
        x, y = pair
        self.first = x if x < y else y
        self.second = y if x < y else x
        self.even = x if y%2 else y
        self.odd = y if y%2 else x
        self.even_over = True if overcrossing == -1 else False

    def __repr__(self):
        return str((self.first, self.second))

    def enter(self, N):
        if N == self.first: return 0
        elif N == self.second: return 1
        else: raise ValueError('%d is not a label of this vertex'%N)

    def exit(self, N):
        if N == self.first: return 2
        elif N == self.second: return 3
        else: raise ValueError('%d is not a label of this vertex'%N)

class DTcode:
    """
    Represents the DTcode of a link projection.
    Instantiate either from a list of tuples or
    an alphabetical code.
    """
    def __init__(self, code):
        if isinstance(code,str):
            code = self.convert_alpha(code)
        overcrossings = [sign(x) for comp in code for x in comp]
        evens = evens = [abs(x) for comp in code for x in comp]
        self.size = size = len(evens)
        odds =  range(1, 2*size, 2)
        self.pairs = pairs = zip(odds, evens)
        # Build a lookup table for vertices.
        # (DT codes are 1-based; we just waste the 0 entry.) 
        self.lookup = lookup = [None for n in range(1+2*size)]
        for pair, overcrossing in zip(pairs, overcrossings):
            m, n = pair
            V = DTvertex(pair, overcrossing)
            lookup[m] = lookup[n] = V
        # Build the fatgraph determined by the DT code.
        self.fat_graph = FatGraph()
        N = start = 1
        last_odd = -1
        V = self[1]
        for component in code:
            last_odd += 2*len(component)
            # Walk around a component, adding edges.
            while N <= last_odd:
                W = self[N + 1]
                self.fat_graph.add_edge((V, V.exit(N)),
                                        (W, W.enter(N+1)))
                N += 1
                V = W
            # Close it up and go to the next one.
            S = self[start]
            self.fat_graph.add_edge((V, V.exit(N)),
                                    (S, S.enter(start)))
            start = N
            V = W

    def __getitem__(self, n):
        """
        We can lookup a vertex by either label
        """
        return self.lookup[n]

    def flip(self, vertex):
        self.fat_graph.reorder(vertex, (0,3,2,1))
