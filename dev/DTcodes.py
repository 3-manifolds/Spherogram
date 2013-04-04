from spherogram import FatGraph

def sign(x):
    return 1 if x > 0 else -1 if x < 0 else 0

class DTvertex:
    """
    A vertex of the 4-valent graph described by a DT code.
    Instantiated with an even-odd pair, in either order.
    """
    def __init__(self, pair):
        x, y = pair
        self.small = x if x < y else y
        self.big = y if x < y else x
        self.even = x if y%2 else y
        self.odd = y if y%2 else x
    def __repr__(self):
        return str((self.small, self.big))
        
class DTcode:
    """
    Represents the DTcode of a link projection.
    Instantiate either from a list of tuples or
    an alphabetical code.
    """
    def __init__(self, code):
        if isinstance(code,str):
            code = self.convert_alpha(code)
        # a negative sign means that the even strand goes over
        self.overs = [sign(x) for comp in code for x in comp]
        self.evens = evens = [abs(x) for comp in code for x in comp]
        self.size = size = len(evens)
        self.odds =  range(1, 2*size, 2)
        self.pairs = pairs = zip(self.odds, self.evens)
        self.lookup = lookup = [None for n in range(1+2*size)]
        # just waste the 0 entry, since DT codes are 1-based
        for pair in pairs:
            m, n = pair
            V = DTvertex(pair)
            lookup[m] = lookup[n] = V
        # the length is the number of odds on a component
        self.lengths = lengths = [len(comp) for comp in code]
        edges = []
        
    
    def __getitem__(self, n):
        """
        We can lookup a vertex by either label
        """
        return self.lookup[n]
