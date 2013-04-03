def sign(x):
    return 1 if x > 0 else -1 if x < 0 else 0

class DTcode:
    """
    Represents the DTcode of a link projection.
    Instantiate either from a list of tuples or
    an alphabetical code.
    """
    def __init__(self, code):
        if isinstance(code,str):
            code = self.convert_alpha(code)
        self.overs = [sign(x) for comp in code for x in comp]
        self.lengths = [len(comp) for comp in code]
        self.evens = [abs(x) for comp in code for x in comp]
        self.odds =  range(1, 2*len(self.evens), 2)
        self.pairs = zip(self.odds, self.evens)
