"""
Switching crossings to make each twist region consistent.  
"""

from .links import CrossingStrand
from ..graphs import CyclicList

class  TwistRegionCap:
    def __init__(self, crossing_strand):
        self.cs = crossing_strand

    def sign(self):
        return self.cs.strand_index % 2

    def swap_crossing(self):
        self.cs.crossing.rotate_by_90()
    

class TwistRegionCrossing(TwistRegionCap):
    """
    A crossing together with a chosen bigon.  Recorded by the CrossingStrand
    which gives the side of said bigon which is "most clockwise" when viewed from
    the crossing.  
    """
    def __init__(self, crossing):
        if isinstance(crossing, CrossingStrand):
            self.cs = crossing
        else:
            assert is_end_of_twist_region(crossing)
            neighbors = CyclicList(C for (C, i) in crossing.adjacent)
            for i in range(4):
                if neighbors[i] == neighbors[i+1]:
                    self.cs = CrossingStrand(crossing, i)

    def next(self):
        # Move along edge to next crossing and rotate to preferred side of bigon.
        cs = self.cs.opposite().rotate(1)
        C, e = cs
        if C.adjacent[e][0] != C.adjacent[e+1][0]:
            return TwistRegionCap(cs)
        return TwistRegionCrossing(cs)

    def __repr__(self):
        return "<Twist: %s>" % (self.cs, )

class TwistRegion:
    def __init__(self, crossing):
        C = TwistRegionCrossing(crossing)
        crossings = [C]
        while isinstance(C, TwistRegionCrossing):
            C = C.next()
            crossings.append(C)
                        
        self.crossings = crossings

    def signs(self):
        return [C.sign() for C in self.crossings]
    
    def make_consistent(self):
        sign = self.crossings[0].sign()
        for C in self.crossings:
            if C.sign() != sign:
                C.swap_crossing()
        
    def __len__(self):
        return len(self.crossings)
    
def is_end_of_twist_region(crossing):
    return len(set([C for (C, i) in crossing.adjacent])) == 3

def make_twist_regions_consistent(link):
    """
    Changes crossings so that no bigon permits a Type II move cancelling
    the pair of crossings at the end of the bigon.  The system is that the 
    end crossing with the lowest label "wins" in determining the type of the twist
    region.  The code assumes that no link component is effectively a meridian
    loop for another component; said differently, no two bigons share a
    common edge.

    Note, this implementation fails if given something like a (2, n) torus link.
    """
    ans = []
    for C in link.crossings:
        if is_end_of_twist_region(C):
            twist = TwistRegion(C)
            twist.make_consistent()
    link._rebuild()
