from snappy import *
from spherogram import *
from sage.all import *

def braid_closure_via_tangles(braid):
    n = max(-min(braid),max(braid))+1
    l = len(braid)
    B = IdentityBraid(n)
    for i in range(l):
        v = braid[i]
        nextpiece = IdentityBraid(abs(v)-1)|RationalTangle(sign(v))|IdentityBraid(n-abs(v)-1)
        B = B*nextpiece
    return B.braid_closure()

def braid_closure_via_crossings(braid):
    """
    Compute the braid closure of a word given in the form of a list of integers,
    where 1, 2, 3, etc correspond to the generators sigma_1, sigma_2, sigma_3, 
    and so on, and negative numbers to their inverses.
    """
    l = len(braid)
    crossings = [Crossing(i) for i in range(l)]
    for i in range(l):
        foundleftstrand = False
        foundrightstrand = False
        j = (i+1)%l

        while ((not foundleftstrand) or (not foundrightstrand)):
            if (not foundleftstrand) and (abs(braid[j]) == abs(braid[i])-1):
                crossings[i][3] = crossings[j][1]
                foundleftstrand = True
            if (not foundrightstrand) and (abs(braid[j]) == abs(braid[i])+1):
                crossings[i][0] = crossings[j][2]
                foundrightstrand = True
            if abs(braid[j]) == abs(braid[i]):
                if (not foundleftstrand):
                    crossings[i][3] = crossings[j][2]
                    foundleftstrand = True
                if (not foundrightstrand):
                    crossings[i][0] = crossings[j][1]
                    foundrightstrand = True
            j = (j+1)%l
            
    for i in range(l):
        if braid[i]>0:
            crossings[i].rotate_by_90()
    return Link(crossings)
