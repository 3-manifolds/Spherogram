import spherogram

"""
This file contains functions to distinguish tangle diagrams up to planar 
isotopy. It can allow for roots, as well as distinguishing the shadows or 
the tangle diagrams with the information of whether a crossing is over or under.
"""

def rotate_list(L, s):
    n = len(L)
    return [ L[(i + s) % n] for i in range(n) ]

def flip(L):
    half = len(L)/2
    for i in range(half):
        L[i],L[i+half] = L[i+half],L[i]

def clear_orientations(crossings):
    for c in crossings:
        c.directions.clear()
        c.sign = 0


def isosig(tangle,root = None, over_or_under = False):
    strands, loops, orientations, crossing_order, over_or_under_data = all_cross_strands(tangle)
    isosig_strands = []
    isosig_loops = []
    if root != None:
        want_root = True
        root_position = crossing_order.index(root.crossing)
        root_tuple = tuple(root)
        root_tuple_op = tuple(root.rotate(2))
    else:
        want_root = False
    found_root = False
    root_isosig = None
    for strand,end in strands:
        isst = map(lambda cs: crossing_order.index(cs[0]),strand)
        if (want_root) and (not found_root) and (root_position in isst):
            root_strand_indices = [i for i in range(len(isst)) if isst[i]==root_position]
            for root_strand_index in root_strand_indices:
                if root_tuple == strand[root_strand_index]:
                    found_root = True
                    root_isosig = (root_position,(strands.index((strand,end)),root_strand_index),-1,'s')
                elif root_tuple_op == strand[root_strand_index]:
                    found_root = True
                    root_isosig = (root_position,(strands.index((strand,end)),root_strand_index),1,'s')
        isosig_strands.append((tuple(isst),end))
    for loop in loops:
        isl = map(lambda cs: crossing_order.index(cs[0]),loop)
        if (want_root) and (not found_root) and (root_position in isl):
            root_loop_indices = [i for i in range(len(isl)) if isl[i]==root_position]
            for root_loop_index in root_loop_indices:
                if root_tuple == loop[root_loop_index]:
                    found_root = True
                    root_isosig = (root_position,(loops.index(loop),root_loop_index),-1,'l')
                elif root_tuple_op == loop[root_loop_index]:
                    found_root = True
                    root_isosig = (root_position,(loops.index(loop),root_loop_index),1,'l')
        isosig_loops.append(tuple(isl))
    isosig_over_or_under = []
    if over_or_under:
        isosig_over_or_under = [over_or_under_data[c] for c in crossing_order]
    isosig_orientations = [orientations[c] for c in crossing_order]
    return (len(tangle.crossings),tangle.n),tuple(isosig_strands),tuple(isosig_loops),tuple(isosig_orientations), tuple(isosig_over_or_under), root_isosig

def min_isosig(tangle,root=None,over_or_under=False):
    if root != None:
        cs_name = cslabel(root)
    isosigs = []
    for i in range(tangle.n*2):
        rotated_tangle = tangle.circular_rotate(i)
        if root != None:
            rotated_root = crossing_strand_from_name(rotated_tangle,cs_name)
        else:
            rotated_root = None
        isosigs.append(isosig(rotated_tangle, root=rotated_root,over_or_under=over_or_under))
    return min(isosigs)


def isosig_with_gluings(tangle, gluings, root=None):
    return (isosig(tangle,root = root),tuple(gluings))

def min_isosig_with_gluings(tangle, gluings, root = None):
    if root != None:
        cs_name = cslabel(root)
    isosigs = []
    for i in range(tangle.n*2):
        rotated_tangle = tangle.circular_rotate(i)
        if root != None:
            rotated_root = crossing_strand_from_name(rotated_tangle,cs_name)
        else:
            rotated_root = None
        #permuting the indices in the gluings
        perm = range(len(tangle.adjacent))
        perm[len(perm)/2:] = reversed(perm[len(perm)/2:])
        perm = rotate_list(perm,i)
        perm[len(perm)/2:] = reversed(perm[len(perm)/2:])
        rotated_gluings = []
        for g in gluings:
            new_g = [perm[g[0]],perm[g[1]]]
            new_g.sort()
            rotated_gluings.append(tuple(new_g))
        rotated_gluings.sort()
        isosigs.append(isosig_with_gluings(rotated_tangle, rotated_gluings,root=rotated_root))        


    return min(isosigs)


def crossing_orientations(strands):
    """
    Given the strands, compute the orientations (+1 or -1) for each crossing
    in the tangle
    """
    orientations = {}
    over_or_under = {}
    css_seen = []
    for strand in strands:
        for cs in strand:
            for seen_cs in css_seen: 
                if cs[0] == seen_cs[0]:
                    orientation = (cs[1]-seen_cs[1])%4
                    if orientation == 3:
                        orientation = -1
                    orientations[cs[0]]=orientation
                    over_or_under[cs[0]] = (cs[1]%2)
                    break
            css_seen.append(cs) #didn't find cs
    return orientations, over_or_under

"""
Give a list of the crossing strands encountered starting at 
a strand on the boundary of a tangle and moving to the other
end of that strand.
"""
def cross_strand(tangle,i):
    if i >= 2*tangle.n:
        raise Exception("Not a valid start position for strand")
    cs = tangle.adjacent[i]
    strand = [cs]
    while (cs[0],(cs[1]+2)%4) not in tangle.adjacent:
        cs = cs[0].adjacent[(cs[1]+2)%4]
        strand.append(cs)
    return strand

def loop_strand(cs):
    """
    Get the closed loop starting at crossing strand cs
    """
    strand = [cs]
    cs = cs[0].adjacent[(cs[1]+2)%4]
    while cs not in strand:
        strand.append(cs)
        cs = cs[0].adjacent[(cs[1]+2)%4]
    return strand

def all_cross_strands(tangle):
    """
    Returns all the strands but without duplicate in the opposite direction,
    starting at position 0 and going clockwise, and then components that 
    don't intersect the boundary.
    """
    other_ends_seen = [] #use to eliminate duplicate strand in other direction
    strands = []
    strands_with_ends = []
    loops = []
    clockwise_order = range(tangle.n)
    clockwise_order.extend(reversed(range(tangle.n,tangle.n*2)))
    for i in clockwise_order:
        if i not in other_ends_seen:
            strand = cross_strand(tangle,i)
            cs = strand[-1]
            end = tangle.adjacent.index((cs[0],(cs[1]+2)%4))
            if end not in other_ends_seen:
                strands.append(strand)
                strands_with_ends.append((strand,end))
                other_ends_seen.append(end)
    orientations, over_or_under = crossing_orientations(strands)
    cs_seen = [cs for strand in strands for cs in strand]
    seen_once = set(cs[0] for cs in cs_seen)
    for crossing in orientations:
        seen_once.remove(crossing)
    for strand in strands:
        for cs in strand:
            if cs[0] in seen_once:
                loop = loop_strand((cs[0],(cs[1]+1)%4))
                loops.append(loop)
                cs_seen.extend(loop)
                for loop_cs in loop:
                    if loop_cs[0] in seen_once:
                        for seen_cs in cs_seen:
                            if loop_cs[0] == seen_cs[0]:
                                orientation = (loop_cs[1]-seen_cs[1])%4
                                if orientation == 3:
                                    orientation = -1
                                orientations[loop_cs[0]]=orientation
                                over_or_under[loop_cs[0]] = loop_cs[1]%2
                                seen_once.remove(loop_cs[0])
                                break
                    else:
                        seen_once.add(loop_cs[0])
    while len(orientations)<len(tangle.crossings):
        for loop in loops:
            for cs in loop:
                if cs[0] in seen_once:
                     loop = loop_strand((cs[0],(cs[1]+1)%4))
                     loops.append(loop)
                     cs_seen.extend(loop)
                     for loop_cs in loop:
                         if loop_cs[0] in seen_once:
                             for seen_cs in cs_seen:
                                 if loop_cs[0] == seen_cs[0]:
                                     orientation = (loop_cs[1]-seen_cs[1])%4
                                     if orientation == 3:
                                         orientation = -1
                                     orientations[loop_cs[0]]=orientation
                                     over_or_under[loop_cs[0]]= loop_cs[1]%2
                                     seen_once.remove(loop_cs[0])
                                     break
                         else:
                             seen_once.add(loop_cs[0])
    crossings = map(lambda x: x[0], cs_seen)
    crossing_order = []
    for c in crossings:
        if c not in crossing_order:
            crossing_order.append(c)
    return strands_with_ends,loops,orientations,crossing_order, over_or_under

def cslabel(cs):
    """
    Label of crossing strand, without frills
    """
    return (cs[0].label,cs[1])

def crossing_strand_from_name(link,csname):
    """
    Find crossing strand object from it's name in the format of cslabel above
    """
    for c in link.crossings:
        if c.label == csname[0]:
            return c.crossing_strands()[csname[1]]
    raise Exception("Crossing not found")

def crossing_from_name(link,crossname):
    for c in link.crossings:
        if c.label == crossname:
            return c
    raise Exception("Crossing not found")

def analogous_crossing_strand(link,cs):
    """
    Get a crossing in link with the same label as cs
    """
    return crossing_strand_from_name(link,cslabel(cs))

def cut_strand(link, cs):
    """
    Cut a link along a strand to get a tangle with 1 strand
    """
    link_copy = link.copy()
    cs_copy = crossing_strand_from_name(link_copy,cslabel(cs))
    op = cs_copy.opposite()
    cs_copy.crossing.adjacent[cs_copy.strand_index] = None
    op.crossing.adjacent[op.strand_index] = None
    T = spherogram.Tangle(1, link_copy.crossings, [(cs_copy.crossing,cs_copy.strand_index),(op.crossing,op.strand_index)])
    return T


def link_isosig(link, root = None, over_or_under = False ):
    """
    A list of data which encodes a planar isotopy class of links instead of
    tangles.  This is just the minimal isosig gotten by cutting all possible
    strands to get a tangle with 1 strand.
    """
    isosigs = [cut_strand(link,cs).isosig(root,over_or_under) for cs in link.crossing_strands()]
    return min(isosigs)
