"""
Computing Seifert matrices and putting Links in braid closure form. 

Code written by Malik Obeidin.  
"""

from .ordered_set import OrderedSet
from .simplify import reverse_type_II
from .links_base import Link  # Used for testing only
from .. import ClosedBraid    # Used for testing only
from itertools import combinations

def cyclic_permute(l,n):
    return [l[(i+n)%len(l)] for i in range(len(l))]

def seifert_circles(link):
    """
    Returns the circles in the diagram created by Seifert's algorithm

    >>> L = Link('L12n111')
    >>> [len(circle) for circle in seifert_circles(L)]
    [12, 4, 5, 3]
    """
    ceps = OrderedSet(link.crossing_entries())
    circles = []
    while ceps:
        start_cep = ceps.pop()
        circle = [start_cep]
        cep = start_cep.other().next()
        while cep != start_cep:
            circle.append(cep)
            ceps.remove(cep)
            cep = cep.other().next()
        circles.append(circle)
    return circles

def seifert_crossing_entry(crossing_strand):
    d = (crossing_strand.strand_index,(crossing_strand.strand_index+2)%4)
    if d in crossing_strand.crossing.directions:
        return crossing_strand
    else:
        entries = crossing_strand.crossing.entry_points()
        entries.remove(crossing_strand.rotate(2))
        return entries[0]

def admissible_moves(link):
    circles = seifert_circles(link)
    cs_to_seifert_circle = {}
    pairs = []
    seifert_circle_pairs = []
    for cs in link.crossing_strands():
        cep = seifert_crossing_entry(cs)
        for c in circles:
            if cep in c:
                break
        cs_to_seifert_circle[cs] = circles.index(c)
    for face in link.faces():
        for cs1, cs2 in combinations(face,2):
            circle1, circle2 = cs_to_seifert_circle[cs1], cs_to_seifert_circle[cs2]
            if circle1 != circle2:
                pairs.append((cs1,cs2))
                seifert_circle_pairs.append((circle1,circle2))
    return pairs, seifert_circle_pairs

def seifert_tree(link):
    """
    The oriented tree corresponding to the complementary regions of
    the Seifert circles.

    >>> T = seifert_tree(Link('K5a2'))
    >>> T == [(frozenset([0]), frozenset([0, 1])), (frozenset([0, 1]), frozenset([1]))]
    True
    """
    circles = seifert_circles(link)
    edges = [[set([n]),set([n])] for n in range(len(circles))]
    for c in link.crossings:
        under, over = c.entry_points()
        under_circle, over_circle = -1,-1
        sign = c.sign
        for n, circle in enumerate(circles):
            if under in circle:
                under_circle = n
            if over in circle:
                over_circle = n
            if under_circle>0 and over_circle>0:
                break
        if sign == -1:
            connect_head_to_tail(edges[under_circle],edges[over_circle])
        else:
            connect_head_to_tail(edges[over_circle],edges[under_circle])
    for e1, e2 in combinations(edges,2): #connect all vertices which intersect
        for i in range(2):
            for j in range(2):
                if len(e1[i].intersection(e2[j])) > 1:
                    connect_vertices(e1,i,e2,j)
    
    return [(frozenset(e[0]),frozenset(e[1])) for e in edges]

def remove_admissible_move(link):
    """
    Performs a Reidemester II move to remove one branching point of the Seifert
    tree.  The goal is to turn the Seifert tree into a chain.
    """
    circles = seifert_circles(link)
    moves, circle_pairs = admissible_moves(link)
    tree = seifert_tree(link)
    tails = [e[0] for e in tree]
    heads = [e[1] for e in tree]
    found_move = False
    for e1, e2 in combinations(tree,2):
        if e1[0] == e2[0]: #edges start at same point
            circles = set([tree.index(e1), tree.index(e2)])
            found_move = True
        elif e1[1] == e2[1]: #edges end at same point
            circles = set([tree.index(e1), tree.index(e2)])
            found_move = True
        if found_move:
            move_possible = False
            for n, pair in enumerate(circle_pairs):
                if set(pair) == circles:
                    cs1, cs2 = moves[n]
                    move_possible = True
                    break
            if move_possible:
                label1 = 'n'+str(cs1.crossing.label)
                label2 = 'n'+str(cs2.crossing.label)
                reverse_type_II(link,cs1,cs2,label1,label2)
                link._rebuild(same_components_and_orientations=True)
                break
            else:
                found_move = False
    return found_move

def isotope_to_braid(link):
    """
    Performs Reidemester II moves until the Seifert tree becomes a chain, i.e.
    the Seifert circles are all nested and compatibly oriented, following
    P. Vogel, "Representation of links by braids, a new algorithm"
    """
    while remove_admissible_move(link):
        pass

def is_chain(tree):
    tails = [e[0] for e in tree]
    heads = [e[1] for e in tree]
    if len(set(tails)) == len(tails) and len(set(heads)) == len(heads):
        return True
    else:
        return False

def connect_head_to_tail(e1,e2):
    e1[1] = e1[1]|e2[0]
    e2[0] = e1[1]|e2[0]

def connect_vertices(e1, v1, e2, v2):
    e1[v1] = e1[v1] | e2[v2]
    e2[v2] = e1[v1] | e2[v2]

def straighten_arrows(arrows):
    totally_straightened = False
    while not totally_straightened:
        totally_straightened = True
        for arrow in arrows:
            tail, head = arrow[0], arrow[1]
            if tail < head: #need to move tail down
                diff = head-tail
                same_start_strand = [x for x in arrows if x[2] == arrow[2] and x[0]>=tail]
                for other_arrow in same_start_strand:
                    other_arrow[0] += diff
                one_strand_behind = [x for x in arrows if x[2] == arrow[2]-1 and x[1]>=tail]
                for other_arrow in one_strand_behind:
                    other_arrow[1] += diff
                totally_straightened = False
            elif head < tail: #need to move head down
                diff = tail-head
                same_end_strand = [x for x in arrows if x[2] == arrow[2] and x[1]>=head]
                for other_arrow in same_end_strand:
                    other_arrow[1] += diff
                one_strand_ahead = [x for x in arrows if x[2] == arrow[2]+1 and x[0]>=head]
                for other_arrow in one_strand_ahead:
                    other_arrow[0] += diff
                totally_straightened = False

def braid_arrows(link):
    """
    Helper function to determine positions of all the crossings in a braid
    description of the link.
    """
    link_copy = link.copy()
    isotope_to_braid(link_copy)
    circles = seifert_circles(link_copy)
    tree = seifert_tree(link_copy)
    tails = [e[0] for e in tree]
    heads = [e[1] for e in tree]
    for t in tails:
        if t not in heads:
            start = tails.index(t)
            break

    ordered_strands = [circles[start]]
    for i in range(len(circles)-1):
        new_tail = tree[start][1]
        start = tails.index(new_tail)
        ordered_strands.append(circles[start])
    positions_in_next_strand = []

    for i in range(len(ordered_strands)-1):
        for n, cep in enumerate(ordered_strands[i]):
            found_next = False
            for m, next_cep in enumerate(ordered_strands[i+1]):
                if cep.crossing == next_cep.crossing:
                    ordered_strands[i+1] = cyclic_permute(ordered_strands[i+1],m)
                    found_next = True
                    break
            if found_next:
                break

    for i in range(len(ordered_strands)-1):
        positions = {}
        for n, cep in enumerate(ordered_strands[i]):
            for m, next_cep in enumerate(ordered_strands[i+1]):
                if cep.crossing == next_cep.crossing:
                    positions[n] = ((m,cep.strand_index%2))
                    break
        positions_in_next_strand.append(positions)

    ordered_strands = ordered_strands[::-1]
    arrows = [[i,positions[i][0],n,positions[i][1]] for n, positions in enumerate(positions_in_next_strand) for i in positions]
    straighten_arrows(arrows)
    arrows = sorted(arrows,key=lambda x: x[0])
    for arrow in arrows:
        arrow.pop(1) #start and end positions are now the same
    return arrows

def braid_word(link):
    """
    Return a list of integers which defines a braid word whose closure is the 
    given link.  The natural numbers 1, 2, 3, etc are the generators and the
    negatives are the inverses.

    >>> w = braid_word(Link('13n1234'))
    >>> M = ClosedBraid(w).exterior()
    >>> M.identify()
    [K13n1234(0,0)]

    Implementation follows P. Vogel, "Representation of links by
    braids, a new algorithm".
    """
    arrows = braid_arrows(link)
    word = []
    for position, strand, over_or_under in arrows:
        if over_or_under != 0:
            word.append(-strand-1)
        else:
            word.append(strand+1)
    return word

def seifert_matrix(link, return_matrix_of_types=False):
    """
    Returns the Seifert matrix of a link by first making it isotopic to a braid
    closure.

    >>> L = Link('K8n1')
    >>> seifert_matrix(L)  # doctest: +NORMALIZE_WHITESPACE
    [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
     [0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [-1, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
     [0, -1, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0],
     [0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, -1, 1, 0, 0, 1, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0]]

    Uses the algorithm described in:

    J. Collins, "An algorithm for computing the Seifert matrix of a link 
    from a braid representation." (2007).
    """
    arrows = braid_arrows(link)
    strands = set([x[1] for x in arrows])
    grouped_by_strand = [[x for x in arrows if x[1] == strand] for strand in strands]
    hom_gens = [[(group[i][0],group[i+1][0],group[i][2],group[i+1][2])
                 for i in range(len(group)-1)] for group in grouped_by_strand]
    num_gens = sum(map(len,hom_gens))
    matrix = [[0]*num_gens for i in range(num_gens)]
    entries = [(i,j) for i in range(len(hom_gens)) for j in range(len(hom_gens[i]))]
    type_matrix = [[0]*num_gens for i in range(num_gens)]
    for n, strand in enumerate(hom_gens):
        #diagonal entries
        for m, gen in enumerate(strand):
            if gen[2] == gen[3]: #same sign, otherwise entry is zero
                if gen[2] == 0: #both right handed
                    matrix[entries.index((n,m))][entries.index((n,m))] = -1
                    type_matrix[entries.index((n,m))][entries.index((n,m))] = 1
                else: #both left handed
                    matrix[entries.index((n,m))][entries.index((n,m))] = 1
                    type_matrix[entries.index((n,m))][entries.index((n,m))] = 2

        
        #two gens on same strand, one after the other
        for m, gen in enumerate(strand[:-1]):
            if gen[3] == 0: #shared crossing is right handed
                matrix[entries.index((n,m+1))][entries.index((n,m))] = 1
                type_matrix[entries.index((n,m+1))][entries.index((n,m))] = 3
            else: #shared crossing if left handed
                matrix[entries.index((n,m))][entries.index((n,m+1))] = -1
                type_matrix[entries.index((n,m))][entries.index((n,m+1))] = 4

        #two gens on adjacent strand, "staggered"
        if n != len(hom_gens)-1:
            next_strand = hom_gens[n+1]
            for m, gen in enumerate(strand):
                for l, next_gen in enumerate(next_strand):
                    if next_gen[0] < gen[0] < next_gen[1] < gen[1]:
                        matrix[entries.index((n+1,l))][entries.index((n,m))] = 1
                        type_matrix[entries.index((n+1,l))][entries.index((n,m))] = 5
                    elif gen[0] < next_gen[0] < gen[1] < next_gen[1]:
                        matrix[entries.index((n+1,l))][entries.index((n,m))] = -1
                        type_matrix[entries.index((n+1,l))][entries.index((n,m))] = 6
    
    if return_matrix_of_types:
        return matrix, type_matrix
    else:
        return matrix



