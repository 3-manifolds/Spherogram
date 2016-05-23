import spherogram
from spherogram.links.tangles import Tangle, OneTangle, MinusOneTangle
import networkx as nx
from random import randint,choice,sample
from spherogram.links.random_links import map_to_link, random_map

"""
This file contains some unfinished code to work with tangles, including 
getting tangles out of link diagrams, as well as some Conway mutation.
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


def add_random_crossing(self,label):
    """
    Randomly chooses position on boundary of the tangle and splits into a new
    crossing.
    """

    tangle_copy = self.copy()
    adj = tangle_copy.adjacent
    adj[len(adj)/2:] = reversed(adj[len(adj)/2:])
    new_crossing = spherogram.Crossing(label)
    old_position = randint(0,len(adj)-1)
    old_crossing, old_strand = adj.pop(old_position)
    new_strand = randint(0,3)
    old_crossing[old_strand]=new_crossing[new_strand]
    for i in range(1,4):
        adj.insert(old_position,(new_crossing,(new_strand-i)%4))
    adj[len(adj)/2:] = reversed(adj[len(adj)/2:])    
    tangle_copy.crossings.append(new_crossing)
    tangle_copy.n = self.n+1
    return tangle_copy

def random_tree(size):
    """
    Repeatedly splits crossings, starting with a single crossing, resulting
    in a random 4-valent tree of crossings.
    """

    T = OneTangle()
    T.crossings[0].label = 0
    for i in range(1,size):
        T = T.add_random_crossing(i)
    return T

def random_tree_link(size):
    """
    Generate two random tree tangles and glues together.  Gives a link of size
    2*size
    """

    return random_tree(size).circular_sum(random_tree(size),0)


def random_tree_knot(size,simplify=None,prime_decomp = False):
    found_nontrivial = False
    while not found_nontrivial:
        link = random_tree_link(size)
        knot = link.sublink(max(link.link_components,key=len))
        knot.simplify(mode=simplify)
        if len(knot) > 0:
            found_nontrivial=True
    
    if prime_decomp:
        cant_deconnect = False
        while cant_deconnect:
            ds = knot.deconnect_sum()
            knot = max(ds,key=len)
            cant_deconnect = (len(ds)>1)
    return knot

def all_circular_sums(self,other):
    """
    All possible tangle sums as above
    """

    if len(self.adjacent) != len(other.adjacent):
        raise Exception("Tangles do not have the same number of strands")
    return [self.circular_sum(other,n) for n in range(len(other.adjacent))]



def random_root(self):
    return choice(choice(self.crossings).crossing_strands())


def isosig_with_gluings(self, gluings, root=None):
    return (self.isosig(root = root),tuple(gluings))

def min_isosig_with_gluings(self, gluings, root = None):
    if root != None:
        cs_name = cslabel(root)
    isosigs = []
    for i in range(self.n*2):
        rotated_tangle = self.circular_rotate(i)
        if root != None:
            rotated_root = crossing_strand_from_name(rotated_tangle,cs_name)
        else:
            rotated_root = None
        #permuting the indices in the gluings
        perm = range(len(self.adjacent))
        perm[len(perm)/2:] = reversed(perm[len(perm)/2:])
        perm = rotate_list(perm,i)
        perm[len(perm)/2:] = reversed(perm[len(perm)/2:])
        rotated_gluings = []
        for g in gluings:
            new_g = [perm[g[0]],perm[g[1]]]
            new_g.sort()
            rotated_gluings.append(tuple(new_g))
        rotated_gluings.sort()
        isosigs.append(rotated_tangle.isosig_with_gluings(rotated_gluings,root=rotated_root))        


    return min(isosigs)

Tangle.all_circular_sums = all_circular_sums
Tangle.add_random_crossing = add_random_crossing

def cycle_basis(G):
    """
    Uses networkx's cycle basis function on dual graph and converts to 
    form with spherogram objects
    """

    Gx = G.to_networkx()
    vert_cycles = nx.cycle_basis(Gx)
    return [edge_cycle(vert_cycle,G) for vert_cycle in vert_cycles]

def is_trivial(four_cycle):
    """
    A trivial four cycle in the dual graph bounds a single quadrilateral
    """
    crossings = map(lambda x: map(lambda y: y.crossing,x.interface), four_cycle)
    return len(set(crossings[0])&set(crossings[1])&set(crossings[2])&set(crossings[3])) != 0


def all_four_cycles_at_vertex(G, start_vertex):
    adjacent = G.children(start_vertex)
    four_cycles = set()
    for v in adjacent:
        for w in adjacent:
            if v==w:
                continue
            new_adj = G.children(v).intersection(G.children(w))
            new_adj.remove(start_vertex)
            for opposite_vertex in new_adj:
                for e1 in G.edges_between(start_vertex,v):
                    for e2 in G.edges_between(v, opposite_vertex):
                        for e3 in G.edges_between(opposite_vertex, w):
                            for e4 in  G.edges_between(w,start_vertex):
                                four_cycle = (e1,e2,e3,e4)
                                if is_trivial(four_cycle):
                                    continue
                                else:
                                    four_cycles.add(four_cycle)

    return four_cycles

def unknot_search(num_attempts, backtrack_height, num_mutations):
    c = spherogram.Crossing(0)
    c[0]=c[1]
    c[2]=c[3]
    U = spherogram.Link([c])
    for i in range(num_attempts):
        print(i)
        Uc = U.copy()
        Uc.backtrack(backtrack_height)
        for i in range(num_mutations):
            Uc = random_mutate(Uc)
        Uc.simplify(mode='level')
        if len(Uc)>0:
            return Uc
    return None


def get_four_cycle(G, start_vertex):
    """
    Returns the first nontrivial 4-cycle found in the graph G
    """
    adjacent = G.children(start_vertex)
    for v in adjacent:
        for w in adjacent:
            if v==w:
                continue
            new_adj = G.children(v).intersection(G.children(w))
            new_adj.remove(start_vertex)
            for opposite_vertex in new_adj:
                for e1 in G.edges_between(start_vertex,v):
                    for e2 in G.edges_between(v, opposite_vertex):
                        for e3 in G.edges_between(opposite_vertex, w):
                            for e4 in  G.edges_between(w,start_vertex):
                                four_cycle = (e1,e2,e3,e4)
                                if is_trivial(four_cycle):
                                    continue
                                else:
                                    return four_cycle
    return []

def all_four_cycles(G):
    """
    Returns all possible four cycles in G
    """
    four_cycles = [x for v in G.vertices for x in all_four_cycles_at_vertex(G,v)]
    four_cycles_no_duplicates = []
    for fc in four_cycles:
        seen_before = False
        for seen_fc in four_cycles_no_duplicates:
            if set(fc) == set(seen_fc):
                seen_before = True
                break
        if not seen_before:
            four_cycles_no_duplicates.append(fc)
    return four_cycles_no_duplicates

def edge_cycle(vert_list,G):
    """
    Converts from list of vertices of dual graph to list of edges.  
    If multiple edges, just chooses one.
    """
    edges = list(G.edges)
    cycle = []
    for i in range(len(vert_list)-1):
        face_pair = [vert_list[i],vert_list[i+1]]
        for edge in edges:
            if set(face_pair) == set(edge.incident_to()):
                cycle.append(edge)
                break
    face_pair = [vert_list[0],vert_list[-1]]
    for edge in edges:
        if set(face_pair) == set(edge.incident_to()):
            cycle.append(edge)
            break
    return cycle

def crossing_ball(crossing,radius):
    """
    Returns the crossings within distance r of the crossing, in the form
    of a dictionary, where the values are the distances to the center crossing,
    and a list of the crossing strands along the boundary.
    """
    distances = {crossing: 0}
    opposite_positions = [cs.opposite() for cs in crossing.crossing_strands() if cs.opposite().crossing != crossing]
    for i in range(1,radius):
        new_opposites = []
        for cs in opposite_positions:
            if cs.crossing not in distances:
                distances[cs.crossing] = i
            new_opposites.append(cs.rotate(1).opposite())
            new_opposites.append(cs.rotate(2).opposite())
            new_opposites.append(cs.rotate(3).opposite())
        op_repeats = [x for x in new_opposites if x.crossing not in distances]
        opposite_positions = []
        for cs in op_repeats:
            if cs not in opposite_positions:
                opposite_positions.append(cs)
    return distances, map(lambda x: x.opposite(), opposite_positions)


def boundary_components(link,crossing,radius):
    crossings, adjacent = crossing_ball(crossing,radius)
    crossings = list(crossings)
    G = underlying_graph(link)
    for c in crossings:
        G.remove_node(c)
    return nx.connected_components(G)

def underlying_graph(link):
    G = nx.Graph()
    edges = [(c,adj) for c in link.crossings for adj in map(lambda x: x[0],c.adjacent)]
    G.add_edges_from(edges)
    return G

def trace_boundary_component(start_cs,full_boundary):
    boundary_comp = [start_cs]
    cs = start_cs.next_corner()
    i = 0
    while cs != start_cs:
        while cs.rotate(1) not in full_boundary:
            print(cs)
            cs = cs.next_corner()
            i += 1
            if i>100: raise Exception()
        cs = cs.rotate(1)
        boundary_comp.append(cs)
    boundary_comp.pop(-1) #code aboves adds the start_cs twice
    return boundary_comp

    

def tangle_neighborhood(link,crossing,radius,return_gluings=True,hull=False):
    """
    Splits a link into two tangles along a ball around a crossing of the given
    radius.  Destroys original link.  This might not generate an actual tangle;
    the graph metric ball will usually have multiple boundary components.
    """

    crossings, adjacent = crossing_ball(crossing,radius)
    crossings = list(crossings)

    opposites = list(reversed(map(lambda x: x.opposite(),adjacent)))
    outside_crossings = [c for c in link.crossings if c not in crossings]
    if len(outside_crossings) == 0:
        raise Exception("Neighborhood is entire link")
    n = len(adjacent)/2

    if hull:
        comps = list(boundary_components(link,crossing,radius))
        largest_comp = max(comps)
        sides = dict([(cslabel(cross_strand), cross_strand) for cross_strand in adjacent])
        c = largest_comp.pop()
        cs = choice(c.crossing_strands())
        exit_strand = meander(cs,sides)[1]
        exit_strand = exit_strand[0].crossing_strands()[exit_strand[1]]
        main_boundary_comp = trace_boundary_component(exit_strand,adjacent)
        print('main_boundary_comp' + str(main_boundary_comp))
        print('all comps: '+str(comps))
        comps.remove(largest_comp) #remove largest component
        for comp in comps:
            print('crossings: ' +str(crossings))
            print('filling in comp:'+str(comp))
            print('adjacent: '+str(adjacent))
            c = comp.pop()
            cs = choice(c.crossing_strands())
            
            print('cs: ' + str(cs))
            exit_strand = meander(cs,sides)[1] #meander until you hit boundary
            exit_strand = exit_strand[0].crossing_strands()[exit_strand[1]]
            print('exit_strand: ' +str(exit_strand))
            bound_comp = trace_boundary_component(exit_strand,adjacent)
            print('traced component: ' + str(bound_comp))
            if exit_strand not in main_boundary_comp:
                for x in bound_comp:
                    adjacent.remove(x)
                print('updated adjacent: ' +str(adjacent))
                crossings.append(c)
                crossings.extend(list(comp))
            
    adjacent[n:] = reversed(adjacent[n:])
    opposites[n:] = reversed(opposites[n:])
    gluings = []
    seen_cs = []
    #figure out which strands are glued to each other
    for cs in adjacent:
        if cs in seen_cs: continue
        next_cross = cs.next_corner()
        while next_cross not in adjacent:
            next_cross = next_cross.next_corner()
        if next_cross != cs:
            gluings.append((adjacent.index(cs),adjacent.index(next_cross)))
            seen_cs.append(next_cross)
    clear_orientations(crossings)
    clear_orientations(outside_crossings)
    gluings.sort()
    if return_gluings:
        return Tangle(n,crossings,adjacent),Tangle(n,outside_crossings,opposites),gluings
    else:
        return Tangle(n,crossings,adjacent),Tangle(n,outside_crossings,opposites)

def mutate(link,four_cycle):
    link_copy = link.copy()
    four_cycle_copy = [corresponding_edge(link_copy,edge) for edge in four_cycle]
    T1,T2 = tangle_cut(link_copy,four_cycle_copy)
    #print(len(T1.crossings),len(T2.crossings))
    return T1.circular_sum(T2,2)

def all_mutants(link):
    G = link.dual_graph()
    return [mutate(link,four_cycle) for four_cycle in all_four_cycles(G)]

def random_mutate(link):
    G = link.dual_graph()
    #afc = all_four_cycles(G)
    #if len(afc)==0:
    #    return link.copy()
    #four_cycle = choice(afc)
    v = choice(tuple(G.vertices))
    all_fcs = tuple(all_four_cycles_at_vertex(G,v))
    if len(all_fcs) == 0:
        return link.copy()
    four_cycle = choice(all_fcs)
    return mutate(link,four_cycle)


def corresponding_edge(new_link,edge):
    G = new_link.dual_graph()
    for new_edge in G.edges:
        if str(new_edge) == str(edge):
            return new_edge

def tangle_cut(link, cycle):
    """
    Creates two Tangle objects from a given cycle (with no self intersections)
    in the dual graph of a link (inside and outside).  Cycle is given
    as a list of (oriented) edges in the dual graph. Make sure crossings 
    are uniquely labeled. Destroys original link.
    """

    sides = {} #keeps track of which side each crossing strand is as dictionary
    sides[cslabel(cycle[0].interface[0])] = 0 #start by assigning to first pair
    sides[cslabel(cycle[0].interface[1])] = 1 #sides
    side0 = [cycle[0].interface[0]]
    side1 = [cycle[0].interface[1]]
    for i in range(len(cycle)-1):
#        print(i)
        edge1 = cycle[i]
        edge2 = cycle[i+1]
        edge1_cs1, edge1_cs2 = edge1.interface #crossing strands on each end
        edge2_cs1, edge2_cs2 = edge2.interface
        edge1_cs1_side = sides[cslabel(edge1_cs1)]
        edge1_cs2_side = 1-edge1_cs1_side

        wrong_face = False
        #starting with edge1_cs1, travel around face
        #if I run into one of edge2_cs*, they are on different sides
        #if I don't, try again with the other
        while True:
#            print('sides:')
#            print(sides)
#            print('edge1:')
#            print(edge1_cs1,edge1_cs2)
#            print('edge2:')
#            print(edge2_cs1,edge2_cs2)
            edge1_cs1 = edge1_cs1.next_corner()
#            print(type(edge1_cs1.crossing.label))
#            print(edge1_cs1)
            if edge1_cs1 == edge2_cs1:
                sides[cslabel(edge2_cs1)] = edge1_cs2_side
                sides[cslabel(edge2_cs2)] = edge1_cs1_side
                if edge1_cs2_side == 1:
                    side1.append(edge2_cs1)
                    side0.append(edge2_cs2)
                else:
                    side0.append(edge2_cs1)
                    side1.append(edge2_cs2)                    
#                print('case 1')
                break
            if edge1_cs1 == edge2_cs2:
                sides[cslabel(edge2_cs1)] = edge1_cs1_side
                sides[cslabel(edge2_cs2)] = edge1_cs2_side
                if edge1_cs2_side == 1:
                    side0.append(edge2_cs1)
                    side1.append(edge2_cs2)
                else:
                    side1.append(edge2_cs1)
                    side0.append(edge2_cs2)                    
#                print('case 2')
                break
            if edge1_cs1.opposite() == edge1_cs2: #returned without seeing them
                wrong_face = True
#                print('case 3')
                break

        # Try again with other side
        while wrong_face:
            edge1_cs1_side = sides[cslabel(edge1_cs1)]
            edge1_cs2_side = 1-edge1_cs1_side
            edge1_cs2 = edge1_cs2.next_corner()
            if edge1_cs2 == edge2_cs1:
#                print('case 4')
                sides[cslabel(edge2_cs1)] = edge1_cs1_side
                sides[cslabel(edge2_cs2)] = edge1_cs2_side
                if edge1_cs2_side == 1:
                    side0.append(edge2_cs1)
                    side1.append(edge2_cs2)
                else:
                    side1.append(edge2_cs1)
                    side0.append(edge2_cs2)                    
                break
            if edge1_cs2 == edge2_cs2:
 #               print('case 5')
                sides[cslabel(edge2_cs1)] = edge1_cs2_side
                sides[cslabel(edge2_cs2)] = edge1_cs1_side
                if edge1_cs2_side == 1:
                    side1.append(edge2_cs1)
                    side0.append(edge2_cs2)
                else:
                    side0.append(edge2_cs1)
                    side1.append(edge2_cs2)
                break
            if edge1_cs2.opposite() == edge1_cs1: #returned without seeing them again, must be error
                raise Exception("Neither side worked")
    
    crossing_sides = fill_in_crossings(link,sides)
    n = len(cycle)
    side0[n/2:] = reversed(side0[n/2:]) #flip to use as adjacent in tangle
    side1[n/2:] = reversed(side1[n/2:])
    crossings0 = [crossing_from_name(link,c) for c in crossing_sides if crossing_sides[c] == 0]
    crossings1 = [crossing_from_name(link,c) for c in crossing_sides if crossing_sides[c] == 1]
    
    #clear crossing info
    clear_orientations(crossings0)
    clear_orientations(crossings1)

    #One of the tangles is the 'outside', and needs to be flipped
    #Just check side0
    side0_needs_flip = False
    c,i = side0[0]
    while True:
        next_cep = c.crossing_strands()[(i+1)%4]
        c,i = next_cep.crossing,next_cep.strand_index
#        print(c,i)
#        print(side0)
        if (c,i) in side0:
#            print('hit_end')
#            print(c,i)
#            print(side0[1])
#            print((c,i)==side0[1])
            side0_needs_flip = (c,i) != (side0[1])
            break
        c,i = next_cep.opposite().crossing,next_cep.opposite().strand_index
    if side0_needs_flip:
#        print('flipped side 0')
        flip(side0)
    else:
#        print('flipped side 1')
        flip(side1)
    return Tangle(n/2,crossings0,side0),Tangle(n/2,crossings1,side1)

def fill_in_crossings(link,sides):
    """
    Given boundary as above , fill in crossings on either side from sides.
    Returns a dictionary with the side (0 or 1) of each crossing.
    """

    crossing_sides = dict([(x[0], sides[x]) for x in sides])
    crossing_labels = map(lambda c: c.label,link.crossings)
    crossings_to_sort = set(crossing_labels)-set(x[0] for x in sides)
    while len(crossings_to_sort)>0:
        start_crossing = crossings_to_sort.pop()
        accumulated_crossings = [start_crossing]
        m,end_side = meander(crossing_strand_from_name(link,(start_crossing,randint(0,3))),sides)
        accumulated_crossings.extend(map(lambda x: x.label,m))
        for c in accumulated_crossings:
            crossing_sides[c]=end_side
    return crossing_sides


def meander(cs,sides):
    """
    Wander randomly starting at crossing strand cs until you hit
    a boundary strand in sides.Assumes the crossing of cs is not on the side.
    Returns a set of all the crossings encountered along the way,
    including ones on the boundary, and which side (0 or 1) is hit
    """
    crossings_encountered = [cs.crossing]
    end_side = 0
    while True:
        cs = cs.opposite().rotate(randint(1,3))
        if cslabel(cs) in sides: # hit the side
            end_side = sides[cslabel(cs)]
            break
        crossings_encountered.append(cs.crossing)
    return set(crossings_encountered),end_side
    
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



def crossing_by_label(label,link):
    for c in link.crossings:
        if c.label == label:
            return c

def all_neighborhoods(link,radius):
    nhds = []
    i = 0
    for c in link.crossings:
        print(i)
        i += 1
        link_copy = link.copy()
        c_copy = crossing_by_label(c.label,link_copy)
        T1, T2 = tangle_neighborhood(link_copy,c_copy,radius)
        nhds.append(T1)
    return nhds

def neighborhood_distribution(link,radius):
    nhds = all_neighborhoods(link,radius)
    nhd_classes = []
    for nhd in nhds:
        already_found = False
        for nhd_class in nhd_classes:
            if nhd.is_isotopic(nhd_class[0]):
                nhd_class[1] += 1
                already_found = True
                break
        if not already_found:
            nhd_classes.append([nhd,1])
    return nhd_classes
