import spherogram
from spherogram.links.tangles import Tangle, OneTangle, MinusOneTangle
import networkx as nx
from random import randint,choice,sample
from spherogram.links.random_links import map_to_link, random_map

def rotate_list(L, s):
    n = len(L)
    return [L[(i + s) % n] for i in range(n)]


def flip(L):
    half = len(L) // 2
    for i in range(half):
        L[i],L[i+half] = L[i+half],L[i]


def clear_orientations(crossings):
    for c in crossings:
        c.directions.clear()
        c.sign = 0


def circular_rotate(self, n):
    """
    Rotate a tangle in a circular fashion.
    """
    tangle_copy = self.copy()
    adj = tangle_copy.adjacent
    #reverse second half
    adj[len(adj)//2:] = reversed(adj[len(adj)//2:])
    rotated_adj = rotate_list(adj, n)
    #undo reversal of second half
    rotated_adj[len(rotated_adj)//2:] = reversed(rotated_adj[len(rotated_adj)//2:])
    tangle_copy.adjacent = rotated_adj
    return tangle_copy


def add_random_crossing(self,label):
    """
    Randomly chooses position on boundary of the tangle and splits into a new
    crossing.
    """
    tangle_copy = self.copy()
    adj = tangle_copy.adjacent
    adj[len(adj)//2:] = reversed(adj[len(adj)//2:])
    new_crossing = spherogram.Crossing(label)
    old_position = randint(0,len(adj)-1)
    old_crossing, old_strand = adj.pop(old_position)
    new_strand = randint(0,3)
    old_crossing[old_strand] = new_crossing[new_strand]
    for i in range(1,4):
        adj.insert(old_position,(new_crossing,(new_strand-i) % 4))
    adj[len(adj)//2:] = reversed(adj[len(adj)//2:])
    tangle_copy.crossings.append(new_crossing)
    tangle_copy.n = self.n+1
    return tangle_copy


"""
Repeatedly splits crossings, starting with a single crossing, resulting
in a random 4-valent tree of crossings.
"""
def random_tree(size):
    T = OneTangle()
    T.crossings[0].label = 0
    for i in range(1,size):
        T = T.add_random_crossing(i)
    return T


"""
Generate two random tree tangles and glues together.  Gives a link of size
2*size
"""
def random_tree_link(size):
    return random_tree(size).circular_sum(random_tree(size),0)


def random_tree_knot(size, simplify=None, prime_decomp=False):
    found_nontrivial = False
    while not found_nontrivial:
        link = random_tree_link(size)
        knot = link.sublink(max(link.link_components,key=len))
        knot.simplify(mode=simplify)
        if len(knot) > 0:
            found_nontrivial = True

    if prime_decomp:
        cant_deconnect = False
        while cant_deconnect:
            ds = knot.deconnect_sum()
            knot = max(ds,key=len)
            cant_deconnect = (len(ds) > 1)
    return knot


"""
Glue two tangles together.  There are many ways to do this.  Choice
is given by the choice if integer n
"""
def circular_sum(self,other,n):
    if len(self.adjacent) != len(other.adjacent):
        raise Exception("Tangles do not have the same number of strands")
    return (self*(other.circular_rotate(n))).denominator_closure()


"""
All possible tangle sums as above
"""
def all_circular_sums(self,other):
    if len(self.adjacent) != len(other.adjacent):
        raise Exception("Tangles do not have the same number of strands")
    return [self.circular_sum(other,n) for n in range(len(other.adjacent))]


"""
def is_isotopic(self,other):
    if len(self.crossings) != len(other.crossings): return False
    if self.n != other.n: return False
    self_strands = self.all_cross_strands()
    self_orientations = crossing_orientations(self_strands)
    self_strands_lengths = map(len,self_strands)
    for i in range(self.n*2):
        other_strands = other.circular_rotate(i).all_cross_strands()
        other_strands_lengths = map(len,other_strands)
        if other_strands_lengths == self_strands_lengths:
            other_orientations = crossing_orientations(other_strands)
            crossing_function = []
            orientations_match = True
            print('self_orientations: ')
            print(self_orientations)
            print('other_orientations:')
            print(other_orientations)
            print('self_strands:')
            print(self_strands)
            print('other_strands:')
            print(other_strands)
            print('self.crossings:')
            print(self.crossings)
            print('other.crossings:')
            print(other.crossings)
            print('self.adjacent:')
            print(self.adjacent)
            print('other.adjacent:')
            print(other.adjacent)
            for j in range(len(self_strands)):
                for k in range(self_strands_lengths[j]):
                    cself = self_strands[j][k][0]
                    cother = other_strands[j][k][0]
                    if self_orientations[cself] != other_orientations[cother]:
                        orientations_match = False
                        break
                    crossing_function.append((cself,cother))
                if not orientations_match:
                    break
            if orientations_match:
#                print(crossing_function)
                if _is_injection(crossing_function):
                    return True #Necessarily surjective by construction
    return False
"""


def unrooted_isotopic(self, other, rotation=False):
    """
    Determine if two diagrams are isotopic to each (ignoring over/under)
    """
    if len(self.crossings) != len(other.crossings):
        return False
    if self.n != other.n:
        return False
    for i in range(2*self.n):
        if self.isosig() == other.circular_rotate(i).isosig():
            if rotation:
                return True,2*self.n-i
            else:
                return True
    return False

def random_root(self):
    return choice(choice(self.crossings).crossing_strands())


def isosig(self, root=None, over_or_under=False):
    strands, loops, orientations, crossing_order, over_or_under_data = self.all_cross_strands()
    isosig_strands = []
    isosig_loops = []
    if root is not None:
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
            root_strand_indices = [i for i in range(len(isst)) if isst[i] == root_position]
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
            root_loop_indices = [i for i in range(len(isl)) if isl[i] == root_position]
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
    return (len(self.crossings),self.n),tuple(isosig_strands),tuple(isosig_loops),tuple(isosig_orientations), tuple(isosig_over_or_under), root_isosig

def min_isosig(self,root=None,over_or_under=False):
    if root is not None:
        cs_name = cslabel(root)
    isosigs = []
    for i in range(self.n*2):
        rotated_tangle = self.circular_rotate(i)
        if root is not None:
            rotated_root = crossing_strand_from_name(rotated_tangle,cs_name)
        else:
            rotated_root = None
        isosigs.append(rotated_tangle.isosig(root=rotated_root,over_or_under=over_or_under))
    return min(isosigs)


def isosig_with_gluings(self, gluings, root=None):
    return (self.isosig(root=root), tuple(gluings))


def min_isosig_with_gluings(self, gluings, root=None):
    if root is not None:
        cs_name = cslabel(root)
    isosigs = []
    for i in range(self.n*2):
        rotated_tangle = self.circular_rotate(i)
        if root is not None:
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


"""
Given the strands, compute the orientations (+1 or -1) for each crossing
in the tangle
"""
def crossing_orientations(strands):
    orientations = {}
    over_or_under = {}
    css_seen = []
    for strand in strands:
        for cs in strand:
            for seen_cs in css_seen:
                if cs[0] == seen_cs[0]:
                    orientation = (cs[1]-seen_cs[1]) % 4
                    if orientation == 3:
                        orientation = -1
                    orientations[cs[0]] = orientation
                    over_or_under[cs[0]] = (cs[1] % 2)
                    break
            css_seen.append(cs) #didn't find cs
    return orientations, over_or_under


"""
Helper function, determines if a list of pairs defines an injection
"""
def _is_injection(pairs):
    for pair1 in pairs:
        for pair2 in pairs:
            if (pair1[0] == pair2[0]) and (pair1[1] != pair2[1]):
                return False
            if (pair1[1] == pair2[1]) and (pair1[0] != pair2[0]):
                return False
    return True


"""
Give a list of the crossing strands encountered starting at
a strand on the boundary of a tangle and moving to the other
end of that strand.
"""
def cross_strand(self, i):
    if i >= 2*self.n:
        raise Exception("Not a valid start position for strand")
    cs = self.adjacent[i]
    strand = [cs]
    while (cs[0], (cs[1] + 2) % 4) not in self.adjacent:
        cs = cs[0].adjacent[(cs[1] + 2) % 4]
        strand.append(cs)
    return strand


"""
Get the closed loop starting at crossing strand cs
"""
def loop_strand(cs):
    strand = [cs]
    cs = cs[0].adjacent[(cs[1] + 2) % 4]
    while cs not in strand:
#        print(strand)
        strand.append(cs)
        cs = cs[0].adjacent[(cs[1] + 2) % 4]
    return strand


"""
Returns all the strands but without duplicate in the opposite direction,
starting at position 0 and going clockwise, and then components that
don't intersect the boundary.
"""
def all_cross_strands(self):
    other_ends_seen = [] #use to eliminate duplicate strand in other direction
    strands = []
    strands_with_ends = []
    loops = []
    clockwise_order = range(self.n)
    clockwise_order.extend(reversed(range(self.n,self.n*2)))
    for i in clockwise_order:
        if i not in other_ends_seen:
            strand = self.cross_strand(i)
            cs = strand[-1]
            end = self.adjacent.index((cs[0],(cs[1]+2) % 4))
            if end not in other_ends_seen:
                strands.append(strand)
                strands_with_ends.append((strand,end))
                other_ends_seen.append(end)
    orientations, over_or_under = crossing_orientations(strands)
    cs_seen = [cs for strand in strands for cs in strand]
    seen_once = {cs[0] for cs in cs_seen}
    for crossing in orientations:
        seen_once.remove(crossing)
    for strand in strands:
        for cs in strand:
            if cs[0] in seen_once:
                loop = loop_strand((cs[0],(cs[1]+1) % 4))
                loops.append(loop)
                cs_seen.extend(loop)
                for loop_cs in loop:
                    if loop_cs[0] in seen_once:
                        for seen_cs in cs_seen:
                            if loop_cs[0] == seen_cs[0]:
                                orientation = (loop_cs[1]-seen_cs[1]) % 4
                                if orientation == 3:
                                    orientation = -1
                                orientations[loop_cs[0]] = orientation
                                over_or_under[loop_cs[0]] = loop_cs[1] % 2
                                seen_once.remove(loop_cs[0])
                                break
                    else:
                        seen_once.add(loop_cs[0])
    while len(orientations) < len(self.crossings):
        for loop in loops:
            for cs in loop:
                if cs[0] in seen_once:
                    loop = loop_strand((cs[0],(cs[1]+1) % 4))
                    loops.append(loop)
                    cs_seen.extend(loop)
                    for loop_cs in loop:
                        if loop_cs[0] in seen_once:
                            for seen_cs in cs_seen:
                                if loop_cs[0] == seen_cs[0]:
                                    orientation = (loop_cs[1]-seen_cs[1]) % 4
                                    if orientation == 3:
                                        orientation = -1
                                    orientations[loop_cs[0]] = orientation
                                    over_or_under[loop_cs[0]] = loop_cs[1] % 2
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


Tangle.circular_rotate = circular_rotate
Tangle.circular_sum = circular_sum
Tangle.all_circular_sums = all_circular_sums
Tangle.add_random_crossing = add_random_crossing
Tangle.unrooted_isotopic = unrooted_isotopic
Tangle.cross_strand = cross_strand
Tangle.all_cross_strands = all_cross_strands
Tangle.isosig = isosig
Tangle.min_isosig = min_isosig
Tangle.isosig_with_gluings = isosig_with_gluings
Tangle.min_isosig_with_gluings = min_isosig_with_gluings


"""
Uses networkx's cycle basis function on dual graph and converts to form with spherogram objects
"""
def cycle_basis(G):
    Gx = G.to_networkx()
    vert_cycles = nx.cycle_basis(Gx)
    return [edge_cycle(vert_cycle,G) for vert_cycle in vert_cycles]


def is_trivial(four_cycle):
    crossings = map(lambda x: map(lambda y: y.crossing,x.interface), four_cycle)
    return bool(len(set(crossings[0]) & set(crossings[1]) & set(crossings[2]) & set(crossings[3])))


def all_four_cycles_at_vertex(G, start_vertex):
    adjacent = G.children(start_vertex)
    four_cycles = set()
    for v in adjacent:
        for w in adjacent:
            if v == w:
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
    c[0] = c[1]
    c[2] = c[3]
    U = spherogram.Link([c])
    for i in range(num_attempts):
        print(i)
        Uc = U.copy()
        Uc.backtrack(backtrack_height)
        for i in range(num_mutations):
            Uc = random_mutate(Uc)
        Uc.simplify(mode='level')
        if len(Uc) > 0:
            return Uc
    return None


# Returns first nontrivial 4-cycle found
def get_four_cycle(G, start_vertex):
    adjacent = G.children(start_vertex)
    for v in adjacent:
        for w in adjacent:
            if v == w:
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


"""
Converts from list of vertices of dual graph to list of edges.  If multiple edges,
just chooses one.
"""
def edge_cycle(vert_list, G):
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


"""
Returns the crossings within distance r of the crossing, in the form
of a dictionary, where the values are the distances to the center crossing,
and a list of the crossing strands along the boundary.
"""
def crossing_ball(crossing,radius):
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


def boundary_components(link, crossing, radius):
    crossings, adjacent = crossing_ball(crossing, radius)
    crossings = list(crossings)
    G = underlying_graph(link)
    for c in crossings:
        G.remove_node(c)
    return nx.connected_components(G)


def boundary_comp_dist(samples,size,radius,edge_conn=2):
    dist = []
    for i in range(samples):
        print(i)
        K = map_to_link(random_map(size, edge_conn))
        c = choice(K.crossings)
        dist.append(num_boundary_components(K, c, radius))
    return Counter(dist)


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
            if i > 100:
                raise Exception()
        cs = cs.rotate(1)
        boundary_comp.append(cs)
    boundary_comp.pop(-1)  # code aboves adds the start_cs twice
    return boundary_comp


"""
Splits a link into two tangles along a ball around a crossing of the given
radius.  Destroys original link.
"""
def tangle_neighborhood(link,crossing,radius,return_gluings=True,hull=False):
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
        sides = {cslabel(cross_strand): cross_strand for cross_strand in adjacent}
        c = largest_comp.pop()
        cs = choice(c.crossing_strands())
        exit_strand = meander(cs,sides)[1]
        exit_strand = exit_strand[0].crossing_strands()[exit_strand[1]]
        main_boundary_comp = trace_boundary_component(exit_strand,adjacent)
        print('main_boundary_comp' + str(main_boundary_comp))
        print('all comps: '+str(comps))
        comps.remove(largest_comp) #remove largest component
        for comp in comps:
            print('crossings: ' + str(crossings))
            print('filling in comp:' + str(comp))
            print('adjacent: ' + str(adjacent))
            c = comp.pop()
            cs = choice(c.crossing_strands())

            print('cs: ' + str(cs))
            exit_strand = meander(cs,sides)[1] #meander until you hit boundary
            exit_strand = exit_strand[0].crossing_strands()[exit_strand[1]]
            print('exit_strand: ' + str(exit_strand))
            bound_comp = trace_boundary_component(exit_strand,adjacent)
            print('traced component: ' + str(bound_comp))
            if exit_strand not in main_boundary_comp:
                for x in bound_comp:
                    adjacent.remove(x)
                print('updated adjacent: ' + str(adjacent))
                crossings.append(c)
                crossings.extend(list(comp))

    adjacent[n:] = reversed(adjacent[n:])
    opposites[n:] = reversed(opposites[n:])
    gluings = []
    seen_cs = []
    #figure out which strands are glued to each other
    for cs in adjacent:
        if cs in seen_cs:
            continue
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
    # print(len(T1.crossings),len(T2.crossings))
    return T1.circular_sum(T2,2)


def mutate_reflect(link,four_cycle,other_reflection=False):
    link_copy = link.copy()
    four_cycle_copy = [corresponding_edge(link_copy,edge) for edge in four_cycle]
    if other_reflection:
        f = four_cycle_copy.pop(0)
        four_cycle_copy.append(f)
    T1,T2 = tangle_cut(link_copy,four_cycle_copy)
    crossing_labels = [c.label for c in T2.crossings]
#    print(len(T1.crossings),len(T2.crossings))
    new_link = T1.circular_sum(tangle_reflect(T2),0)
    for cl in crossing_labels:
        crossing_from_name(new_link,cl).rotate(1)
    new_link._rebuild()
    return new_link


def tangle_reflect(T):
    for c in T.crossings:
        c1adj = c.adjacent[1][:]
        c[1] = c.adjacent[3]
        c[3] = c1adj
    for i in range(T.n):
        a = T.adjacent.pop(0)
        T.adjacent.append(a)
#    [c.info() for c in T.crossings]
#    print(T.adjacent)
    return T


def all_rotation_mutants(link):
    G = link.dual_graph()
    return [mutate(link,four_cycle) for four_cycle in all_four_cycles(G)]


def all_mutants(link):
    G = link.dual_graph()
    fcs = all_four_cycles(G)
    mutants = [mutate(link,four_cycle) for four_cycle in fcs]
    mutants.extend( [mutate_reflect(link,four_cycle) for four_cycle in fcs] )
    mutants.extend( [mutate_reflect(link,four_cycle) for four_cycle in fcs] )
    return mutants


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


"""
Creates two Tangle objects from a given cycle (with no self intersections)
in the dual graph of a link (inside and outside).  Cycle is given
as a list of (oriented) edges in the dual graph. Make sure crossings are uniquely labeled. Destroys original link.
"""
def tangle_cut(link, cycle):
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

    # One of the tangles is the 'outside', and needs to be flipped
    # Just check side0
    side0_needs_flip = False
    c, i = side0[0]
    while True:
        next_cep = c.crossing_strands()[(i+1) % 4]
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


"""
Given boundary as above , fill in crossings on either side from sides. Returns
a dictionary with the side (0 or 1) of each crossing
"""
def fill_in_crossings(link,sides):
    crossing_sides = {x[0]:sides[x] for x in sides}
    crossing_labels = map(lambda c: c.label,link.crossings)
    crossings_to_sort = set(crossing_labels)-{x[0] for x in sides}
    while len(crossings_to_sort) > 0:
        start_crossing = crossings_to_sort.pop()
        accumulated_crossings = [start_crossing]
        m,end_side = meander(crossing_strand_from_name(link,(start_crossing,randint(0,3))),sides)
        accumulated_crossings.extend(map(lambda x: x.label,m))
        for c in accumulated_crossings:
            crossing_sides[c] = end_side
    return crossing_sides


"""
Wander randomly starting at crossing strand cs
until you hit a boundary strand in sides.
Assumes the crossing of cs is not on the side.
Returns a set of all the crossings encountered along the way, including ones on the boundary,
and which side (0 or 1) is hit
"""
def meander(cs, sides):
    crossings_encountered = [cs.crossing]
    end_side = 0
    while True:
        cs = cs.opposite().rotate(randint(1, 3))
        if cslabel(cs) in sides:  # hit the side
            end_side = sides[cslabel(cs)]
            break
        crossings_encountered.append(cs.crossing)
    return set(crossings_encountered), end_side


"""
Label of crossing strand, without frills
"""
def cslabel(cs):
    return (cs[0].label,cs[1])


"""
Find crossing strand object from it's name in the format of cslabel above
"""
def crossing_strand_from_name(link,csname):
    for c in link.crossings:
        if c.label == csname[0]:
            return c.crossing_strands()[csname[1]]
    raise Exception("Crossing not found")

def crossing_from_name(link,crossname):
    for c in link.crossings:
        if c.label == crossname:
            return c
    raise Exception("Crossing not found")


"""
T = spherogram.links.tangles.OneTangle()
Ts = T*T

K1 = Ts.denominator_closure()
print(Ts.adjacent)
Tsr = Ts.tangle_rotate(1)
print(Tsr.adjacent)
K2 = Tsr.denominator_closure()
"""

"""
c0,c1,c2,c3 = [spherogram.Crossing(i) for i in range(4)]
c0[0]=c1[1]
c0[1]=c1[0]
c2[3]=c3[0]
c2[0]=c3[3]
c0[2]=c2[1]
c0[3]=c3[2]
c1[3]=c2[2]
c1[2]=c3[1]
fig8 = spherogram.Link([c0,c1,c2,c3])
G8 = fig8.dual_graph()
cycle8 = cycle_basis(G8)[0]
for V in G8.vertices:
    print(V,V[:])
for e in G8.edges:
    print(e,e.interface)
T1, T2 = tangle_neighborhood(fig8,cycle8)
print('T1:')
print(T1.crossings)
print(T1.adjacent)

print('T2:')
print(T2.crossings)
print(T2.adjacent)
"""


def hamilton(G):
    F = [(G,[G.nodes()[0]])]
    n = G.number_of_nodes()
    while F:
        graph,path = F.pop()
        confs = []
        for node in graph.neighbors(path[-1]):
            conf_p = path[:]
            conf_p.append(node)
            conf_g = nx.Graph(graph)
            conf_g.remove_node(path[-1])
            confs.append((conf_g,conf_p))
        for g,p in confs:
            if len(p) == n:
                return p
            else:
                F.append((g, p))
    return None


def hamiltonian_cycle(G):
    return edge_cycle(hamilton(G.to_networkx()),G)

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


def neighborhood_distribution(link, radius):
    nhds = all_neighborhoods(link, radius)
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


from collections import Counter


def isosig_dist(num_samples, size, radius, edge_conn=2):
    nhds = []
    for i in range(num_samples):
        print(i)
        link = map_to_link(random_map(size,edge_conn))
        c = choice(link.crossings)
        root = choice(c.crossing_strands())
        T1, T2, gluings = tangle_neighborhood(link,c,radius)
        root = crossing_strand_from_name(T1,cslabel(root))
        nhds.append(T1.min_isosig(root))
    return Counter(nhds)

def isosig_dist_with_gluings(num_samples,size,radius,edge_conn=2):
    nhds = []
    for i in range(num_samples):
        print(i)
        link = map_to_link(random_map(size,edge_conn))
        c = choice(link.crossings)
        root = choice(c.crossing_strands())
        T1, T2, gluings = tangle_neighborhood(link,c,radius,hull=True)
        root = crossing_strand_from_name(T1,cslabel(root))
        nhds.append(T1.min_isosig_with_gluings(gluings,root))
    return Counter(nhds)


def unrooted_isosig_dist(num_samples,size,radius,edge_conn=2):
    nhds = []
    for i in range(num_samples):
        print(i)
        link = map_to_link(random_map(size,edge_conn))
        c = choice(link.crossings)
        T1, T2 = tangle_neighborhood(link,c,radius)
        nhds.append(T1.min_isosig())
    return Counter(nhds)

def nhd_size_dist(num_samples,link_size,radius):
    nhds = []
    for i in range(num_samples):
        print(i)
        link = map_to_link(random_map(link_size,2))
        c = link.crossings[0]
        T1, T2 = tangle_neighborhood(link,c,radius)
        nhds.append((len(T1.crossings),T1.n))
    return Counter(nhds)


def neighborhood_distribution_different_links(num_samples,size,radius):
    nhd_classes = []
    nhds = []
    for i in range(num_samples):
        print(i)
        link = map_to_link(random_map(size,2))
        c = link.crossings[0]
        T1, T2 = tangle_neighborhood(link,c,radius)
        nhds.append(T1)
    i = 0
    for nhd in nhds:
        print(i)
        i += 1
        already_found = False
        for nhd_class in nhd_classes:
            if nhd.is_isotopic(nhd_class[0]):
                nhd_class[1] += 1
                already_found = True
                break
        if not already_found:
            nhd_classes.append([nhd, 1])
    return nhd_classes


def all_neighborhood_volumes(link,radius):
    vols = []
    for c in link.crossings:
        link_copy = link.copy()
        c_copy = crossing_by_label(c.label,link_copy)
        T1, T2 = tangle_neighborhood(link_copy,c_copy,radius)
        double = T1.circular_sum(T1,0)
        v = double.exterior().volume()
        vols.append((double,v))
    return vols


def all_nhd_vol_dists(link,radius,tolerance):
    nhds = []
    for c in link.crossings:
        link_copy = link.copy()
        c_copy = crossing_by_label(c.label,link_copy)
        T1, T2 = tangle_neighborhood(link_copy,c_copy,radius)
        double_volumes = map(lambda x: x.exterior().volume(),T1.all_circular_sums(T1))
        double_volumes.sort()
        for i in range(len(double_volumes)):
            if double_volumes[i] < tolerance:
                double_volumes[i] = 0
        new_dist = True
        for dist in nhds:
            if close_float_sets(double_volumes,dist[0],tolerance):
                new_dist = False
                dist[1] += 1
                break
        if new_dist:
            nhds.append([double_volumes,1,T1])
    return nhds


def close_float_sets(L1, L2, tolerance):
    """
    Asks if two sorted lists of floats are near each other
    """
    if len(L1) != len(L2):
        return False
    for i in range(len(L1)):
        if abs(L1[i] - L2[i]) > tolerance:
            return False
    return True


K = map_to_link(random_map(100,2))
Kc = K.copy()
c = Kc.crossings[0]
v = Kc.exterior().volume()
G = K.dual_graph()
cycle = max(cycle_basis(G),key=len)
T1, T2, gluings = tangle_neighborhood(Kc,c,1)
print(gluings)
"""
v1 = T1.circular_sum(T1,0).exterior().volume()
v2 = T2.circular_sum(T2,0).exterior().volume()
avg_vol = v1/2+v2/2
tsums = T1.all_circular_sums(T2)
vols = map(lambda x: x.exterior().volume(), tsums)

print(v,avg_vol)
"""

a, b, c, d, e, f, g, h = (spherogram.Crossing(x) for x in 'abcdefgh')
a[0] = e[2]
a[1] = b[3]
a[3] = e[3]
b[0] = f[2]
b[1] = c[3]
b[2] = c[2]
c[0] = g[2]
c[1] = d[3]
d[0] = h[2]
d[1] = h[1]
e[1] = f[3]
f[0] = g[0]
f[1] = g[3]
g[1] = h[3]
crossings = [a, b, c, d, e, f, g, h]
adj = [(e, 0), (h, 0), (a, 2), (d, 2)]
testT = spherogram.Tangle(2, crossings, adj)
