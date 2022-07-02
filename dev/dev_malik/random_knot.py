from spherogram import Crossing, Link
from spherogram.graphs import FatGraph, CyclicList
from random import choice, randint
import networkx as nx
from itertools import combinations, product
import pickle    


def start_string():
    crossings = [Crossing('0')]
    crossings[0][2]=crossings[0][3]
    final_cs = crossings[0].crossing_strands()[0]
    loose_cs = crossings[0].crossing_strands()[1]
    return crossings, loose_cs, final_cs


def random_knot(n, method='close_under', alternate=False,
                bias=False, num_opportunities=1):
    if method == 'close_under':
        crossings, loose_cs, final_cs = random_open_string(n, alternate=alternate, bias=bias)
        connect_loose_strands(crossings, loose_cs, final_cs)
        return Link(crossings)
    elif method == 'close_on_opportunity':
        crossings = [Crossing('0')]
        crossings[0][2]=crossings[0][3]
        final_cs = crossings[0].crossing_strands()[0]
        loose_cs = crossings[0].crossing_strands()[1]
        i = 0
        while num_opportunities > 0:
            available = available_strands(loose_cs)
            strand_to_cross = choice(available)
            if alternate:
                over_or_under = 1-(i%2) 
            else:
                over_or_under = randint(0,1)
            loose_cs = cross_strand(crossings, loose_cs, 
                                strand_to_cross, str(i+1), over_or_under)
            same_face = set(available_strands(loose_cs)) == set(available_strands(final_cs))
            if same_face:
                num_opportunities -= 1
            i += 1
            if i >= n:
                raise Exception('Not closeable after n crossings')
        connect_loose_strands(crossings, loose_cs, final_cs)
        return Link(crossings)                


def random_open_string(n, alternate=False, bias=False):
    crossings = [Crossing('0')]
    crossings[0][2]=crossings[0][3]
    final_cs = crossings[0].crossing_strands()[0]
    loose_cs = crossings[0].crossing_strands()[1]
    for i in range(n):
        available = available_strands(loose_cs)
        if bias:
            strand_to_cross = choice(bias_middle(available))
        else:
            strand_to_cross = choice(available)
        if alternate:
            over_or_under = 1-(i%2)            
        else:
            over_or_under = randint(0,1)
        loose_cs = cross_strand(crossings, loose_cs, 
                                strand_to_cross, str(i+1), over_or_under) 
    return crossings, loose_cs, final_cs

def bias_middle(start_list):
    biased_list = []
    num_copies = range(len(start_list)/2)
    if len(start_list)%2 == 1:
        num_copies.append(len(start_list)/2)
    num_copies.extend(reversed(range(len(start_list)/2)))
    for i, x in enumerate(num_copies):
        for n in range(x+1):
            biased_list.append(start_list[i])
    return biased_list

def distance_drift(n):
    crossings, loose_cs, final_cs = start_string()
    dists = []
    for i in range(n):
        print(i)
        available = available_strands(loose_cs)
        strand_to_cross = choice(available)
        over_or_under = randint(0,1)
        loose_cs = cross_strand(crossings, loose_cs, 
                                strand_to_cross, str(i+1), over_or_under)
        dists.append(distance(crossings, loose_cs, final_cs))
    return dists

def distance(crossings, loose_cs, final_cs):
    G, loose_face, final_face = open_string_dual_graph(crossings, loose_cs, 
                                                       final_cs)
    path = nx.shortest_path(G,loose_face,final_face)
    return len(path)

def volume_evolution(n):
    crossings = [Crossing('0')]
    crossings[0][2]=crossings[0][3]
    final_cs = crossings[0].crossing_strands()[0]
    loose_cs = crossings[0].crossing_strands()[1]
    vols = []
    for i in range(n):
        print(i)
        available = available_strands(loose_cs)
        if final_cs in available:
            available.remove(final_cs)
        strand_to_cross = choice(available)
        loose_cs = cross_strand(crossings, loose_cs, 
                                strand_to_cross, str(i+1), randint(0,1))
        vols.append(open_string_volume(crossings))
    return vols


def open_string_volume(crossings):
    crossings_copy = pickle.loads(pickle.dumps(crossings))
    loose = []
    for c in crossings_copy:
        for i in range(4):
            if c.adjacent[i] is None:
                loose.append(c.crossing_strands()[i])
    connect_loose_strands(crossings_copy, loose[0], loose[1])
    return Link(crossings_copy).exterior().volume()


def alexander_evolution(n):
    crossings = [Crossing('0')]
    crossings[0][2]=crossings[0][3]
    final_cs = crossings[0].crossing_strands()[0]
    loose_cs = crossings[0].crossing_strands()[1]
    alex_polys = []
    for i in range(n):
        print(i)
        available = available_strands(loose_cs)
        if final_cs in available:
            available.remove(final_cs)
        strand_to_cross = choice(available)
        loose_cs = cross_strand(crossings, loose_cs, 
                                strand_to_cross, str(i+1), randint(0,1))
        alex_polys.append(open_string_alexander(crossings))
    return alex_polys


def function_evolution(n, function, simplify='level', skip_first=0):
    crossings = [Crossing('0')]
    crossings[0][2]=crossings[0][3]
    final_cs = crossings[0].crossing_strands()[0]
    loose_cs = crossings[0].crossing_strands()[1]
    values = []
    for i in range(n):
        print(i)
        available = available_strands(loose_cs)
        if final_cs in available:
            available.remove(final_cs)
        strand_to_cross = choice(available)
        loose_cs = cross_strand(crossings, loose_cs, 
                                strand_to_cross, str(i+1), randint(0,1))
        if i >= skip_first:
            values.append(open_string_evaluation(crossings, function, simplify))
    return values

def turn_list_from_open_string(crossings, loose_cs, final_cs):
    turn_list = []
    while len(crossings) > 1:
#        print(turn_list)
        old_position_crossing, old_position_index = loose_cs.crossing.adjacent[(loose_cs.strand_index-1)%4]
        old_crossing_sign = loose_cs.strand_index%2
        if old_crossing_sign == 0:
            old_crossing_sign = -1
        else:
            old_crossing_sign = 1
        old_cs_position = old_position_crossing.crossing_strands()[old_position_index]
        loose_cs = go_back(crossings, loose_cs, final_cs)
        turn_index = available_strands(loose_cs).index(old_cs_position)+1
        turn_list.insert(0,turn_index*old_crossing_sign)
    return turn_list


def fix_turn_list(turn_list):
    return turn_list_from_open_string(*open_string_from_turn_list(turn_list))


def function_evolution_deterministic(n, function, simplify='level',
                                     start_turn_list=[1]):
    crossings, loose_cs, final_cs = open_string_from_turn_list(start_turn_list)
    turn_list = start_turn_list
    values = []
    for i in range(n):
        print(i, turn_list)
        values.append(open_string_evaluation(crossings, function, simplify))
        crossings, loose_cs, final_cs, turn_list = next_open_string_with_over_under(crossings, loose_cs, final_cs, turn_list)
    return values


def trivial_jones_search(n, simplify='global', start_turn_list=[1]):
    crossings, loose_cs, final_cs = open_string_from_turn_list(start_turn_list)
    turn_list = start_turn_list
    trivial_jones = []
    for i in range(n):
        K = open_string_evaluation(crossings, lambda K: K, simplify)
        if K:
            if len(K)>17:
                p = K.jones_polynomial()
                if p.is_one():
                    trivial_jones.append(K)
        print(i,turn_list, K)
        crossings, loose_cs, final_cs, turn_list = next_open_string_with_over_under(crossings, loose_cs, final_cs, turn_list)
    return trivial_jones


def knot_search(turn_list, max_steps):
    crossings, loose_cs, final_cs = open_string_from_turn_list(turn_list)
    knots = []
    for i in range(max_steps):
        print(turn_list)
        crossings, loose_cs, final_cs, turn_list = next_open_string(crossings,
                                                                    loose_cs,
                                                                    final_cs,
                                                                    turn_list)
        if set(available_strands(loose_cs)) == set(available_strands(final_cs)):
            knots.append(tuple(turn_list))
    return knots


def open_string_evaluation(crossings, function, simplify='level'):
    crossings_copy = pickle.loads(pickle.dumps(crossings))
    loose = []
    for c in crossings_copy:
        for i in range(4):
            if c.adjacent[i] is None:
                loose.append(c.crossing_strands()[i])
    connect_loose_strands(crossings_copy, loose[0], loose[1])
    K = Link(crossings_copy)
    if simplify:
        K.simplify(mode=simplify)
    if len(K) > 0:
        return function(K)


def open_string_alexander(crossings):
    crossings_copy = pickle.loads(pickle.dumps(crossings))
    loose = []
    for c in crossings_copy:
        for i in range(4):
            if c.adjacent[i] is None:
                loose.append(c.crossing_strands()[i])
    connect_loose_strands(crossings_copy, loose[0], loose[1])
    K = Link(crossings_copy)
    K.simplify(mode='level')
    return K.alexander_polynomial() if len(K) > 0 else 1


def open_string_from_turn_list(turn_list):
    crossings = [Crossing('0')]
    crossings[0][2]=crossings[0][3]
    final_cs = crossings[0].crossing_strands()[0]
    loose_cs = crossings[0].crossing_strands()[1]
    for i, turn in enumerate(turn_list):
        available = available_strands(loose_cs)
#        print('available')
#        print(available)
#        print(loose_cs, final_cs)
#        print('crossings before')
#        [c.info() for c in crossings]
        if final_cs in available:
            available.remove(final_cs)
        strand_to_cross = available[(abs(turn)-1) % len(available)]
 #       print('strand_to_cross')
#        print(strand_to_cross)

#        print('available after')
#        print(available)
        loose_cs = cross_strand(crossings, loose_cs, 
                                strand_to_cross, str(i+1), (turn>0)) 
#        print('crossings after')
#        [c.info() for c in crossings]

    return crossings, loose_cs, final_cs

def knot_tree(depth):
    crossings = [Crossing('0')]
    crossings[0][2]=crossings[0][3]
    final_cs = crossings[0].crossing_strands()[0]
    loose_cs = crossings[0].crossing_strands()[1]
    tree = nx.Graph()
    tree.add_node('')
    leaves = ['']
    for i in range(depth):
        print(i)
        leaves = split_leaves(tree, leaves)
    return tree, [map(int, leaf.split()) for leaf in leaves]

def split_leaves(tree, leaves):
    new_leaves = []
    for leaf in leaves:
        turn_list = map(int, leaf.split())
        crossings, loose_cs, final_cs = open_string_from_turn_list(turn_list)
        available = available_strands(loose_cs)
        for i in range(1, len(available)+1):
            new_leaf = leaf + ' ' + str(i)
            tree.add_edge(leaf, new_leaf)
            new_leaves.append(new_leaf)

    return new_leaves


def all_possible_next_crossings(crossings, loose_cs, final_cs):
    available = available_strands(loose_cs)
    open_strands = []
    for strand_to_cross in available:
        loose_cs = cross_strand(crossings, loose_cs, 
                                strand_to_cross, str(i+1), 0)

    return crossings, loose_cs, final_cs

def go_back(crossings, loose_cs, final_cs):
    new_loose_cs = loose_cs.rotate(2).opposite()
    new_loose_cs.crossing.adjacent[new_loose_cs.strand_index] = None
    crossings.remove(loose_cs.crossing)
    connect_crossing_strands(loose_cs.rotate(1).opposite(), loose_cs.rotate(3).opposite())
    return new_loose_cs

def next_open_string_with_over_under(crossings, loose_cs, final_cs, turn_list):
    can_move_right = False
    can_switch_to_under = (turn_list[-1] > 0)
    start_size = len(crossings)
    while not can_move_right or can_switch_to_under:
#        print(turn_list, len(crossings))
#        print(crossings)
        if len(crossings) == 1:
            turn_list = [1]*start_size
            crossings, loose_cs, final_cs = open_string_from_turn_list(turn_list)
            return crossings, loose_cs, final_cs, turn_list
        loose_cs = go_back(crossings, loose_cs, final_cs)
        turn_list_end = turn_list.pop()
        last_strand_position = abs(turn_list_end)
        can_switch_to_under = (turn_list_end > 0)
        available = available_strands(loose_cs)
        if can_switch_to_under:
            strand_to_cross = available[last_strand_position-1]
            loose_cs = cross_strand(crossings, loose_cs, strand_to_cross, 0, 1)
            turn_list.append(-last_strand_position)
            break
        if last_strand_position < len(available):
            can_move_right = True
            strand_to_cross = available[last_strand_position]
            loose_cs = cross_strand(crossings, loose_cs, strand_to_cross, 0, 0)
            turn_list.append(last_strand_position+1)
    while len(crossings) < start_size:
        available = available_strands(loose_cs)
        strand_to_cross = available[0]
        loose_cs = cross_strand(crossings, loose_cs, strand_to_cross, 0, 0)
        turn_list.append(1)
    return crossings, loose_cs, final_cs, turn_list


def next_open_string(crossings, loose_cs, final_cs, turn_list):
    can_move_right = False
    start_size = len(crossings)
    while not can_move_right:
#        print(turn_list, len(crossings))
#        print(crossings)
        if len(crossings) == 1:
            turn_list = [1]*start_size
            crossings, loose_cs, final_cs = open_string_from_turn_list(turn_list)
            return crossings, loose_cs, final_cs, turn_list
        loose_cs = go_back(crossings, loose_cs, final_cs)
        last_strand_position = turn_list.pop()
        available = available_strands(loose_cs)
        if last_strand_position < len(available):
            can_move_right = True
            strand_to_cross = available[last_strand_position]
            loose_cs = cross_strand(crossings, loose_cs, strand_to_cross, 0, 0)
            turn_list.append(last_strand_position+1)
    while len(crossings) < start_size:
        available = available_strands(loose_cs)
        strand_to_cross = available[0]
        loose_cs = cross_strand(crossings, loose_cs, strand_to_cross, 0, 0)
        turn_list.append(1)
    return crossings, loose_cs, final_cs, turn_list

def hard_unknot_search(turn_list, max_steps):
    crossings, loose_cs, final_cs = open_string_from_turn_list(turn_list)
    unknots = []
    i = 0
    while i < max_steps:
        i += 1
        print(i)
        crossings, loose_cs, final_cs, turn_list = next_open_string_with_over_under(crossings, loose_cs, final_cs, turn_list)
        if set(available_strands(loose_cs)) == set(available_strands(final_cs)):
            print(turn_list)
            K = open_string_evaluation(crossings, lambda x: x, simplify=None)
            K.simplify()
            if len(K)>3:
                K.simplify(mode='global')
                if len(K)>3:
                    if K.alexander_polynomial().is_one():
                        print(turn_list)
                        unknots.append(tuple(turn_list))
            
    return unknots


def knot_search(turn_list, max_steps):
    crossings, loose_cs, final_cs = open_string_from_turn_list(turn_list)
    knots = []
    for i in range(max_steps):
        print(turn_list)
        crossings, loose_cs, final_cs, turn_list = next_open_string(crossings,
                                                                    loose_cs,
                                                                    final_cs,
                                                                    turn_list)
        if set(available_strands(loose_cs)) == set(available_strands(final_cs)):
            knots.append(tuple(turn_list))
    return knots


def knot_search_by_num_crossings(num_crossings):
    turn_list = [1]*(num_crossings-1)
    crossings, loose_cs, final_cs = open_string_from_turn_list(turn_list)
    knots = []
    while len(crossings) == num_crossings:
        crossings, loose_cs, final_cs, turn_list = next_open_string(crossings,
                                                                    loose_cs,
                                                                    final_cs,
                                                                    turn_list)
        if set(available_strands(loose_cs)) == set(available_strands(final_cs)):
            knots.append(tuple(turn_list))
    return knots

def all_possible_crossing_switching(knot):
    all_crossing_switches = product((0,1), repeat=len(knot))
    switched_knots = []
    for crossing_switch in all_crossing_switches:
        K = knot.copy()
        for i, c in enumerate(K.crossings):
            c.rotate(crossing_switch[i])
        K._rebuild()
        switched_knots.append(K)
    return switched_knots

def all_possible_crossing_switching_in_place(knot, func_to_evaluate):
    all_crossing_switches = product((0,1), repeat=len(knot))
    switched_knots_evaluations = []
    for crossing_switch in all_crossing_switches:
        for i, c in enumerate(knot.crossings):
            c.rotate(crossing_switch[i])
        knot._rebuild()
        switched_knots_evaluations.append(func_to_evaluate(knot))
        for i, c in enumerate(knot.crossings):
            c.rotate(crossing_switch[i])
        knot._rebuild()

    return switched_knots_evaluations


def build_knot_database(filename):
    f = open(filename, 'a+')
    turn_list = eval(f.readlines()[-1])
    crossings, loose_cs, final_cs = open_string_from_turn_list(turn_list)
    num_crossings = len(crossings)
    while len(crossings) == num_crossings:
        crossings, loose_cs, final_cs, turn_list = next_open_string(crossings,
                                                                    loose_cs,
                                                                    final_cs,
                                                                    turn_list)
        if set(available_strands(loose_cs)) == set(available_strands(final_cs)):
            f.write(str(turn_list) + '\n')
    f.close()

def knot_from_turn_list(turn_list):
    crossings, loose_cs, final_cs = open_string_from_turn_list(turn_list)
#    randomly_flip_crossings(crossings, loose_cs, final_cs)
#    print(open_string_faces(crossings, loose_cs, final_cs))
    connect_loose_strands(crossings, loose_cs, final_cs)
    return Link(crossings)

#def randomly_flip_crossings(crossings, loose_cs, final_cs):
#    for c in crossings:
#        if c != loose_cs.crossing and c != final_cs.crossing:
#            c.rotate(randint(0,1))

def cross_strand(crossings, loose_cs, strand_to_cross, new_label, over_vs_under):
    opposite = strand_to_cross.opposite()
    new_crossing = Crossing(new_label)
    crossings.append(new_crossing)
#    cs0, cs1, cs2, cs3 = new_crossing.crossing_strands()
#    connect_crossing_strands(cs0, loose_cs)
#    connect_crossing_strands(cs1, strand_to_cross)
#    connect_crossing_strands(cs3, opposite)

    css = new_crossing.crossing_strands()
    connect_crossing_strands(css[0+over_vs_under], loose_cs)
    connect_crossing_strands(css[1+over_vs_under], strand_to_cross)
    connect_crossing_strands(css[(3+over_vs_under)%4], opposite)

    return css[2+over_vs_under]


def connect_crossing_strands(cs1, cs2):
    cs1[0][cs1[1]] = cs2[0][cs2[1]]


def available_strands(loose_cs):
    start_cs = next_available_corner(loose_cs)
    available = [start_cs]
    next_cs = next_available_corner(start_cs)
    while next_cs != start_cs:
        available.append(next_cs)
        next_cs = next_available_corner(next_cs)
    return available


def next_available_corner(cs):
    rotated = cs.rotate()
    if rotated[0].adjacent[rotated[1]] is None:
        rotated_twice = rotated.rotate()
        if rotated_twice[0].adjacent[rotated_twice[1]] is None:
            return rotated_twice.rotate().opposite()
        else:
            return rotated_twice.opposite()
    else:
        return rotated.opposite()


def connect_loose_strands(crossings, loose_cs, final_cs):
    G, loose_face, final_face = open_string_dual_graph(crossings, loose_cs, 
                                                       final_cs)
    path = nx.shortest_path(G,loose_face,final_face)
    for i in range(len(path)-1):
        face = path[i]
        next_face = path[i+1]
        edges = edges_between(face, next_face)
        strand_to_cross = edges[0][0]
        loose_cs = cross_strand(crossings, loose_cs, strand_to_cross, 'u'+str(i), 0)
    connect_crossing_strands(loose_cs, final_cs)

def open_string_faces(crossings, loose_cs, final_cs):
    crossing_strands = set([cs for c in crossings for cs in c.crossing_strands()])
    faces = []
    crossing_strands.remove(loose_cs)
    loose_face = available_strands(loose_cs)
    crossing_strands.remove(final_cs)
    final_face = available_strands(final_cs)
    while crossing_strands:
        cs = crossing_strands.pop()
        face = available_strands(cs)
        faces.append(tuple(face))
        crossing_strands = crossing_strands - set(face)
    for face in faces:
        if set(face) == set(loose_face):
            loose_face = face
        if set(face) == set(final_face):
            final_face = face
    return faces, loose_face, final_face

def open_string_dual_graph(crossings, loose_cs, final_cs):
    faces, loose_face, final_face = open_string_faces(crossings, loose_cs,
                                                       final_cs)
    G = nx.Graph()
    for face in faces:
        for other_face in faces:
            if edges_between(face, other_face):
                G.add_edge(face,other_face)    
    return G, loose_face,  final_face

def edges_between(face1, face2):
    edges = []
    for cs in face1:
        if cs.opposite() in face2:
            edges.append([cs, cs.opposite()])
    return edges



"""
def random_dual_graph(n):
#    G = FatGraph()
#    G.add_edge( ('0',0) , ('1',0) )
#    current_face = '1'
#    current_corner = G.incidence_dict[current_face][0]
#    face_with_loose_edge = '0'

    G = FatGraph()
    G.add_edge( ('0',0) , ('1',0) )
    G.add_edge( ('0',1) , ('2',0) )
    G.add_edge( ('0',2) , ('3',1) )
    G.add_edge( ('2',1) , ('1',1) )
    G.add_edge( ('2',2) , ('3',0) )
    current_face = '0'
    current_corner = G.incidence_dict[current_face][1]

    for i in range(n):
        print('Vertex ordering')
        print(adjacent_vertices_in_order(G))
        next_edge = choice(list(G.incident(current_face)))
        print('next_edge')
        print(next_edge)
        print('current_face, current_corner')
        print(current_face, current_corner)
        current_face, current_corner = split_face(G, next_edge, 
                                                  current_corner, current_face)

    return G

def split_face(G, next_edge, current_corner, current_face):
    new_0_face = current_face + '0'
    new_1_face = current_face + '1'

    G.add_vertex(new_0_face)
    G.add_vertex(new_1_face)
    print('Vertices0')
    print(adjacent_vertices_in_order(G))

    current_face_edges = G.incidence_dict[current_face]

    vertex_1_edges = []
    while current_corner != next_edge:
        vertex_1_edges.append(current_corner)
        current_corner = current_face_edges.succ(current_corner)
    vertex_0_edges = [e for e in current_face_edges if e not in vertex_1_edges]
    vertex_0_edges.remove(next_edge)

#    min_0_slot = 0
    for e in vertex_0_edges:
        current_face_side = e.index(current_face)
        current_face_slot = e.slots[current_face_side]
#        if current_face_slot < min_0_slot:
#            min_0_slot = current_face_slot
        other_face, other_face_slots = e[1-current_face_side], e.slots[1-current_face_side]
        G.add_edge( (new_0_face,current_face_slot) , (other_face,other_face_slots) )
    print('Vertices1')
    print(adjacent_vertices_in_order(G))

#    min_1_slot = 0
    for e in vertex_1_edges:
        current_face_side = e.index(current_face)
        current_face_slot = e.slots[current_face_side]
#        if current_face_slot < min_1_slot:
#            min_1_slot = current_face_slot
        other_face, other_face_slots = e[1-current_face_side], e.slots[1-current_face_side]
        G.add_edge( (new_1_face,current_face_slot) , (other_face,other_face_slots) )
    print('Vertices2')
    print(adjacent_vertices_in_order(G))
    other_face = next_edge[1-next_edge.index(current_face)]
#    other_slot = next_edge.slots[1-next_edge.index(current_face)]
#   The edge immediately after where triangle will be inserted
    other_edge = G.incidence_dict[other_face].succ(next_edge)

#    print('Slots before')
#    print(G.incidence_dict)
#    print(all_occupied_slots(G))
#    print('min0: ' +str(min_0_slot) + ' ' + str(new_0_face))
#    print('min1: ' +str(min_1_slot) + ' ' + str(new_1_face))
#    print('other: ' +str(other_slot) + ' ' + str(other_face))

#    make_room(G, other_face, other_slot)
#    make_room(G, new_0_face, min_0_slot-2)
#    make_room(G, new_1_face, min_1_slot-2)

#    print('Slots after')
#    print(G.incidence_dict)
#    print(all_occupied_slots(G))
#    print('min0: ' +str(min_0_slot) + ' ' + str(new_0_face))
#    print('min1: ' +str(min_1_slot) + ' ' + str(new_1_face))
#    print('other: ' +str(other_slot) + ' ' + str(other_face))

    print('new edges 0')
    print(vertex_0_edges)
    print(G.incidence_dict[new_0_face])
    print('new edges 1')
    print(G.incidence_dict[new_1_face])
    print(vertex_1_edges)

    G.remove_vertex(current_face)

    cyclify(G)

    #Traverse around face to find where to connect new_0_face and new_1_face
    #to other_face
    if G.incidence_dict[new_0_face]:
        seeking_0_edge, seeking_0_face = (other_edge, other_face)
        while new_0_face != seeking_0_face:
            seeking_0_edge, seeking_0_face = next_corner(G, seeking_0_edge, seeking_0_face)
        first_0_edge = seeking_0_edge
    else:
        first_0_edge = None

    if G.incidence_dict[new_1_face]: #other direction
        pred_edge = G.incidence_dict[other_face].pred(other_edge)
        seeking_1_edge, seeking_1_face = (pred_edge, other_face)
        while new_1_face != seeking_1_face:
            seeking_1_edge, seeking_1_face = next_corner(G, seeking_1_edge, seeking_1_face)
        first_1_edge = G.incidence_dict[new_1_face].succ(seeking_1_edge)
    else:
        first_1_edge = None

    print('first_0_edge')
    print(first_0_edge)
    print('first_1_edge')
    print(first_1_edge)
    print('other_edge')
    print(other_edge)
    print(G.incidence_dict)

 #   G.add_edge( (new_0_face,min_0_slot-2) , (new_1_face,min_1_slot-1) )
    insert_edge(G, new_0_face, first_0_edge, new_1_face, first_1_edge)
    if first_0_edge is None:
        first_0_edge = G.edges_between(new_0_face, new_1_face).pop()
    if first_1_edge is None:
        first_1_edge = G.edges_between(new_0_face, new_1_face).pop()

    print('Vertices3')
    print(adjacent_vertices_in_order(G))
    print(G.incidence_dict)

#    G.add_edge( (new_0_face,min_0_slot-1) , (other_face,other_slot) )
    insert_edge(G, new_0_face, first_0_edge, other_face, other_edge)
    print('Vertices4')
    print(adjacent_vertices_in_order(G))
    print(G.incidence_dict)

#    G.add_edge( (new_1_face,min_1_slot-2) , (other_face,other_slot+1) )
    edge_between_0_and_1 = G.edges_between(new_0_face,new_1_face).pop()
    insert_edge(G, new_1_face, edge_between_0_and_1, other_face, other_edge)
    print('Vertices5')
    print(adjacent_vertices_in_order(G))



    cyclify(G)
    new_next_corner = G.incidence_dict[new_1_face].pred(edge_between_0_and_1)

    return other_face, new_next_corner

def edge_to_close(face_list):
    t1, t2 = filter(lambda f: len(f) == 3, face_list)
    t1 = set([x[0] for x in t1])
    t2 = set([x[0] for x in t2])
    return t1.intersection(t2)

def close_up(G):
    s = edge_to_close(faces(G))
    if s:
        G.remove_edge(s.pop())
    else:
        raise Exception("Not closeable")
        

def knot(G):
    face_list = faces(G)
    edge_labels = {e: i for i,e in enumerate(G.edges)}
    print(edge_labels)
    crossings = [Crossing(i) for i in range(len(face_list))]
    PD = [ [edge_labels[e[0]] for e in face] for face in face_list]
    return PD

def cyclify(G):
    for v in G.incidence_dict:
        G.incidence_dict[v] = CyclicList(G.incidence_dict[v])

def insert_edge(G, vertex, edge, other_vertex, other_edge):

    if G.incidence_dict[vertex]:
        vertex_insertion_point = edge.slots[edge.index(vertex)]
        move_right = 0
        for e in G.incidence_dict[vertex]:
            if e == edge:
                move_right = 1
            e.slots[e.index(vertex)] += move_right
    else:
        vertex_insertion_point = 0

    if G.incidence_dict[other_vertex]:
        other_vertex_insertion_point = other_edge.slots[other_edge.index(other_vertex)]
        move_right = 0
        for e in G.incidence_dict[other_vertex]:
            if e == other_edge:
                move_right = 1
            e.slots[e.index(other_vertex)] += move_right
    else:
        other_vertex_insertion_point = 0

    G.add_edge( (vertex, vertex_insertion_point), (other_vertex, other_vertex_insertion_point))


def make_room(G, vertex, insertion_edge):
    for edge in G.incident(vertex):
        if edge.slots[edge.index(vertex)] >= insertion_point:
            edge.slots[edge.index(vertex)] += 1

def get_edge_at_slot(G, face, slot):
    for edge in G.incident(face):
        if edge.slots[edge.index(face)] == slot:
            return edge


def next_corner(G, edge, head):
    next_edge = G.incidence_dict[head].succ(edge)
    next_head = next_edge[1-next_edge.index(head)]
    return next_edge, next_head

def faces(G):
    diedges = set([(edge, head) for edge in G.edges for head in edge])
    face_list = []
    while diedges:
        start_edge, start_head = diedges.pop()
        face = [(start_edge, start_head)]
        edge, head = next_corner(G, start_edge, start_head)
        while (edge, head) != (start_edge, start_head):
            face.append((edge,head))
            diedges.remove((edge,head))
            edge, head = next_corner(G,edge, head)
        face_list.append(face)
    return face_list

def occupied_slots(G,v):
    return [e.slots[e.index(v)] for e in G.incidence_dict[v]]

def all_occupied_slots(G):
    return [(v,occupied_slots(G,v)) for v in G.vertices]

def has_overfilled_slot(G):
    slots = all_occupied_slots(G)
    for v,s in slots:
        if len(set(s)) < len(s):
            return True
    return False

def adjacent_vertices_in_order(G):
    l = []
    for v in G.incidence_dict:
        l.append((v,[e[1-e.index(v)] for e in G.incidence_dict[v]]))
    return l


def old_split_face(G, next_edge, current_corner, current_face):
    new_0_face = current_face + '0'
    new_1_face = current_face + '1'

    G.add_vertex(new_0_face)
    G.add_vertex(new_1_face)
    print('Vertices0')
    print(adjacent_vertices_in_order(G))

    current_face_edges = G.incidence_dict[current_face]

    vertex_1_edges = []
    while current_corner != next_edge:
        vertex_1_edges.append(current_corner)
        current_corner = current_face_edges.succ(current_corner)
    vertex_0_edges = [e for e in current_face_edges if e not in vertex_1_edges]
    vertex_0_edges.remove(next_edge)

    min_0_slot = 0
    for e in vertex_0_edges:
        current_face_side = e.index(current_face)
        current_face_slot = e.slots[current_face_side]
        if current_face_slot < min_0_slot:
            min_0_slot = current_face_slot
        other_face, other_face_slots = e[1-current_face_side], e.slots[1-current_face_side]
        G.add_edge( (new_0_face,current_face_slot) , (other_face,other_face_slots) )
    print('Vertices1')
    print(adjacent_vertices_in_order(G))

    min_1_slot = 0
    for e in vertex_1_edges:
        current_face_side = e.index(current_face)
        current_face_slot = e.slots[current_face_side]
        if current_face_slot < min_1_slot:
            min_1_slot = current_face_slot
        other_face, other_face_slots = e[1-current_face_side], e.slots[1-current_face_side]
        G.add_edge( (new_1_face,current_face_slot) , (other_face,other_face_slots) )
    print('Vertices2')
    print(adjacent_vertices_in_order(G))
    other_face = next_edge[1-next_edge.index(current_face)]
    other_slot = next_edge.slots[1-next_edge.index(current_face)]

    print('Slots before')
    print(G.incidence_dict)
    print(all_occupied_slots(G))
    print('min0: ' +str(min_0_slot) + ' ' + str(new_0_face))
    print('min1: ' +str(min_1_slot) + ' ' + str(new_1_face))
    print('other: ' +str(other_slot) + ' ' + str(other_face))

    make_room(G, other_face, other_slot)
    make_room(G, new_0_face, min_0_slot-2)
    make_room(G, new_1_face, min_1_slot-2)

    print('Slots after')
    print(G.incidence_dict)
    print(all_occupied_slots(G))
    print('min0: ' +str(min_0_slot) + ' ' + str(new_0_face))
    print('min1: ' +str(min_1_slot) + ' ' + str(new_1_face))
    print('other: ' +str(other_slot) + ' ' + str(other_face))


    G.add_edge( (new_0_face,min_0_slot-2) , (new_1_face,min_1_slot-1) )
    print('Vertices3')
    print(adjacent_vertices_in_order(G))
    print(G.incidence_dict)

    G.add_edge( (new_0_face,min_0_slot-1) , (other_face,other_slot) )
    print('Vertices4')
    print(adjacent_vertices_in_order(G))
    print(G.incidence_dict)

    G.add_edge( (new_1_face,min_1_slot-2) , (other_face,other_slot+1) )
    print('Vertices5')
    print(adjacent_vertices_in_order(G))

    G.remove_vertex(current_face)
    new_next_corner = get_edge_at_slot(G, other_face, other_slot)
    
    cyclify(G)
    return other_face, new_next_corner
"""
