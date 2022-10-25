from spherogram import Crossing, Link
from spherogram.graphs import FatGraph, CyclicList
from random import choice, randint
import networkx as nx
from itertools import combinations


def random_knot(n):
    crossings, loose_cs, final_cs = random_open_string(n)
    flip_crossings(crossings, loose_cs, final_cs)
#    print(open_string_faces(crossings, loose_cs, final_cs))
    connect_loose_strands(crossings, loose_cs, final_cs)
    return Link(crossings)

def random_open_string(n):
    crossings = [Crossing('0')]
    crossings[0][2]=crossings[0][3]
    final_cs = crossings[0].crossing_strands()[0]
    loose_cs = crossings[0].crossing_strands()[1]
    for i in range(n):
        available = available_strands(loose_cs)
#        print('available')
#        print(available)
#        print(loose_cs, final_cs)
#        print('crossings before')
#        [c.info() for c in crossings]
        if final_cs in available:
            available.remove(final_cs)
        strand_to_cross = choice(available)
 #       print('strand_to_cross')
#        print(strand_to_cross)

#        print('available after')
#        print(available)
        loose_cs = cross_strand(crossings, loose_cs,
                                           strand_to_cross, str(i+1))
#        print('crossings after')
#        [c.info() for c in crossings]

    return crossings, loose_cs, final_cs

def flip_crossings(crossings, loose_cs, final_cs):
    for c in crossings:
        if c != loose_cs.crossing and c != final_cs.crossing:
            c.rotate(randint(0,1))

def cross_strand(crossings, loose_cs, strand_to_cross, new_label):
    opposite = strand_to_cross.opposite()
    new_crossing = Crossing(new_label)
    crossings.append(new_crossing)
    cs0, cs1, cs2, cs3 = new_crossing.crossing_strands()
    connect_crossing_strands(cs0, loose_cs)
    connect_crossing_strands(cs1, strand_to_cross)
    connect_crossing_strands(cs3, opposite)
    return cs2


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
        loose_cs = cross_strand(crossings, loose_cs, strand_to_cross, 'u'+str(i))
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
